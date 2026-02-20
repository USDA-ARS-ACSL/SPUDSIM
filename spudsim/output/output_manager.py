"""
Streamlined output manager for SPUDSIM model runs
"""

import os
import csv
import json
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime, timezone
import logging

from .config import OutputConfig, OutputFormat

logger = logging.getLogger(__name__)


class OutputManager:
    """
    Manages streamlined data output for SPUDSIM model runs
    
    Features:
    - Multiple output formats (CSV, JSON, HDF5, NetCDF)
    - Configurable output frequency and variables
    - Efficient buffering for large datasets
    - Automatic compression for supported formats
    """
    
    def __init__(self, config: OutputConfig, run_id: Optional[str] = None):
        """
        Initialize output manager
        
        Args:
            config: Output configuration
            run_id: Unique identifier for this model run
        """
        self.config = config
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(config.output_dir) / self.run_id
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize buffers for each format
        self._buffers: Dict[OutputFormat, List] = {fmt: [] for fmt in config.formats}
        
        # Track metadata
        self.metadata: Dict[str, Any] = {
            'run_id': self.run_id,
            'start_time': datetime.now(timezone.utc).isoformat(),
            'config': config.to_dict()
        }
        
        logger.info(f"OutputManager initialized for run {self.run_id}")
    
    def write_timestep(self, timestep: int, data: Dict[str, Any]) -> None:
        """
        Write data for a single timestep
        
        Args:
            timestep: Timestep number
            data: Dictionary of variable names to values
        """
        # Filter variables if specified
        if self.config.variables:
            data = {k: v for k, v in data.items() if k in self.config.variables}
        
        # Round floating point values
        data = self._round_data(data)
        
        # Add timestep to data
        record = {'timestep': timestep, **data}
        
        # Buffer the data
        for fmt in self.config.formats:
            self._buffers[fmt].append(record)
    
    def write_batch(self, timesteps: List[int], data_list: List[Dict[str, Any]]) -> None:
        """
        Write multiple timesteps at once (more efficient)
        
        Args:
            timesteps: List of timestep numbers
            data_list: List of data dictionaries
        """
        for timestep, data in zip(timesteps, data_list):
            self.write_timestep(timestep, data)
    
    def flush(self) -> None:
        """Flush all buffered data to disk"""
        for fmt in self.config.formats:
            if self._buffers[fmt]:
                self._write_format(fmt)
                self._buffers[fmt].clear()
        
        logger.info(f"Flushed all output buffers for run {self.run_id}")
    
    def finalize(self) -> None:
        """Finalize output and write metadata"""
        # Flush remaining data
        self.flush()
        
        # Write metadata
        if self.config.include_metadata:
            self.metadata['end_time'] = datetime.now(timezone.utc).isoformat()
            self._write_metadata()
        
        logger.info(f"Finalized output for run {self.run_id}")
    
    def _write_format(self, fmt: OutputFormat) -> None:
        """Write buffered data in specific format"""
        if fmt == OutputFormat.CSV:
            self._write_csv()
        elif fmt == OutputFormat.JSON:
            self._write_json()
        elif fmt == OutputFormat.HDF5:
            self._write_hdf5()
        elif fmt == OutputFormat.NETCDF:
            self._write_netcdf()
    
    def _write_csv(self) -> None:
        """Write data to CSV format"""
        data = self._buffers[OutputFormat.CSV]
        if not data:
            return
        
        filepath = self.output_dir / f"output_{self.run_id}.csv"
        mode = 'a' if self.config.append_mode and filepath.exists() else 'w'
        
        # Get all fieldnames from data
        fieldnames = ['timestep'] + [k for k in data[0].keys() if k != 'timestep']
        
        with open(filepath, mode, newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            
            # Write header only if new file or not appending
            if mode == 'w':
                writer.writeheader()
            
            writer.writerows(data)
        
        logger.debug(f"Wrote {len(data)} records to {filepath}")
    
    def _write_json(self) -> None:
        """Write data to JSON format"""
        data = self._buffers[OutputFormat.JSON]
        if not data:
            return
        
        filepath = self.output_dir / f"output_{self.run_id}.json"
        
        if self.config.append_mode and filepath.exists():
            # Read existing data and append
            with open(filepath, 'r') as f:
                existing = json.load(f)
            data = existing + data
        
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2 if not self.config.compress else None)
        
        logger.debug(f"Wrote {len(data)} records to {filepath}")
    
    def _write_hdf5(self) -> None:
        """Write data to HDF5 format"""
        try:
            import h5py
            import numpy as np
        except ImportError:
            logger.warning("h5py not installed, skipping HDF5 output")
            return
        
        data = self._buffers[OutputFormat.HDF5]
        if not data:
            return
        
        filepath = self.output_dir / f"output_{self.run_id}.h5"
        
        mode = 'a' if self.config.append_mode and filepath.exists() else 'w'
        
        with h5py.File(filepath, mode) as f:
            # Convert data to structured array
            if data:
                keys = list(data[0].keys())
                
                # Create datasets for each variable
                for key in keys:
                    values = [record[key] for record in data]
                    
                    if key not in f:
                        # Create new dataset
                        maxshape = (None,) if self.config.append_mode else None
                        f.create_dataset(
                            key, 
                            data=values,
                            compression='gzip' if self.config.compress else None,
                            maxshape=maxshape
                        )
                    else:
                        # Append to existing dataset
                        dset = f[key]
                        old_size = dset.shape[0]
                        dset.resize(old_size + len(values), axis=0)
                        dset[old_size:] = values
        
        logger.debug(f"Wrote {len(data)} records to {filepath}")
    
    def _write_netcdf(self) -> None:
        """Write data to NetCDF format"""
        try:
            import netCDF4 as nc
            import numpy as np
        except ImportError:
            logger.warning("netCDF4 not installed, skipping NetCDF output")
            return
        
        data = self._buffers[OutputFormat.NETCDF]
        if not data:
            return
        
        filepath = self.output_dir / f"output_{self.run_id}.nc"
        
        mode = 'a' if self.config.append_mode and filepath.exists() else 'w'
        
        with nc.Dataset(filepath, mode) as ds:
            if mode == 'w':
                # Create dimensions
                ds.createDimension('time', None)  # Unlimited dimension
            
            time_size = len(ds.dimensions['time']) if mode == 'a' else 0
            
            # Add data for each variable
            if data:
                keys = list(data[0].keys())
                
                for key in keys:
                    values = [record[key] for record in data]
                    
                    if key not in ds.variables:
                        # Create new variable
                        var = ds.createVariable(
                            key, 
                            'f8',  # Double precision
                            ('time',),
                            compression='zlib' if self.config.compress else None
                        )
                        var[:] = values
                    else:
                        # Append to existing variable
                        var = ds.variables[key]
                        var[time_size:] = values
        
        logger.debug(f"Wrote {len(data)} records to {filepath}")
    
    def _write_metadata(self) -> None:
        """Write metadata to JSON file"""
        filepath = self.output_dir / f"metadata_{self.run_id}.json"
        with open(filepath, 'w') as f:
            json.dump(self.metadata, f, indent=2)
        logger.debug(f"Wrote metadata to {filepath}")
    
    def _round_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Round floating point values to configured precision"""
        rounded = {}
        for key, value in data.items():
            if isinstance(value, float):
                rounded[key] = round(value, self.config.decimal_places)
            else:
                rounded[key] = value
        return rounded
    
    def add_metadata(self, key: str, value: Any) -> None:
        """Add custom metadata"""
        self.metadata[key] = value
