"""
Configuration management for streamlined data output
"""

from dataclasses import dataclass, field
from typing import List, Optional
from enum import Enum


class OutputFormat(Enum):
    """Supported output formats"""
    CSV = "csv"
    JSON = "json"
    HDF5 = "hdf5"
    NETCDF = "netcdf"


class OutputFrequency(Enum):
    """Output frequency options"""
    DAILY = "daily"
    HOURLY = "hourly"
    FINAL = "final"


@dataclass
class OutputConfig:
    """Configuration for model output"""
    
    # Output directory
    output_dir: str = "output"
    
    # Output formats to use
    formats: List[OutputFormat] = field(default_factory=lambda: [OutputFormat.CSV])
    
    # Output frequency
    frequency: OutputFrequency = OutputFrequency.DAILY
    
    # Variables to output (None = all variables)
    variables: Optional[List[str]] = None
    
    # Compression for supported formats
    compress: bool = True
    
    # Precision for floating point values
    decimal_places: int = 4
    
    # Include metadata in output
    include_metadata: bool = True
    
    # Separate files per variable vs combined
    separate_variables: bool = False
    
    # Append to existing files or overwrite
    append_mode: bool = False
    
    def to_dict(self) -> dict:
        """Convert configuration to dictionary"""
        return {
            'output_dir': self.output_dir,
            'formats': [f.value for f in self.formats],
            'frequency': self.frequency.value,
            'variables': self.variables,
            'compress': self.compress,
            'decimal_places': self.decimal_places,
            'include_metadata': self.include_metadata,
            'separate_variables': self.separate_variables,
            'append_mode': self.append_mode
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'OutputConfig':
        """Create configuration from dictionary"""
        config = cls()
        config.output_dir = data.get('output_dir', config.output_dir)
        config.formats = [OutputFormat(f) for f in data.get('formats', ['csv'])]
        config.frequency = OutputFrequency(data.get('frequency', 'daily'))
        config.variables = data.get('variables')
        config.compress = data.get('compress', True)
        config.decimal_places = data.get('decimal_places', 4)
        config.include_metadata = data.get('include_metadata', True)
        config.separate_variables = data.get('separate_variables', False)
        config.append_mode = data.get('append_mode', False)
        return config
