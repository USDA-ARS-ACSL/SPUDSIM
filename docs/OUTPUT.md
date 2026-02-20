# SPUDSIM - Streamlined Data Output

This document describes the streamlined data output system for SPUDSIM model runs.

## Overview

The SPUDSIM data output system provides:

- **Multiple output formats**: CSV, JSON, HDF5, NetCDF
- **Configurable output frequency**: Daily, hourly, or final only
- **Variable selection**: Choose which variables to output
- **Efficient buffering**: Minimize I/O operations for large simulations
- **Automatic compression**: Reduce file sizes for supported formats
- **Metadata tracking**: Store run configuration and metadata

## Quick Start

### Basic Usage

```python
from spudsim.output import OutputManager, OutputConfig
from spudsim.output.config import OutputFormat

# Configure output
config = OutputConfig(
    output_dir="output",
    formats=[OutputFormat.CSV, OutputFormat.JSON],
    compress=True
)

# Create output manager
output_mgr = OutputManager(config, run_id="my_run")

# Write timestep data
for day in range(1, 121):
    data = {
        'leaf_area_index': 3.5,
        'biomass_kg_ha': 5000,
        'soil_water_mm': 150
    }
    output_mgr.write_timestep(day, data)

# Finalize and write all data
output_mgr.finalize()
```

## Configuration Options

### OutputConfig Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_dir` | str | "output" | Directory for output files |
| `formats` | List[OutputFormat] | [CSV] | Output formats to generate |
| `frequency` | OutputFrequency | DAILY | How often to output data |
| `variables` | List[str] | None | Variables to output (None = all) |
| `compress` | bool | True | Enable compression for supported formats |
| `decimal_places` | int | 4 | Precision for floating point values |
| `include_metadata` | bool | True | Generate metadata file |
| `separate_variables` | bool | False | Create separate files per variable |
| `append_mode` | bool | False | Append to existing files |

### Output Formats

- **CSV**: Human-readable, widely compatible
- **JSON**: Structured data, easy to parse
- **HDF5**: Efficient for large datasets, supports compression
- **NetCDF**: Standard for scientific data, self-describing

## Best Practices

### 1. Use Batch Writing for Large Datasets

```python
# More efficient
output_mgr.write_batch(list(range(1, 1000)), data_list)
```

### 2. Flush Periodically

```python
for day in range(1, 365):
    output_mgr.write_timestep(day, data)
    if day % 30 == 0:  # Flush monthly
        output_mgr.flush()
```

### 3. Filter Variables

```python
config = OutputConfig(
    variables=['biomass_kg_ha', 'tuber_yield_kg_ha']
)
```

### 4. Choose Appropriate Formats

- **CSV**: For small datasets, data sharing, spreadsheet analysis
- **JSON**: For web applications, APIs
- **HDF5**: For large datasets, numerical analysis, compression
- **NetCDF**: For scientific data sharing, CF compliance

## Output File Structure

```
output/
├── 20240520_143022/
│   ├── output_20240520_143022.csv
│   ├── output_20240520_143022.json
│   └── metadata_20240520_143022.json
```

## API Reference

See inline documentation in:
- `spudsim/output/config.py` - Configuration classes
- `spudsim/output/output_manager.py` - Output manager implementation
