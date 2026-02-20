# SPUDSIM

USDA-ARS ACSL potato simulation model with streamlined data output.

## Features

### Streamlined Data Output

SPUDSIM includes a modern, efficient data output system:

- **Multiple formats**: CSV, JSON, HDF5, NetCDF
- **Configurable**: Choose variables, frequency, and formats
- **Efficient**: Buffered I/O and compression support
- **Metadata**: Automatic tracking of run configuration

## Installation

### Requirements

- Python 3.7+
- Optional: h5py (for HDF5), netCDF4 (for NetCDF)

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/USDA-ARS-ACSL/SPUDSIM.git
cd SPUDSIM

# Install in development mode
pip install -e .

# Optional: Install additional format support
pip install -e ".[all]"  # Install all optional dependencies
# Or install specific formats:
pip install h5py  # For HDF5 support
pip install netCDF4  # For NetCDF support
```

## Quick Start

### Running an Example

```bash
python examples/example_usage.py
```

This will:
1. Run a simulated potato growth model for 120 days
2. Generate output in CSV and JSON formats
3. Create metadata files with run information

### Using in Your Code

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
output_mgr = OutputManager(config, run_id="my_simulation")

# Run your model and output data
for timestep in range(1, 121):
    # Your model calculations here
    data = {
        'leaf_area_index': 3.5,
        'biomass_kg_ha': 5000,
        'tuber_yield_kg_ha': 12000
    }
    output_mgr.write_timestep(timestep, data)

# Finalize
output_mgr.finalize()
```

## Documentation

- [Data Output Guide](docs/OUTPUT.md) - Comprehensive guide to the output system
- See `examples/` directory for usage examples
- Run tests: `python -m pytest tests/`

## Project Structure

```
SPUDSIM/
├── spudsim/              # Main package
│   ├── __init__.py
│   └── output/           # Data output module
│       ├── __init__.py
│       ├── config.py     # Configuration classes
│       └── output_manager.py  # Output manager
├── tests/                # Unit tests
│   └── test_output.py
├── examples/             # Usage examples
│   └── example_usage.py
├── docs/                 # Documentation
│   └── OUTPUT.md
└── README.md
```

## Testing

Run the test suite:

```bash
python -m pytest tests/ -v
```

Or run tests individually:

```bash
python tests/test_output.py
```

## Examples

The `examples/` directory contains:

- `example_usage.py` - Comprehensive examples of the output system

Run examples:

```bash
cd examples
python example_usage.py
```

## Output Formats

### CSV
- Human-readable, widely compatible
- Good for small to medium datasets
- Easy to open in Excel or other tools

### JSON  
- Structured data, easy to parse
- Good for web applications and APIs
- Human-readable with formatting

### HDF5
- Efficient for large datasets
- Supports compression and chunking
- Industry standard for scientific computing

### NetCDF
- Self-describing format
- CF-compliant for climate/weather data
- Standard in earth sciences

## Configuration Options

The output system is highly configurable:

```python
config = OutputConfig(
    output_dir="output",           # Output directory
    formats=[OutputFormat.CSV],    # Output formats
    frequency=OutputFrequency.DAILY,  # Output frequency
    variables=['var1', 'var2'],    # Variables to output
    compress=True,                  # Enable compression
    decimal_places=4,               # Floating point precision
    include_metadata=True           # Generate metadata files
)
```

See [OUTPUT.md](docs/OUTPUT.md) for complete documentation.

## Contributing

Contact Dr. David Fleisher (david.fleisher@usda.gov) for questions or contributions.

## License

This software is in the public domain as a work of the United States Government.

## Citation

When using SPUDSIM in your research, please cite:

USDA-ARS Adaptive Cropping Systems Laboratory. (2024). SPUDSIM: Potato Simulation Model. 
Retrieved from https://github.com/USDA-ARS-ACSL/SPUDSIM

## Contact

For assistance with SPUDSIM:
- Dr. David Fleisher - david.fleisher@usda.gov
- USDA-ARS Adaptive Cropping Systems Laboratory

## Version

Current version: 2.2.0
