"""
Example usage of SPUDSIM streamlined output system
"""

import sys
from pathlib import Path
import random

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spudsim.output import OutputManager, OutputConfig
from spudsim.output.config import OutputFormat, OutputFrequency


def simulate_model_run(num_days: int = 120) -> None:
    """
    Simulate a simple potato growth model run
    
    Args:
        num_days: Number of days to simulate
    """
    
    # Configure output
    config = OutputConfig(
        output_dir="output",
        formats=[OutputFormat.CSV, OutputFormat.JSON],
        frequency=OutputFrequency.DAILY,
        compress=True,
        decimal_places=3,
        include_metadata=True
    )
    
    # Create output manager
    output_mgr = OutputManager(config, run_id="example_run")
    
    # Add some metadata
    output_mgr.add_metadata('crop', 'potato')
    output_mgr.add_metadata('variety', 'Russet Burbank')
    output_mgr.add_metadata('location', 'Nebraska')
    output_mgr.add_metadata('planting_date', '2024-05-01')
    
    print(f"Starting model simulation for {num_days} days...")
    print(f"Output directory: {output_mgr.output_dir}")
    
    # Simulate daily timesteps
    for day in range(1, num_days + 1):
        # Simulate some model variables (these would come from actual model)
        data = {
            'day_of_year': day,
            'leaf_area_index': min(5.0, day * 0.05 + random.uniform(-0.1, 0.1)),
            'biomass_kg_ha': day * 50 + random.uniform(-10, 10),
            'soil_water_mm': 150 + random.uniform(-20, 20),
            'canopy_temp_c': 20 + random.uniform(-5, 5),
            'tuber_yield_kg_ha': max(0, (day - 40) * 80 + random.uniform(-50, 50)) if day > 40 else 0,
            'nitrogen_uptake_kg_ha': day * 2 + random.uniform(-0.5, 0.5),
            'evapotranspiration_mm': 5 + random.uniform(-2, 2)
        }
        
        # Write data for this timestep
        output_mgr.write_timestep(day, data)
        
        # Flush every 10 days to avoid large memory buffers
        if day % 10 == 0:
            output_mgr.flush()
            print(f"  Day {day}: LAI={data['leaf_area_index']:.2f}, Yield={data['tuber_yield_kg_ha']:.1f} kg/ha")
    
    # Finalize output
    output_mgr.finalize()
    
    print(f"\nSimulation complete!")
    print(f"Output files written to: {output_mgr.output_dir}")
    print("\nOutput files created:")
    for file in sorted(output_mgr.output_dir.iterdir()):
        size_kb = file.stat().st_size / 1024
        print(f"  - {file.name} ({size_kb:.1f} KB)")


def demonstrate_batch_output():
    """Demonstrate efficient batch output writing"""
    
    print("\n" + "="*60)
    print("Demonstrating batch output (more efficient)")
    print("="*60 + "\n")
    
    config = OutputConfig(
        output_dir="output",
        formats=[OutputFormat.CSV],
        frequency=OutputFrequency.DAILY,
        variables=['biomass_kg_ha', 'tuber_yield_kg_ha']  # Only output specific variables
    )
    
    output_mgr = OutputManager(config, run_id="batch_example")
    
    # Prepare batch data
    timesteps = list(range(1, 121))
    data_list = []
    
    for day in timesteps:
        data = {
            'biomass_kg_ha': day * 50 + random.uniform(-10, 10),
            'tuber_yield_kg_ha': max(0, (day - 40) * 80) if day > 40 else 0,
            'soil_water_mm': 150,  # This will be filtered out
        }
        data_list.append(data)
    
    # Write all data at once (more efficient)
    output_mgr.write_batch(timesteps, data_list)
    output_mgr.finalize()
    
    print(f"Batch output complete: {output_mgr.output_dir}")


def demonstrate_multiple_formats():
    """Demonstrate output to multiple formats simultaneously"""
    
    print("\n" + "="*60)
    print("Demonstrating multiple output formats")
    print("="*60 + "\n")
    
    config = OutputConfig(
        output_dir="output",
        formats=[OutputFormat.CSV, OutputFormat.JSON, OutputFormat.HDF5],
        frequency=OutputFrequency.DAILY,
        compress=True
    )
    
    output_mgr = OutputManager(config, run_id="multi_format_example")
    
    # Simulate some data
    for day in range(1, 31):
        data = {
            'temperature_c': 20 + random.uniform(-5, 5),
            'humidity_pct': 60 + random.uniform(-10, 10),
            'radiation_mj_m2': 15 + random.uniform(-3, 3)
        }
        output_mgr.write_timestep(day, data)
    
    output_mgr.finalize()
    
    print(f"Multi-format output complete: {output_mgr.output_dir}")
    print("\nFormats created:")
    for file in sorted(output_mgr.output_dir.iterdir()):
        print(f"  - {file.name}")


if __name__ == "__main__":
    # Run main example
    simulate_model_run(num_days=120)
    
    # Run additional demonstrations
    demonstrate_batch_output()
    demonstrate_multiple_formats()
    
    print("\n" + "="*60)
    print("All examples completed successfully!")
    print("="*60)
