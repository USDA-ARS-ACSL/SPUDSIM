"""
Unit tests for SPUDSIM output functionality
"""

import unittest
import tempfile
import shutil
import sys
from pathlib import Path
import json
import csv

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spudsim.output import OutputManager, OutputConfig
from spudsim.output.config import OutputFormat, OutputFrequency


class TestOutputConfig(unittest.TestCase):
    """Test output configuration"""
    
    def test_default_config(self):
        """Test default configuration values"""
        config = OutputConfig()
        self.assertEqual(config.output_dir, "output")
        self.assertEqual(config.formats, [OutputFormat.CSV])
        self.assertEqual(config.frequency, OutputFrequency.DAILY)
        self.assertTrue(config.compress)
        self.assertEqual(config.decimal_places, 4)
    
    def test_config_to_dict(self):
        """Test configuration serialization"""
        config = OutputConfig(
            output_dir="test_output",
            formats=[OutputFormat.CSV, OutputFormat.JSON]
        )
        config_dict = config.to_dict()
        
        self.assertEqual(config_dict['output_dir'], "test_output")
        self.assertEqual(config_dict['formats'], ['csv', 'json'])
    
    def test_config_from_dict(self):
        """Test configuration deserialization"""
        data = {
            'output_dir': 'test_output',
            'formats': ['csv', 'json'],
            'frequency': 'daily',
            'decimal_places': 3
        }
        config = OutputConfig.from_dict(data)
        
        self.assertEqual(config.output_dir, 'test_output')
        self.assertEqual(len(config.formats), 2)
        self.assertEqual(config.decimal_places, 3)


class TestOutputManager(unittest.TestCase):
    """Test output manager functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_csv_output(self):
        """Test CSV output generation"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.CSV]
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        
        # Write some data
        for i in range(1, 6):
            data = {
                'temperature': 20.0 + i,
                'humidity': 60.0 + i
            }
            output_mgr.write_timestep(i, data)
        
        output_mgr.finalize()
        
        # Verify CSV file was created
        csv_files = list(output_mgr.output_dir.glob("*.csv"))
        self.assertEqual(len(csv_files), 1)
        
        # Verify CSV content
        with open(csv_files[0], 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            self.assertEqual(len(rows), 5)
            self.assertIn('timestep', rows[0])
            self.assertIn('temperature', rows[0])
    
    def test_json_output(self):
        """Test JSON output generation"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.JSON]
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        
        # Write some data
        for i in range(1, 4):
            data = {'value': i * 10}
            output_mgr.write_timestep(i, data)
        
        output_mgr.finalize()
        
        # Verify JSON file was created
        json_files = list(output_mgr.output_dir.glob("*.json"))
        # Should have output file and metadata file
        self.assertGreaterEqual(len(json_files), 1)
        
        # Find output file (not metadata)
        output_files = [f for f in json_files if 'metadata' not in f.name]
        self.assertEqual(len(output_files), 1)
        
        # Verify JSON content
        with open(output_files[0], 'r') as f:
            data = json.load(f)
            self.assertEqual(len(data), 3)
            self.assertEqual(data[0]['timestep'], 1)
    
    def test_variable_filtering(self):
        """Test filtering specific variables"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.CSV],
            variables=['temperature']  # Only output temperature
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        
        # Write data with multiple variables
        data = {
            'temperature': 25.0,
            'humidity': 70.0,
            'pressure': 1013.0
        }
        output_mgr.write_timestep(1, data)
        output_mgr.finalize()
        
        # Verify only selected variable is in output
        csv_files = list(output_mgr.output_dir.glob("*.csv"))
        with open(csv_files[0], 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            self.assertIn('temperature', rows[0])
            self.assertNotIn('humidity', rows[0])
            self.assertNotIn('pressure', rows[0])
    
    def test_decimal_rounding(self):
        """Test decimal place rounding"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.JSON],
            decimal_places=2
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        
        # Write data with high precision
        data = {'value': 3.14159265359}
        output_mgr.write_timestep(1, data)
        output_mgr.finalize()
        
        # Verify rounding
        json_files = [f for f in output_mgr.output_dir.glob("*.json") if 'metadata' not in f.name]
        with open(json_files[0], 'r') as f:
            data = json.load(f)
            self.assertEqual(data[0]['value'], 3.14)
    
    def test_metadata_generation(self):
        """Test metadata file generation"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.CSV],
            include_metadata=True
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        output_mgr.add_metadata('crop', 'potato')
        output_mgr.add_metadata('location', 'Nebraska')
        output_mgr.finalize()
        
        # Verify metadata file exists
        metadata_files = list(output_mgr.output_dir.glob("metadata_*.json"))
        self.assertEqual(len(metadata_files), 1)
        
        # Verify metadata content
        with open(metadata_files[0], 'r') as f:
            metadata = json.load(f)
            self.assertEqual(metadata['crop'], 'potato')
            self.assertEqual(metadata['location'], 'Nebraska')
            self.assertIn('run_id', metadata)
    
    def test_batch_writing(self):
        """Test batch writing functionality"""
        config = OutputConfig(
            output_dir=self.temp_dir,
            formats=[OutputFormat.CSV]
        )
        
        output_mgr = OutputManager(config, run_id="test_run")
        
        # Prepare batch data
        timesteps = [1, 2, 3, 4, 5]
        data_list = [{'value': i * 10} for i in timesteps]
        
        # Write batch
        output_mgr.write_batch(timesteps, data_list)
        output_mgr.finalize()
        
        # Verify all data was written
        csv_files = list(output_mgr.output_dir.glob("*.csv"))
        with open(csv_files[0], 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            self.assertEqual(len(rows), 5)


if __name__ == '__main__':
    unittest.main()
