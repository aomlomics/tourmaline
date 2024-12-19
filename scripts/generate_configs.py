#!/usr/bin/env python3

import yaml
import itertools
import os
import argparse
from pathlib import Path

def load_yaml(file_path):
    """Load a YAML file."""
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def generate_parameter_combinations(param_space):
    """Generate all combinations of parameters from the parameter space."""
    keys = param_space.keys()
    values = param_space.values()
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

def create_config_file(base_config, params, output_dir, index, run_name_prefix, config_prefix):
    """Create a new config file with the given parameters."""
    # Create a copy of the base config and update with new parameters
    new_config = base_config.copy()
    new_config.update(params)
    
    # Update run_name with prefix and index
    new_config['run_name'] = f"{run_name_prefix}_{index:03d}"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate output filename using the provided config prefix
    output_path = os.path.join(output_dir, f'{config_prefix}_{index:03d}.yaml')
    
    # Write the new config file
    with open(output_path, 'w') as f:
        yaml.dump(new_config, f, default_flow_style=False)
    
    return output_path

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate config files from parameter space')
    parser.add_argument('base_config', help='Path to base config file')
    parser.add_argument('param_space_config', help='Path to parameter space config file')
    args = parser.parse_args()
    
    # Load the parameter space config
    param_space_config = load_yaml(args.param_space_config)
    param_space = param_space_config['parameter_space']
    output_dir = param_space_config.get('output_dir', 'parameter_sweep_configs')
    run_name_prefix = param_space_config.get('run_name_prefix', 'test_data')
    config_prefix = param_space_config.get('config_prefix', 'config-01-qaqc')
    
    # Load the base config
    base_config = load_yaml(args.base_config)
    
    # Generate all parameter combinations
    combinations = generate_parameter_combinations(param_space)
    
    # Generate config files
    generated_files = []
    for i, params in enumerate(combinations):
        output_path = create_config_file(base_config, params, output_dir, i, 
                                       run_name_prefix, config_prefix)
        generated_files.append(output_path)
        
        # Print parameter combination
        print(f"\nGenerated config file {output_path} with parameters:")
        print(f"  run_name: {run_name_prefix}_{i:03d}")
        for key, value in params.items():
            print(f"  {key}: {value}")
    
    print(f"\nTotal config files generated: {len(generated_files)}")

if __name__ == "__main__":
    main()
