#!/usr/bin/env python3
"""
Script to extract the integer from the name of the last column in a TSV file.

Usage:
    python get_last_column_int.py <tsv_file_path>
    
Example:
    python get_last_column_int.py my_file.tsv
"""

import sys
import pandas as pd
import re

def get_last_column_int(tsv_file_path):
    """
    Read a TSV file and return the integer from the name of the last column.
    
    Args:
        tsv_file_path (str): Path to the TSV file
        
    Returns:
        int: The integer from the last column name
        
    Raises:
        ValueError: If the last column name doesn't contain an integer
        FileNotFoundError: If the file doesn't exist
    """
    # Read the TSV file
    df = pd.read_csv(tsv_file_path, sep='\t')
    
    # Get the last column name
    last_column = df.columns[-1]
    
    # Extract integer from the column name
    # This will find the first integer in the string
    match = re.search(r'\d+', str(last_column))
    
    if match:
        return int(match.group())
    else:
        raise ValueError(f"Last column '{last_column}' does not contain an integer")

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_last_column_int.py <tsv_file_path>")
        print("Example: python get_last_column_int.py my_file.tsv")
        sys.exit(1)
    
    tsv_file_path = sys.argv[1]
    
    try:
        result = get_last_column_int(tsv_file_path)
        print(result)
    except FileNotFoundError:
        print(f"Error: File '{tsv_file_path}' not found", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
