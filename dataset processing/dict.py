import os
import sys
import shutil
from collections import OrderedDict

def extract_columns(filename):
    """
    Reads a file and extracts the first column of lines where the sixth column is '2'.
    It starts searching from the end of the file until it finds the first line that starts with a digit.

    Parameters:
    filename (str): The name of the file to read.

    Returns:
    dict: A dictionary with counts of the first column values.
    """
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        num_line_index = None

        # Find the first line starting with a digit
        for i in range(len(lines) - 1, -1, -1):
            if lines[i][0].isdigit():
                num_line_index = i
                break
        
        if num_line_index is None:
            return {}
        
        # Count occurrences of the first column for lines with sixth column as '2'
        column1_count = {}
        for line in lines[num_line_index:]:
            columns = line.split()
            if len(columns) >= 4 and columns[5] == '2':
                column1_value = columns[0]
                column1_count[column1_value] = column1_count.get(column1_value, 0) + 1
        
        return column1_count
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return {}
    except Exception as e:
        print(f"Error: An exception occurred - {str(e)}")
        return {}
'''        
def check_output(result, filename):
    
    if 'N' in result and 'H' in result and 'O' in result and 'C' not in result:
        target_folder = f"n{result['N']}h{result['H']}o{result['O']}"
    elif 'N' in result and 'H' in result and 'C' not in result:
        target_folder = f"n{result['N']}h{result['H']}"
    elif 'N' in result and 'O' in result and 'C' not in result:
        target_folder = f"n{result['N']}o{result['O']}"
    elif 'N' in result and 'C' not in result:
        target_folder = f"n{result['N']}"
    else:
        return
    
    os.makedirs(target_folder, exist_ok=True)
    shutil.copy(filename, target_folder)
    print(f"File '{filename}' copied to '{target_folder}'.")

# Get command line arguments
if len(sys.argv) < 2:
    print("Usage: python3 111.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
result = extract_columns(filename)
check_output(result, filename)

"""

def check_output(result, filename):
    
    if 'N' in result and 'H' in result and 'O' in result and 'C' not in result:
        target_folder = f"n{result['N']}h{result['H']}o{result['O']}"
    elif 'N' in result and 'H' in result and 'C' not in result:
        target_folder = f"n{result['N']}h{result['H']}"
    elif 'N' in result and 'O' in result and 'C' not in result:
        target_folder = f"n{result['N']}o{result['O']}"
    elif 'N' in result and 'C' not in result:
        target_folder = f"n{result['N']}"
    else:
        return
    
    os.makedirs(target_folder, exist_ok=True)
    shutil.copy(filename, target_folder)
    print(f"File '{filename}' copied to '{target_folder}'.")
'''    
    
def check_output(result, filename):
    if result.get('H') == 2 and result.get('O') == 1 and 'C' not in result and 'N' not in result:
        target_folder = "h2o"
    elif result.get('H') == 1 and 'O' not in result and 'C' not in result and 'N' not in result:
        target_folder = "h"
    else:
        return

    os.makedirs(target_folder, exist_ok=True)
    if os.path.exists(filename):
        shutil.copy(filename, target_folder)
        print(f"File '{filename}' copied to '{target_folder}'.")
    else:
        print(f"Error: File '{filename}' does not exist.")
        
   
