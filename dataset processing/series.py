import sys
from collections import Counter

def process_file(filename, n):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Find the first and second occurrences of lines starting with a digit
    start_idx = -1
    end_idx = -1
    for i, line in enumerate(lines):
        if line[0].isdigit():
            if start_idx == -1:
                start_idx = i
            else:
                end_idx = i
                break

    if start_idx == -1 or end_idx == -1:
        print("Could not find the required lines starting with a digit.")
        return

    # Extract the relevant lines, excluding the first line and the last n lines
    relevant_lines = lines[start_idx+2:end_idx-n]

    # Count occurrences of content before the first space in each relevant line
    counter = Counter()
    for line in relevant_lines:
        first_space_idx = line.find(' ')
        if first_space_idx != -1:
            content_before_space = line[:first_space_idx]
            counter[content_before_space] += 1

    return counter

if __name__ == "__main__":
    # Example usage: python script.py filename.extxyz 3
    if len(sys.argv) != 3:
        print("Usage: python script.py <filename.extxyz> <n>")
        sys.exit(1)

    filename = sys.argv[1]
    n = int(sys.argv[2])
    result = process_file(filename, n)

    # Print the result
    if result:
        for item, count in result.items():
            print(f"{item}{count}")
