import pandas as pd
import shutil
import os

# Set the paths
xlsx_path = r'E:\oc20\oc\all\surface_cp\u20\111.xlsx'  # Replace with the path to your Excel file
source_folder = r'E:\oc20\oc\all\surface_cp'  # Replace with the source folder path containing the files to be copied
destination_folder = r'E:\oc20\oc\all\surface_cp\u20'  # Replace with the destination folder path

# Read the first column of the Excel file
df = pd.read_excel(xlsx_path)
file_names = df.iloc[:, 0].astype(str).tolist()  # Convert the first column to a list of strings

# Copy files with the same name
for file_name in file_names:
    source_file = os.path.join(source_folder, file_name)
    destination_file = os.path.join(destination_folder, file_name)
    if os.path.isfile(source_file):  # Check if the source file exists
        shutil.copy(source_file, destination_file)
        #shutil.move(source_file, destination_file)