import pandas as pd
import os
import re
import gzip
import shutil


# reading the excel file with the information for all mocks
df = pd.read_excel("Data/MOCK_databases_communities_description.xlsx")
print("Total number of commercial mocks used:", df.shape[0])


# filtering and writing mocks on an new excel file
platfomrs = ["ONT", "Illumina", "PacBio"]
writer = pd.ExcelWriter("./Data/MOCK_communities_description_filtered.xlsx", engine='xlsxwriter')
df.to_excel(writer, sheet_name="Comercial mocks", index=False)


for platform in platfomrs:

    mocks = df[df.Platform == platform]
    print(f"Number of mocks sequenced by {platform}:", mocks.shape[0])
    mocks.to_excel(writer, sheet_name=platform, index=False)


writer.close()


# filtering the amplicon sequencing data files by platform
destination_path = "./Data/filtered_mocks"
if "filtered_mocks" not in os.listdir("Data"):
    os.makedirs(destination_path)


for platform in platfomrs:

    df_filtered = pd.read_excel(io="Data/MOCK_communities_description_filtered.xlsx",
         sheet_name=platform)
    mock_names = df_filtered["SRA(NCBI)"].values.tolist()

    if platform not in os.listdir("Data/filtered_mocks"):
        new_directory = os.makedirs(os.path.join(destination_path,f"{platform}"))

    for directory in os.listdir("Data"):

        if re.findall(r'[a-zA-Z\d]+',f"{directory}")[0] in mock_names:
            source_path = os.path.join("Data",directory)
            destination = os.path.join(f"Data/filtered_mocks/{platform}",directory)
            shutil.move(source_path, destination)

    
# extracting the gz files
for platform in platfomrs:
    extract_to_directory = f"Data/filtered_mocks/{platform}"

    for file in os.listdir(f"Data/filtered_mocks/{platform}"):
        # checking if the file has the termination .gz
        if "gz" in re.findall(r'[a-zA-Z\d]+',f"{file}"):
            gz_file_path = os.path.join(f"Data/filtered_mocks/{platform}",file)

            # Open the .gz file and extract its contents
            with gzip.open(gz_file_path, 'rb') as f_in:
                # remove the '.gz' extension
                extracted_file_path = gz_file_path[:-3] 
                
                # Extract the contents of the .gz file
                with open(extracted_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    print(f"{file} successfuly extracted")
