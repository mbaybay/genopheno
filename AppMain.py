__author__ = 'Lakshmi Arbatti'

import os
import pandas as pd
from GeneExtractionClass import GeneExtraction
import numpy as np

# Local path and data directory
DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Data")
# The folder which stores the DBSNp queries zip files
db_snps_path = "SNP_results"
# This is the consolidated database built with every snp and its corresponding gene info available
db_snps_output = "snp_database.csv"
# This is the folder in which you have the consolidated genes-snps file
main_gene_file = os.path.join(DATA_PATH,db_snps_path,db_snps_output)
# This is the column in the snp_database.csv file that needs to be split for multiple genes
col_split = "Gene_info"

# The list of humans (user id) that have available defines eye color phenotype
pheno_input_file_name = "phenolt2300.csv"
# The file that contains the binned info of users and eye color
pheno_output_file_name = "final_eyecolor_phenotypes.csv"
# This is the folder which contains all the user data from 23&me and Ancestory.com
# user_folder = "user_subset"
user_folder = "user"
user_path = os.path.join(DATA_PATH,user_folder)
output_inter = os.path.join(DATA_PATH,"output_intermediate")
# out_path_final = os.path.join(DATA_PATH,"output_subset")
out_path = os.path.join(DATA_PATH,"output_temp")
# out_path = os.path.join(DATA_PATH,"output_temp_subset")
output_final = os.path.join(DATA_PATH,"output_final")

gene_extractor = GeneExtraction()

# This function reads each snp per row and converts to 0 or 1 based of if there has been a mutation or not
def snp_count(x):
    if (x in (1,2)):
        ret = 1
    else:
        ret = 0
    return ret

def genotype_phenotype_extraction(usr_file_lst, phenotype, all_gene_data):

    # Get the phenotype information to append to the user name in final.csv
    pheno = "_" + phenotype
    merged = all_gene_data.drop(['Ref', 'Alt'], axis=1)
    print "Total number of users " , len(usr_file_lst)
    len_users = 0
    counter = 0

    for user in usr_file_lst:
        counter+=1
        print "Counter : ", counter, len(usr_file_lst)
        # Dividing into 20s to assist memory issues
        if counter % 20 != 0 and counter <= len(usr_file_lst):
            write = False
            # Extract just the user id from the file name
            short_name = user.split('_y')[0] + '_'
            # Added this loop to take care of other file types other than 23&me and ancestory.ext
            # if "23andme.txt" in user or "ancestry.txt" in user:
            print "Processing file : ", user
            # If present, read the file and build allele transformation info
            full_file_name = os.path.join(user_path,user)
            data = gene_extractor.build_user_data(full_file_name,all_gene_data)

            # Take care of incorrect file formats
            if not data.empty:
                # Number of actual users
                len_users +=1
                data = gene_extractor.allele_transformation(short_name, pheno, data)
                # Merge with the earlier dataset to build the data corpus
                merged = pd.merge(merged,data,on=['Rsid','Gene_info'], how='outer' )
                # Get the snp count for each user
                merged[short_name+'gene'+ pheno] = (merged [short_name + 'snp' + pheno].apply(snp_count)).astype(np.uint8)
                merged[short_name + 'snp' + pheno].fillna(-1,inplace=True)
                merged[short_name+'snp'+ pheno] = merged[short_name+'snp'+ pheno].astype(np.uint8)
                print "Size of merged : ", merged.size
                print "Files merged : ", len_users

        else:
            print "In outside while loop"
            # Set Rsid as index - This is only to avoid Pandas from adding a row number. Increases the size of the output
            merged = merged.set_index(['Rsid'])
            # if its a multiple of 20, write the file into a temporary location
            print "File name being written : " , os.path.join(out_path, "temp_file" + pheno + "_"+str(len_users) + ".csv")
            merged.to_csv(os.path.join(out_path, "temp_file" + pheno + "_"+str(len_users) + ".csv"))
            write = True
            # Reset the merge object and start afresh
            merged = all_gene_data.drop(['Ref', 'Alt'], axis=1)
    # To take care of corner cases, in case the file count is less than 20 and there are corrupt files in the list
    if(not write):
        # Set Rsid as index - This is only to avoid Pandas from adding a row number. Increases the size of the output
        merged = merged.set_index(['Rsid'])
        # if its a multiple of 20, write the file into a temporary location
        print "File name being written : ", os.path.join(out_path, "temp_file" + pheno + "_" + str(len_users) + ".csv")
        merged.to_csv(os.path.join(out_path, "temp_file" + pheno + "_" + str(len_users) + ".csv"))

    print "Final number of users with valid data for phenotype  : ", phenotype , "  is : ",  len_users
    return len_users
# This function reads all the sub-output files and then merges them into a single output file
def merge_files(pheno,len_users):
    print "Merging all the files..."
    lst_merged_users = os.listdir(out_path)
    lst_users = [file_name for file_name in lst_merged_users if pheno in file_name]
    merged = pd.DataFrame()
    for count,user in enumerate(lst_users):
        file_name = os.path.join(out_path,user)
        print "File being merged\n" , file_name
        print "Merging %d of %d" %(count+1,len(lst_users))
        # First time reading a file
        if merged.empty:
            merged = pd.read_csv(file_name,na_values=0)
        else:
            temp = pd.read_csv(file_name,na_values=0)
            merged = pd.merge(merged,temp,on=['Rsid','Gene_info'], how='outer')
        merged.fillna(255,inplace=True)
        df_col_lst = list(merged)
        df_col_lst = df_col_lst[2:]
        merged[df_col_lst] = merged[df_col_lst].astype(np.uint8)
        # print merged.info()
    merged = merged.set_index(['Rsid'])
    print "Writing file : ", os.path.join(output_inter, "Inter_file" + pheno + ".csv")
    merged.to_csv(os.path.join(output_inter, "Inter_file" + pheno + ".csv"))

# Process the merged file
def write_files(suffix,len_users):
    print "Writing final output files..."
    merged_group = pd.DataFrame()
    # Set the chunk size to about 10KMB and read file one chunk at a time
    chunksize = 10 ** 3
    # For the first write
    header = True
    # Breaking down into chunks of 1000MB to support scalability
    for merged in pd.read_csv(os.path.join(output_inter, "Inter_file" + suffix + ".csv"), chunksize=chunksize):
        # Get the list of columns
        df_col_lst = list(merged)
        # Starting from col 3, jump every two columns (this will do the mutated snp count)
        jump_list = df_col_lst[3::2]
        print jump_list

        # Sum the snp mutation count
        merged[jump_list] = (merged[jump_list].replace(255,0))
        merged[jump_list] = merged[jump_list].astype(np.uint8)
        merged['total_user_per_snp'] = merged[jump_list].sum(axis=1)
        merged['total_user_per_snp'] = (merged['total_user_per_snp']).astype(np.uint8)

        # Sets the gene count to 0 or 1 depending on the allele transformation value
        merged['gene_mutation_bool'] = [1 if x != 0  else 0 for x in merged['total_user_per_snp']]
        merged['gene_mutation_bool'] = (merged['gene_mutation_bool']).astype(np.uint8)

        # Delete rows that are all 0s, i.e rows that have no mutation at all - Should we do this or keep the data
        # as negative test cases
        merged = merged[merged['gene_mutation_bool'] != 0]
        merged = merged.replace(255,'NaN')
        # Rearrange columns so that totals are arranged towards the beginning of the table.
        cols = merged.columns.tolist()
        cols = cols[1:2] + cols[0:1] + cols[-1:] + cols[2:-1]
        merged = merged[cols]

        # Set Rsid as index - This is only to avoid Pandas from adding a row number. Increases the size of the output
        merged = merged.set_index(['Rsid'])
        merged_group = pd.concat([merged_group, merged])
        # Save to csv
        try:
            print "Writing file - In Write Files Line 185 : ", os.path.join(output_final, "final_file" + suffix + ".csv")
            # print merged.info()
            with open(os.path.join(output_final, "final_file" + suffix + ".csv"),'a') as final_out:
                merged.to_csv(final_out,header=header)
                header = False

        except:
            print "Error writing final output file. Please check your files again..."

    merged_group = merged_group.reset_index()
    # Drop the two rows and then recalculate and add them back
    merged_group = merged_group.drop(['gene_mutation_bool' , 'total_user_per_snp'] , axis = 1)
    # Aggregate by count of mutations per gene
    merged_group = merged_group.groupby('Gene_info').sum()
    # Count the number of users that have the gene mutated
    merged_group['total_user_per_gene'+ suffix] = (merged_group.astype(bool).sum(axis=1)).astype(np.uint8)
    # This provides the percentage of users per phenotype that had mutation per gene
    merged_group['percent_users' + suffix] = (merged_group['total_user_per_gene' + suffix]/len_users) * 100
    # Rearrange the cols so that percentage figures appear towards the beginning of the table
    cols = merged_group.columns.tolist()
    cols = cols[-1:] + cols[0:-1]
    merged_group = merged_group[cols]
    merged_group.sort_values(by = ['percent_users' + suffix], ascending = [False], inplace = True)
    # Save to csv
    try:
        print "Writing file : ", os.path.join(output_final, "final_grouped" + suffix + ".csv")
        merged_group.to_csv()
    except:
        print "Error writing final grouped file. Please check your files again..."

    return os.path.join(output_final, "final_file" + suffix + ".csv") ,os.path.join(output_final, "final_grouped" + suffix + ".csv")

def file_analysis(file_path,suffix):
    data = pd.read_csv(file_path)
    # row_name = 'percent_users' + suffix
    # print row_name
    # data = data[data[row_name] > 60]
    # data.to_csv("C:\Lakshmi\MSHI\Github\Genotype-Phenotype-Project\GeneExtraction\Data\output_final\\final_grouped_Blue_Green_subset.csv")
    print data.loc[data["Gene_info"] == "4948:OCA2"]

def main():

    # Extract only 23&me and ancestory.com files
    user_lst = os.listdir(user_path)
    user_file_lst = list()
    for file in user_lst:
        if file.endswith("23andme.txt") or file.endswith("ancestry.txt"):
            user_file_lst.append(file)

    gene_extractor = GeneExtraction()

    # 01. Use Shraddha's 'batchquery_processing' logic
    print "Processing batch query..."
    # gene_extractor.batchquery_processing(os.path.join(DATA_PATH,db_snps_path),main_gene_file)

    # 02. Read all snp-gene data and build a data table
    print "Reading all snps and gene data..."
    all_gene_data = gene_extractor.build_gene_data(DATA_PATH,main_gene_file,col_split,sep="|")
     # Drop unnecessary duplicates
    all_gene_data = all_gene_data.drop_duplicates(subset=['Rsid','Gene_info'])
    # Drop rows where there is no gene info
    all_gene_data = all_gene_data[all_gene_data['Gene_info']!='nan']
    # Save for future retrieval
    all_gene_data.to_csv(os.path.join(DATA_PATH,"all_gene.csv"))

    # 03. For each user, read the phenotype file and extract user-id and eye color information and save the output to file
    # NOTE: At this time, this generated file will need to be manually cleaned up and binned into
    # "Blue", "Brown" and "Green" eyed people
    print "Reading user and phenotype information..."
    # # gene_extractor.read_phenotype(DATA_PATH, pheno_input_file_name, pheno_output_file_name)

    # 04. Read the files in users folder and pull out only those that are in the final_eyecolor_phenotypes.csv
    # Read this final users list once again, this is after manual binning of the phenotype
    all_user_pheno_df = pd.read_csv(os.path.join(DATA_PATH, "final_eyecolor_phenotypes.csv"))
    # Add this step to match the user name files in 23&me and Ancestory
    all_user_pheno_df['user_id'] = 'user' + all_user_pheno_df['user_id'].astype(str) + ('_')

    user_file_lst_brown = list()
    user_file_lst_blue_green = list()

    for user in user_file_lst:
        # Extract just the user id from the file name
        short_name = user.split('_')[0]+ '_'
        # Compare it with original dataset
        if all_user_pheno_df['user_id'].str.contains(short_name).any():
            # Get the phenotype information to append to the user name in final.csv
            pheno = (all_user_pheno_df.loc[all_user_pheno_df['user_id'] == short_name].values)[0][1]
            if pheno == "Brown":
                user_file_lst_brown.append(user)
            elif pheno in ("Blue", "Green"):
                user_file_lst_blue_green.append(user)

    lst_colors = ['Brown', 'Blue_Green']

    # 05. Compare if the user id exists in the directory, then read and process it
    for eye_color in lst_colors:
        suffix = "_" + eye_color
        # Process the users with eye_color Brown
        if eye_color == "Brown":
            print "In Brown phenotype loop "
            len_users = genotype_phenotype_extraction(user_file_lst_brown, "Brown", all_gene_data)
            # Merge all the files into one big data file
            # Introducing this step to take care of memory issues. By this step, a large object is released from
            # memory and re-read
            merge_files(suffix)
            # Process the merged file
            final_path, final_grouped_path = write_files(suffix,len_users)
            file_analysis(final_grouped_path, suffix)
        # Process the users with eye_color Blue or Green
        else:
            print "In Blue-green phenotype loop "
            len_users = genotype_phenotype_extraction(user_file_lst_brown, "Blue_Green", all_gene_data)
            # Merge all the files into one big data file
            # Introducing this step to take care of memory issues. By this step, a large object is released from
            # memory and re-read
            merge_files(suffix)
            # Process the merged file
            final_path, final_grouped_path = write_files(suffix,len_users)
            file_analysis(final_grouped_path, suffix)

if __name__  == "__main__":
    main()