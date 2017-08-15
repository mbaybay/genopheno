__author__ = 'Lakshmi Arbatt, Patricia Francis-Lyon'

# take out the rest of the references to Brown and Blue-green phenotypes (use generic phenotypes)
#later take out the gene level stuff or fix it

import os
import pandas as pd
from GeneExtractionClass import GeneExtraction
import numpy as np

# Percent of missing values (NaN, 3) that, if exceeded for a phenotype, will invalidate a given SNP for that phenotype
#  consider reducing MAX_INVALID_PCT for set sizes below 40
MAX_INVALID_PCT = .9
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
#user_folder = "user_subset"  
user_folder = "user"
user_path = os.path.join(DATA_PATH,user_folder)
output_inter = os.path.join(DATA_PATH,"output_intermediate")
# out_path_final = os.path.join(DATA_PATH,"output_subset")
out_path = os.path.join(DATA_PATH,"output_temp")
# out_path = os.path.join(DATA_PATH,"output_temp_subset")
output_final = os.path.join(DATA_PATH,"output_final")

gene_extractor = GeneExtraction()

# This function reads each snp per row and converts to 0 or 1 based on if there has been a mutation or not
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
    # write = True
    print "len(usr_file_lst) : " , len(usr_file_lst)
    print "users : " , usr_file_lst
    for user in usr_file_lst:
        counter+=1
        print "Counter : ", counter, len(usr_file_lst)
        write = False
        # Extract just the user id from the file name
        short_name = user.split('_y')[0] + '_'
        # Added this conditional to take care of other file types other than 23&me and ancestory.ext
        # if "23andme.txt" in user or "ancestry.txt" in user:
        print "Processing file : ", user
        # If present, read the file and build allele transformation info
        full_file_name = os.path.join(user_path,user)
        data = gene_extractor.build_user_data(full_file_name,all_gene_data)

        # Take care of incorrect file formats
        if not data.empty:   #  merge this user into the data corpus
            # Number of actual users
            len_users +=1
            data = gene_extractor.allele_transformation(short_name, pheno, data)#### data is overwritten with modified data
            # Merge with the earlier dataset to build the data corpus
            merged = pd.merge(merged,data,on=['Rsid','Gene_info'], how='outer' )
            # Get the snp count for each user
            merged[short_name+'gene'+ pheno] = (merged [short_name + 'snp' + pheno].apply(snp_count)).astype(np.uint8)
            merged[short_name + 'snp' + pheno].fillna(255,inplace=True)
            merged[short_name+'snp'+ pheno] = merged[short_name+'snp'+ pheno].astype(np.uint8)
            print "Size of merged : ", merged.size
            print "Files merged : ", len_users
        else:    #### this user's data is empty, so is not merged with dataset nor written to any files  ####
            print  short_name, " data is empty"  ####
    print "Final number of users with valid data for phenotype  : ", phenotype , "  is : ",  len_users
    return len_users, merged



def write_final(suffix,len_users,dataset):
    
    # Get the list of columns
    df_col_lst = list(dataset)
    #### Delete rows with more than MAX_INVALID_PCT invalid or missing  (ie: 3 or 255) SNP data
    jump_list_SNPs = df_col_lst[2::2] #### list of user_SNP values {0,1,2,3, 255=NaN}
    given_set = {3, 255}
    ####print dataset[jump_list_SNPs],  "%%%%%%%%%%%%%%%%%%"
    vec =  dataset[jump_list_SNPs].isin(given_set).sum(1)  # count how many  NaNs, 3's  for each SNP (ie: row)
    #print vec, "%%%%%%%%%%%%%%%%%%"
    dataset = dataset[vec <= MAX_INVALID_PCT * dataset[jump_list_SNPs].shape[1] ]  # keep those rows (SNPs) with enough valid data
    ####print dataset[jump_list_SNPs].shape[1], MAX_INVALID_PCT * dataset[jump_list_SNPs].shape[1], "%%%%%%%%%%%%%%"
    ####
    dataset = dataset.replace(255,'NaN')

    dataset = dataset.set_index(['Rsid'])
    ####print dataset
    
    # Save to csv
    try:
        print "Writing file : ", os.path.join(output_final, "final_file" + suffix + ".csv")
        with open(os.path.join(output_final, "final_file" + suffix + ".csv"),'a') as final_out:
            dataset.to_csv(final_out,header=True)
    except IOError as e:
            print "Error writing final output file ", e


def main():

    # Extract only 23&me and ancestory.com files
    user_lst = os.listdir(user_path)
    user_file_lst = list()
    for file in user_lst:
        if file.endswith("23andme.txt") or file.endswith("ancestry.txt"):
            # print "File name: " , file
            user_file_lst.append(file)
    gene_extractor = GeneExtraction()

    # 01. Use Shraddha's 'batchquery_processing' logic
    print "Processing batch query..."
    gene_extractor.batchquery_processing(os.path.join(DATA_PATH,db_snps_path),main_gene_file)

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

    # Compare with the binned user list and build a list of valid users, by Brown and Blue_Green phenotypes
    for user in user_file_lst:
        # Extract just the user id from the file name
        short_name = user.split('_')[0]+ '_'
        # print "Short Name : ", short_name
        # Compare it with original dataset
        if all_user_pheno_df['user_id'].str.contains(short_name).any():
            # Get the phenotype information to append to the user name in final.csv
            pheno = (all_user_pheno_df.loc[all_user_pheno_df['user_id'] == short_name].values)[0][1]
            # print "Pheno : ", pheno
            if pheno == "Brown":
                user_file_lst_brown.append(user)
            elif pheno in ("Blue", "Green"):
                user_file_lst_blue_green.append(user)

    lst_pheno = ['Brown', 'Blue_Green']

    # 05. Loop to process users of each phenotype in lst_pheno
    for pheno in lst_pheno:
        suffix = "_" + pheno
        # Process the users with phenotype pheno
        print "In " + pheno + " phenotype loop "
        len_users, dataset = genotype_phenotype_extraction(user_file_lst_brown, pheno, all_gene_data)
        write_final(suffix,len_users,dataset)


if __name__  == "__main__":
    main()