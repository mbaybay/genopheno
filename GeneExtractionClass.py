__author__ = 'Lakshmi Arbatti, Shraddha Lanka, Gaurika Tyagi'

import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import re
import numpy as np

class GeneExtraction:
    def batchquery_processing(self, input_path, output_file_name):
        onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]
        snp_info = pd.DataFrame(columns=['Chrom', 'Pos', 'Rsid', 'Ref', 'Alt', 'Gene_info'])

        for files in onlyfiles:
             file_temp = os.path.join(input_path, files)
             if '.gz' in files:
                data = pd.read_csv(file_temp, compression="gzip", skiprows=12, sep='\t',
                                   names=['Chrom', 'Pos', 'Rsid', 'Ref', 'Alt', 'Qual',
                                          'Filter', 'Gene_info'])
             elif '.vcf' in files:
                data = pd.read_csv(file_temp, skiprows=12, sep='\t',
                                   names=['Chrom', 'Pos', 'Rsid', 'Ref', 'Alt', 'Qual',
                                          'Filter', 'Gene_info'])
             else:
                 continue
             data["Gene_info"] = data['Gene_info'].apply(lambda x: re.search('GENEINFO=(.*);d', x).group(1)
                                    if (re.search('GENEINFO=(.*);d', x)) else "")
             data = data.drop('Qual', 1)
             data = data.drop('Filter', 1)
             snp_info = pd.concat([snp_info, data])

        snp_info.to_csv(output_file_name)


    def read_phenotype(self, DATA_PATH, input_file_name, output_file_name):
        """
         This function will read the phenotype raw file and extract all columns which have eye color in them so that all
         these users' genes can be extracted. It will save the raw eye colors and the corresponding user_ids to a file named
         final_eyecolor_phenotypes.csv
         :param filename: string which stores the name of the raw file which contains all phenotypes and user_ids
         :return: None
         """
        phenotypes = pd.read_csv(os.path.join(DATA_PATH, input_file_name), header=0)
        eye_humans1 = phenotypes[["user_id", "Eye.color"]]
        eye_humans1.columns = ["user_id", "eye_color"]
        eye_humans2 = phenotypes[["user_id", "eye.colour"]]
        eye_humans2.columns = ["user_id", "eye_color"]
        all_ids = pd.concat([eye_humans1, eye_humans2], axis=0, join='outer', join_axes=None)
        all_ids = all_ids.dropna().drop_duplicates()
        all_ids.to_csv(os.path.join(DATA_PATH, output_file_name), index=None)

    def build_user_data(self, person_file_path, all_snps):
        # Create an empty dataframe
        data = pd.DataFrame();
        try:
            if "23andme.txt" in person_file_path:
                # print "In 23 & me loop"
                data_person = pd.read_table(person_file_path, header=None, comment="#", low_memory=False,
                                            error_bad_lines=False, delim_whitespace=True, warn_bad_lines=False)
                data_person.columns = ["Rsid", "chromosome", "position", "genotype"]
                data_person = data_person.drop(["chromosome", "position"], axis=1)
            else:
                data_person = pd.read_table(person_file_path, header=0, comment="#", low_memory=False,
                                            error_bad_lines=False, delim_whitespace=True, warn_bad_lines=False)
                data_person.columns = ["Rsid", "chromosome", "position", "allele1", "allele2"]
                data_person["genotype"] = data_person["allele1"] + data_person["allele2"]
                data_person = data_person.loc[data_person["Rsid"].isin(list(all_snps["Rsid"]))]
                data_person = data_person.drop(["chromosome", "position", "allele2", "allele1"], axis=1)
            # else:
            #     return None
            data_person.columns = ["Rsid", 'Genotype']
            data = pd.merge(data_person, all_snps, how="inner", on=["Rsid"], right_index=True)
        except:
            print "Incorrect file format, skipping file...: " , person_file_path
        return data

    def build_gene_data(self,DATA_PATH,gene_input_file,col_split,sep):
        data = pd.read_csv(os.path.join(DATA_PATH, gene_input_file), usecols=[3, 4, 5, 6])
        data = self.tidy_split(data,col_split,sep)
        return data

    def tidy_split(self,data, column, sep='|', keep=False):
        """
        Split the values of a column and expand so the new DataFrame has one split
        value per row. Filters rows where the column is missing.

        Params
        ------
        df : pandas.DataFrame
            dataframe with the column to split and expand
        column : str
            the column to split and expand
        sep : str
            the string used to split the column's values
        keep : bool
            whether to retain the presplit value as it's own row

        Returns
        -------
        pandas.DataFrame
            Returns a dataframe with the same columns as `df`.
        """
        indexes = list()
        new_values = list()
        # df = data.dropna(subset=[column])
        for i, presplit in enumerate(data[column].astype(str)):
            values = presplit.split(sep)
            if keep and len(values) > 1:
                indexes.append(i)
                new_values.append(presplit)
            for value in values:
                indexes.append(i)
                new_values.append(value)
        new_df = data.iloc[indexes, :].copy()
        new_df[column] = new_values
        return new_df

    def find_mutations(self,row):

        changes = 0
        try:
            items = list(row[1])
            ref = row[2]
            alt = list(row[3])
            # If one of the mutated alleles is not in either ref or alternate, ignore this mutation and mark 3
            if items[0] != ref and items [0] not in alt:
                changes = 3
                return changes
            if items[1] != ref and items [1] not in alt:
                changes = 3
                return changes

            if items[0] != ref:
                changes +=1
            if items[1] != ref:
                changes +=1
        except:
            return changes
        return changes


    def allele_transformation(self, user, pheno, data):
        # data [user+'snp'+pheno] = data.apply(self.find_mutations, axis = 1)
        # data[user + 'snp' + pheno] = pd.to_numeric(data.apply(self.find_mutations, axis=1))
        data[user + 'snp' + pheno] = data.apply(self.find_mutations, axis=1)
        data = data.drop(['Genotype' , 'Ref', 'Alt'], axis = 1)
        # data = data.drop(['Genotype'], axis = 1)
        return data