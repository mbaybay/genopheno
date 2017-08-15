__author__ = 'Shraddha Lanka, Patricia Francis-Lyon'

# Requires folder Data/output_intermediate containing:
#  final_snp_br.csv, final_snp_bg.csv
# And folder Data/output_temp containing: bg_snp_pct.csv, br_snp_pct.csv, brown_recessive_snps.csv
#  brown_dominant_snps.csv, blue_recessive_snps.csv, blue_dominant_snps.csv 
# And folder Data/output_final containing final_file_Brown.csv,
#  final_file_Blue_Green.csv

import pandas as pd
import re
import csv
from sklearn.preprocessing import Imputer
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
import sklearn.metrics as skm
from sklearn import tree
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import GridSearchCV
#from sklearn.linear_model import LogisticRegression


from operator import itemgetter

from subprocess import check_call
import pydotplus
pydotplus.find_graphviz()
#import graphviz


# To include that SNP in the set of SNPs exhibiting between group differences
# Minimum difference in percentage points of mutated alleles between different phenotypes for a SNP in order to include a SNP in set of candidates
PERCENT_PT_DIFF_MIN = 5  #10
# Minimum percent difference in percentage of mutated alleles between different phenotypes for a SNP in order to include a SNP in set of candidates
PERCENT_DIFF_MIN = 20 #  15%, 20%, 25%, 30%, 35%, 10% (this is percent difference, not percentage point difference)


def classify_metric(test_y, test_y_pred):
    """
    A function to print the metrics for classification
    -------
    Params
    test_y (Numpy array): The test data provided
    test_y_pred (Numpy array): The test predicted by the model

    returns
    ----------
    A numpy array containing metrics
    """
    confusion_matrix = skm.confusion_matrix(test_y, test_y_pred)
    sensitivity = (float(confusion_matrix[1][1]) / (float(confusion_matrix[1][1]) + float(confusion_matrix[1][0])))
    specificity = float(confusion_matrix[0, 0]) / float(confusion_matrix[0, 0] + confusion_matrix[0, 1])
    accuracy = skm.accuracy_score(test_y, test_y_pred)

    print "Sensitivity: ", np.round(sensitivity, 3)
    print "Specificity: ", np.round(specificity, 3)
    print "Accuracy: ", np.round(accuracy, 3)

    return sensitivity, specificity, accuracy

def data_processing():
    '''
    Use files in output_final folder to use only columns that are needed
    :return: 
    '''
    print "Reading Blue Green ..."
    data_bg = pd.read_csv('Data/output_final/final_file_Blue_Green.csv')

    colnames = data_bg.columns.values
    drop_cols = [s for s in colnames[2:(len(colnames) - 1)] if "gene" in s]

    data_bg.drop(labels=drop_cols, axis=1, inplace=True)
    data_bg.drop(labels=['total_user_per_snp'], axis=1, inplace=True)#in new version of step 1, column doesn't currently exist so future comment this out

    colnames = data_bg.columns.values
    user_id = []
    for val in colnames:
        if 'user' in val:
            val = re.sub("[^A-Z\d]", "", re.search("^[^_]*", val).group(0).upper())
            val = val.split('_')
            val = int(re.search(r'\d+', val[0]).group())
        # print val
        user_id += [val]

    data_bg.columns = user_id
    print data_bg.head()
    data_bg.to_csv('Data/output_intermediate/final_snp_bg.csv',index=False)

    # Brown
    print " Reading Brown..."
    data_br = pd.read_csv('Data/output_final/final_file_Brown.csv')

    colnames = data_br.columns.values
    drop_cols = [s for s in colnames[2:(len(colnames) - 1)] if "gene" in s]

    data_br.drop(labels=drop_cols, axis=1, inplace=True)
    data_br.drop(labels=['total_user_per_snp'], axis=1, inplace=True)#in new version of step 1, column doesn't currently exist so future comment this out

    colnames = data_br.columns.values
    user_id = []
    for val in colnames:
        if 'user' in val:
            val = re.sub("[^A-Z\d]", "", re.search("^[^_]*", val).group(0).upper())
            val = val.split('_')
            val = int(re.search(r'\d+', val[0]).group())
        # print val
        user_id += [val]

    data_br.columns = user_id
    print data_br.head()
    data_br.to_csv('Data/output_intermediate/final_snp_br.csv', index=False)


def pct_computation():
    '''
    :return: 
    '''
    data_bg = pd.read_csv('Data/output_intermediate/final_snp_bg.csv')
    print data_bg.head()

    bg_rows, bg_cols = data_bg.shape

    num_fm = (data_bg.apply(lambda row: sum(row[0:bg_rows - 1] == 2), axis=1))
    num_nm = (data_bg.apply(lambda row: sum(row[0:bg_rows - 1] == 0), axis=1))
    num_pm = (data_bg.apply(lambda row: sum(row[0:bg_rows - 1] == 1), axis=1))

    num_users = num_fm + num_nm + num_pm
    pct_fm = (num_fm / num_users) * 100

    data_bg['pct_fm'] = pct_fm
    del [pct_fm, num_fm]

    pct_nm = (num_nm / num_users) * 100
    data_bg['pct_nm'] = pct_nm
    del [pct_nm, num_nm]

    pct_pm = (num_pm / num_users) * 100
    data_bg['pct_pm'] = pct_pm
    del [pct_pm, num_pm]

    num_3 = (data_bg.apply(lambda row: sum(row[0:bg_rows - 1] == 3), axis=1))
    data_bg['num_3'] = num_3
    del num_3

    num_nan = (data_bg.apply(lambda row: row[0:bg_rows - 1].isnull().sum(), axis=1))
    data_bg['num_nan'] = num_nan
    del num_nan

    print data_bg.head()
    data_bg.to_csv('Data/output_temp/bg_snp_pct.csv', index=False)

    # Brown
    data_br = pd.read_csv('Data/output_intermediate/final_snp_br.csv')
    print data_br.head()
    br_rows, br_cols = data_br.shape

    num_fm = (data_br.apply(lambda row: sum(row[0:br_rows - 1] == 2), axis=1))
    num_nm = (data_br.apply(lambda row: sum(row[0:br_rows - 1] == 0), axis=1))
    num_pm = (data_br.apply(lambda row: sum(row[0:br_rows - 1] == 1), axis=1))

    num_users = num_fm + num_nm +num_pm
    pct_fm = (num_fm / num_users) * 100

    data_br['pct_fm'] = pct_fm
    del [pct_fm,num_fm]

    pct_nm = (num_nm / num_users) * 100
    data_br['pct_nm'] = pct_nm
    del [pct_nm,num_nm]

    pct_pm = (num_pm / num_users) * 100
    data_br['pct_pm'] = pct_pm
    del [pct_pm,num_pm]

    num_3 = (data_br.apply(lambda row: sum(row[0:br_rows - 1] == 3), axis=1))
    data_br['num_3'] = num_3
    del num_3

    num_nan = (data_br.apply(lambda row: row[0:br_rows - 1].isnull().sum(), axis=1))
    data_br['num_nan'] = num_nan
    del num_nan

    print data_br.head()
    data_br.to_csv('Data/output_temp/br_snp_pct.csv', index=False)


def recessive_check(upper, lower):
    data_bg = pd.read_csv('Data/output_temp/bg_snp_pct.csv')
    data_bg.sort_values(by=['pct_fm'], ascending=False, inplace=True)
    bg_fm_gtupper = data_bg.loc[data_bg['pct_fm'] >= upper]

    if len(bg_fm_gtupper) < int(0.01 * len(data_bg)):
        n_rows = int(len(data_bg) * 0.01)
        print n_rows
        bg_fm_gtupper = data_bg.head(n=n_rows)

    bg_fm_ltlower = data_bg.loc[data_bg['pct_fm'] <= lower]

    if len(bg_fm_ltlower) < int(0.01 * len(data_bg)):
        n_rows = int(len(data_bg) * 0.01)
        print n_rows
        bg_fm_ltlower = data_bg.tail(n=n_rows)
    print bg_fm_ltlower

    # Brown
    data_br = pd.read_csv('Data/output_temp/br_snp_pct.csv')
    data_br.sort_values(by=['pct_fm'], ascending=False, inplace=True)
    br_fm_gtupper = data_br.loc[data_br['pct_fm'] >= upper]

    if len(br_fm_gtupper) < int(0.01 * len(data_br)):
        n_rows = int(len(data_br) * 0.01)
        print n_rows
        br_fm_gtupper = data_br.head(n=n_rows)

    br_fm_ltlower = data_br.loc[data_br['pct_fm'] <= lower]
    if len(br_fm_ltlower) < int(0.01 * len(data_br)):
        n_rows = int(len(data_br) * 0.01)
        print n_rows
        br_fm_ltlower = data_bg.tail(n=n_rows)

    # ex: blue_recessive as those Rsids with  ((bg_fm_pct >= 75 and br_fm_pct <= 25))
    blue_recessive = set(bg_fm_gtupper['Rsid']).intersection(br_fm_ltlower['Rsid'])
    brown_recessive = set(br_fm_gtupper['Rsid']).intersection(bg_fm_ltlower['Rsid'])

    csvfile = 'Data/output_temp/blue_recessive_snps.csv'

    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in blue_recessive:
            writer.writerow([val])

    csvfile = 'Data/output_temp/brown_recessive_snps.csv'

    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in brown_recessive:
            writer.writerow([val])


def dominant_check(upper, lower):
    data_bg = pd.read_csv('Data/output_temp/bg_snp_pct.csv')

    print "Dominant Blue Traits..."
    data_bg.sort_values(by=['pct_nm'], ascending=False, inplace=True)

    bg_nm_gtupper = data_bg.loc[data_bg['pct_nm'] >= upper]

    print "Number of SNPs having No mutations > %d: %d" % (upper, len(bg_nm_gtupper))

    if len(bg_nm_gtupper) < int(0.01 * len(data_bg)):
        n_rows = int(len(data_bg) * 0.01)
        print n_rows
        bg_nm_gtupper = data_bg.head(n=n_rows)

    bg_nm_ltlower = data_bg.loc[data_bg['pct_nm'] <= lower]

    print "Number of SNPs having No mutations < %d: %d" % (lower, len(bg_nm_ltlower))

    if len(bg_nm_ltlower) < int(0.01 * len(data_bg)):
        n_rows = int(len(data_bg) * 0.01)
        print n_rows
        bg_nm_ltlower = data_bg.tail(n=n_rows)

    # Brown

    print "Dominant Brown Traits..."
    data_br = pd.read_csv('Data/output_temp/br_snp_pct.csv')

    data_br.sort_values(by=['pct_nm'], ascending=False, inplace=True)

    br_nm_gtupper = data_br.loc[data_br['pct_nm'] >= upper]

    print "Number of SNPs having No mutations > %d: %d" % (upper, len(br_nm_gtupper))

    if len(br_nm_gtupper) < int(0.01 * len(data_br)):
        n_rows = int(len(data_br) * 0.01)
        print n_rows
        br_nm_gtupper = data_br.head(n=n_rows)

    br_nm_ltlower = data_br.loc[data_br['pct_nm'] <= lower]

    print "Number of SNPs having No mutations < %d: %d" % (lower, len(br_nm_ltlower))

    if len(br_nm_ltlower) < int(0.01 * len(data_br)):
        n_rows = int(len(data_br) * 0.01)
        print n_rows
        br_nm_ltlower = data_br.tail(n=n_rows)

    # ex:  brown_dominant as those Rsids with (bg_nm_pct <= 25 and br_nm_pct >= 75)
    brown_dominant = set(bg_nm_ltlower['Rsid']).intersection(br_nm_gtupper['Rsid'])
    blue_dominant = set(br_nm_ltlower['Rsid']).intersection(bg_nm_gtupper['Rsid'])

    csvfile = 'Data/output_temp/blue_dominant_snps.csv'

    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in blue_dominant:
            writer.writerow([val])

    csvfile = 'Data/output_temp/brown_dominant_snps.csv'

    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in brown_dominant:
            writer.writerow([val])


def combine_gene_SNP(s):   # helper function for tree label
    """return part of Gene_info after colon """
    return s.rsplit(':', 1)[1] #  .rsplit() searches for the splitting string from the end of input string, and the second argument limits how many times it'll split to just once


def all_snps_on_gene():
    blue_dominant = pd.read_csv(filepath_or_buffer='Data/output_temp/blue_dominant_snps.csv', header=None,
                                names=['Rsid'])
    blue_recessive = pd.read_csv(filepath_or_buffer='Data/output_temp/blue_recessive_snps.csv', header=None,
                                 names=['Rsid'])
    brown_dominant = pd.read_csv(filepath_or_buffer='Data/output_temp/brown_dominant_snps.csv', header=None,
                                 names=['Rsid'])
    brown_recessive = pd.read_csv(filepath_or_buffer='Data/output_temp/brown_recessive_snps.csv',
                                  header=None,
                                  names=['Rsid'])
    combined = brown_dominant['Rsid'].tolist() + blue_dominant['Rsid'].tolist() + blue_recessive['Rsid'].tolist() + \
               brown_recessive['Rsid'].tolist()
    print "Reading final file Blue Green..."
    data_bg = pd.read_csv('Data/output_intermediate/final_snp_bg.csv')

    #### I don't think we use this
    pheno = pd.read_csv('Data/final_eyecolor_phenotypes.csv')

    blue_genes = data_bg.loc[data_bg['Rsid'].isin(combined),'Gene_info']
    blue_df= data_bg.loc[data_bg['Gene_info'].isin(blue_genes)]
    blue_df.drop_duplicates(subset='Rsid', inplace=True)  

    # LATER: SHOULD DO THIS ONCE ONLY, IN STEP 1, NOT IN STEP 2 ML, which is separately for blue and brown
    # rename SNP so as to carry gene name: we want part of Gene_info after colon concat with Rsid
    gene_SNP_blue =  blue_df['Gene_info'].map(combine_gene_SNP) + blue_df['Rsid']
    print gene_SNP_blue
    blue_df['Rsid'] = gene_SNP_blue

    blue_df.drop(labels='Gene_info', axis=1, inplace=True)

    blue_df_t = blue_df.transpose()
    blue_df_t = blue_df_t.rename(columns=blue_df_t.iloc[0])
    blue_df_t = blue_df_t[1:]
    blue_df_t['user_id'] = blue_df_t.index

    blue_df_t['user_id'] = blue_df_t['user_id'].apply(pd.to_numeric)

    print "Reading final file Brown..."
    data_br = pd.read_csv('Data/output_intermediate/final_snp_br.csv')

    brown_genes = data_br.loc[data_br['Rsid'].isin(combined), 'Gene_info']

    brown_df = data_br.loc[data_br['Gene_info'].isin(brown_genes)]
    brown_df.drop_duplicates(subset='Rsid', inplace=True)  

    # LATER: SHOULD DO THIS ONCE ONLY, IN STEP 1, NOT IN STEP 2 ML, which is separately for blue and brown
    # rename SNP so as to carry gene name: we want part of Gene_info after colon concat with Rsid
    gene_SNP_brown=  brown_df['Gene_info'].map(combine_gene_SNP) + brown_df['Rsid']
    print gene_SNP_brown
    brown_df['Rsid'] = gene_SNP_brown

    brown_df.drop(labels='Gene_info', axis=1, inplace=True)

    brown_df_t = brown_df.transpose()
    brown_df_t = brown_df_t.rename(columns=brown_df_t.iloc[0])
    brown_df_t = brown_df_t[1:]
    brown_df_t['user_id'] = brown_df_t.index

    brown_df_t['user_id'] = brown_df_t['user_id'].apply(pd.to_numeric)
    blue_df_t['phenotype'] = 'Blue'
    brown_df_t['phenotype'] = 'Brown'
    final_dataset=pd.concat([blue_df_t,brown_df_t])

    #
    final_dataset.to_csv("Data/output_intermediate/final_dataset.csv", index=False)
    return "Data/output_intermediate/final_dataset.csv"
    #



def between_group_differences(pct_pt_diff_min, pct_diff_min):
    '''
    Function to return those SNPs that satisfy a criterion to check for differences between blue and brown SNPs
    :param pct_pt_diff_min: Minimum difference in percentage points of mutated alleles between different phenotypes in order to include a SNP in set of candidates
    :param pct_diff_min: # Minimum percent difference in percentage of mutated alleles between different phenotypes for a SNP in order to include a SNP in set of candidates

    :return: 
    '''
    print "Reading BlueGreen File..."
    data_bg = pd.read_csv('Data/output_temp/bg_snp_pct.csv')
    print "Reading Brown File..."
    data_br = pd.read_csv('Data/output_temp/br_snp_pct.csv')

    print "Common SNPS in blue and brown sets"
    unique_snps_blue=set(data_bg['Rsid'])
    unique_snps_brown = set(data_br['Rsid'])
    common_unique_snps = list(set(unique_snps_blue).intersection(unique_snps_brown))
    data_bg_subset= data_bg.loc[data_bg['Rsid'].isin(common_unique_snps),]
    data_br_subset = data_br.loc[data_br['Rsid'].isin(common_unique_snps),]

    selected_snps=[]
    for snp in common_unique_snps:
        bg_pct_nm= np.round(list(data_bg_subset.loc[data_bg_subset['Rsid']==snp,'pct_nm'])[0],3)
        bg_pct_pm = np.round(list(data_bg_subset.loc[data_bg_subset['Rsid']==snp,'pct_pm'])[0],3)
        bg_pct_fm = np.round(list(data_bg_subset.loc[data_bg_subset['Rsid']==snp,'pct_fm'])[0],3)
        br_pct_nm = np.round(list(data_br_subset.loc[data_br_subset['Rsid'] == snp, 'pct_nm'])[0], 3)
        br_pct_pm = np.round(list(data_br_subset.loc[data_br_subset['Rsid'] == snp, 'pct_pm'])[0], 3)
        br_pct_fm = np.round(list(data_br_subset.loc[data_br_subset['Rsid'] == snp, 'pct_fm'])[0], 3)

        min_nm = min(bg_pct_nm,br_pct_nm)
        min_pm = min(bg_pct_pm,br_pct_pm)
        min_fm = min(bg_pct_fm,br_pct_fm)
        diff_pct_nm = abs(bg_pct_nm-br_pct_nm)
        diff_pct_pm = abs(bg_pct_pm-br_pct_pm)
        diff_pct_fm = abs(bg_pct_fm-br_pct_fm) 

        if diff_pct_nm >=pct_pt_diff_min: # min difference in percentage points
            if min_nm !=0: #employ more percentage criteria
                if 100* diff_pct_nm/float(min_nm)>= pct_diff_min:
                    #print "diff_pct_nm as percent of smaller stat for this SNP: ", 100* diff_pct_nm/float(min_nm)
                    selected_snps+=[snp]
            else: # if the smaller is zero mutated alleles one criteria only, no division by zero
                #print "difference of percentage points ", diff_pct_nm
                selected_snps+=[snp]
        elif diff_pct_pm >=pct_pt_diff_min: 
            if min_pm !=0:
                if 100* diff_pct_pm/float(min_pm)>= pct_diff_min:
                    #print "diff_pct_pm as percent of smaller stat for this SNP: ", 100* diff_pct_pm/float(min_pm)
                    selected_snps+=[snp]
            else: 
                #print "difference of percentage points ", diff_pct_nm
                selected_snps+=[snp]
        elif diff_pct_fm >=pct_pt_diff_min:
            if min_fm !=0:
                if 100* diff_pct_fm/float(min_fm)>= pct_diff_min:
                    #print "diff_pct_fm as percent of smaller stat for this SNP: ", 100* diff_pct_fm/float(min_fm)
                    selected_snps+=[snp]
            else: 
                #print "difference of percentage points ", diff_pct_fm
                selected_snps+=[snp]

    print "selected_snps ", "%%%%%%%%%%%%%%%%%%%"
    print selected_snps

    print "Reading final file Blue Green..."
    data_bg = pd.read_csv('Data/output_intermediate/final_snp_bg.csv')

    snp_data_blue = data_bg.loc[data_bg['Rsid'].isin(selected_snps),]

    gene_SNP_blue = snp_data_blue['Gene_info'].map(combine_gene_SNP) + snp_data_blue['Rsid']

    snp_data_blue['Rsid'] = gene_SNP_blue

    snp_data_blue.drop(labels='Gene_info', axis=1, inplace=True)
    snp_data_blue.drop_duplicates(subset='Rsid',inplace=True)

    blue_df_t = snp_data_blue.transpose()
    blue_df_t = blue_df_t.rename(columns=blue_df_t.iloc[0])
    blue_df_t = blue_df_t[1:]
    blue_df_t['user_id'] = blue_df_t.index

    blue_df_t['user_id'] = blue_df_t['user_id'].apply(pd.to_numeric)



    print "Reading final file Brown..."
    data_br = pd.read_csv('Data/output_intermediate/final_snp_br.csv')

    snp_data_brown = data_br.loc[data_br['Rsid'].isin(selected_snps),]

    #snp_data_brown.drop_duplicates(subset='Rsid', inplace=True)   ##############

    gene_SNP_brown = snp_data_brown['Gene_info'].map(combine_gene_SNP) + snp_data_brown['Rsid']
    print gene_SNP_brown
    print "len(gene_SNP_blue)", len(gene_SNP_blue)  
    print "len(gene_SNP_brown)", len(gene_SNP_brown)
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print snp_data_brown['Rsid']
    snp_data_brown['Rsid'] = gene_SNP_brown
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    print snp_data_brown['Rsid']
    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

    snp_data_brown.drop(labels='Gene_info', axis=1, inplace=True)
    snp_data_brown.drop_duplicates(subset='Rsid', inplace=True)

    brown_df_t = snp_data_brown.transpose()
    brown_df_t = brown_df_t.rename(columns=brown_df_t.iloc[0])
    brown_df_t = brown_df_t[1:]
    brown_df_t['user_id'] = brown_df_t.index

    brown_df_t['user_id'] = brown_df_t['user_id'].apply(pd.to_numeric)

    blue_df_t['phenotype'] = 'Blue'

    brown_df_t['phenotype'] = 'Brown'

    print blue_df_t.shape

    print blue_df_t.head()

    print brown_df_t.shape

    print brown_df_t.head()



    final_dataset = pd.concat([blue_df_t, brown_df_t])

    #
    final_dataset.to_csv("Data/output_intermediate/final_dataset_diff.csv", index=False)

    return "Data/output_intermediate/final_dataset_diff.csv"


def ML_forest(file_path):
    '''
    Using Machine Learning on the final dataset mentioned in the file_path
    :return: 
    '''
    final_dataset = pd.read_csv(file_path)
    final_dataset.drop('user_id', axis=1, inplace=True)
    print "Shape before Drop: ",final_dataset.shape

    final_dataset.dropna(axis=1, thresh=int(0.1*len(final_dataset)), inplace=True)

    print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    #X = final_dataset.iloc[:,0:15] ####
    y = final_dataset.loc[:, 'phenotype']

    print y.unique()

    X_col_names = X.columns.values
    print X_col_names ####

    #### ToDo: should calculate the most frequent from training set only, then impute this value to all ####
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X)
    X = imp.transform(X)


    clf = RandomForestClassifier(n_estimators=3000, max_features='auto', n_jobs=-1) #### added more trees
    #### clf = clf.fit(X, y)#### use only the training set, otherwise don't evaluate the model on the test set
        #### and also don't evaluate any test subset of this set on the most important features

    ####X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)#### need to stratify
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1,stratify=y)
    clf = clf.fit(X_train, y_train)

    print(clf.feature_importances_)

    imp_df = pd.DataFrame({"col_name": X_col_names, "importance": np.round(clf.feature_importances_, 5)})
    print imp_df.sort_values(by='importance', ascending=False)
    final_dataset_t = final_dataset.transpose()


    imp_df.to_csv("feature_importances.csv")

    y_pred = clf.predict(X_test)

    classify_metric(y_test, y_pred)

    cnf_matrix = confusion_matrix(y_test,y_pred)
    print cnf_matrix
    ####perhaps sort columns of final_dataset_diff by order of importance and save that version



# Chris Strelioff's code: replace this as it is deprecated
#http://chrisstrelioff.ws/sandbox/2015/06/25/decision_trees_in_python_again_cross_validation.html
def report(grid_scores, n_top=3):
    """Report top n_top parameters settings, default n_top=3.
    Args
    ----
    grid_scores -- output from grid or random search
    n_top -- how many to report, of top models

    Returns
    -------
    top_params -- [dict] top parameter settings found in
                  search
    """
    top_scores = sorted(grid_scores,
                        key=itemgetter(1),
                        reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print(("Mean validation score: "
               "{0:.3f} (std: {1:.3f})").format(
               score.mean_validation_score,
               np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")

    return top_scores[0].parameters


def ML_dtree(file_path):
    '''
    Perform Machine Learning using Decision Trees on the dataset lacated at file_path
    :return: 
    '''
    final_dataset = pd.read_csv(file_path)
    final_dataset.drop('user_id', axis=1, inplace=True)
    #print "Shape before Drop: ",final_dataset.shape
    final_dataset.dropna(axis=1, thresh=int(0.1*len(final_dataset)), inplace=True)
    #print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    y = final_dataset.loc[:, 'phenotype']
    print y.unique()
    X_col_names = X.columns.values
    #print X_col_names

    #### ToDo: should calculate the most frequent from training set only, then impute this value to all ####
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X)
    X = imp.transform(X)

    # cross validation grid search
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1,stratify=y)
    # set of parameters to tune
    print "setting up parameter grid for GridSearchCV"
    param_grid={"criterion": ["gini", "entropy"],
              #If float then min_samples_split is a percentage and ceil(min_samples_split * n_samples)
                # are the minimum number of samples for each split.
              "min_samples_split": [.01, .015, .02, .025],
              "max_depth": [None, 4, 5], # int or None. None allows a full tree
              #If float then min_samples_leaf is a percentage and ceil(min_samples_leaf * n_samples)
                #  are the minimum number of samples for each leaf.
              "min_samples_leaf": [.0025, .005, .01,.015], # Berry and Linoff .0025 to .01
              #If float then max_features is a percentage and int(min_max_features* n_features)
                # features are considered at each split. If 'auto' then max_features = sqrt(n_features)
                # If None then max_features = n_features
              "max_features":[.3,.4,.5,] #A Complete Tutorial on Tree Based Modeling from Scratch (in R & Python)is 30 to 40%
              }

    clf_grid = tree.DecisionTreeClassifier()
    grid_search = GridSearchCV(clf_grid, param_grid=param_grid, cv=5)# ?future: random_state for reproducibility?
       # default refit= True fits the best estimator with the entire dataset: this allows predictions after fit
       # future: predict_proba to calc AUC?
    grid_search.fit(X_train, y_train) # fits all models in the param_grid
    top_params = report(grid_search.grid_scores_, 5) # deprecated, change this ********************
    print top_params

    clf = grid_search.best_estimator_
    clf = clf.fit(X_train, y_train) 
    pred = clf.predict(X_test) # Call predict on the estimator with the best found parameter
    classify_metric(y_test, pred)
    # get training stats to compare for overfitting
    pred = clf.predict(X_train) # Call predict on the estimator with the best found parameter
    classify_metric(y_train, pred)
    print confusion_matrix(y_train, pred)

    #Go to plot to find Rules from training:
    with open("Data/output_temp/dtree.dot", 'w') as f:
        f = tree.export_graphviz(clf, out_file=f,feature_names=X_col_names)

    graph = pydotplus.graphviz.graph_from_dot_file('Data/output_temp/dtree.dot')
    graph.write_png('Data/output_temp/dtree.png')

    y_pred = clf.predict(X_test)
    classify_metric(y_test, y_pred)
    cnf_matrix = confusion_matrix(y_test, y_pred)
    print cnf_matrix



if __name__=='__main__':
    ## Do once only: set up for subsequent runs to select candidate SNPs
    print "preprocessing..."
    #data_processing()
    print "Computing percentages..."
    #pct_computation()

    # select candidate SNPs using within-group commonalities and between-group differences
    print "Selecting candidate SNPs..."
    file_path =  between_group_differences(PERCENT_PT_DIFF_MIN, PERCENT_DIFF_MIN)

    ## Alt select candidate SNPs: use dominant-recessive analysis
    # print "Checking Dominant Traits..."
    # #dominant_check(65,35)
    # print "Checking recessive traits..."
    # #recessive_check(65,35)
    # print "Merging Phenotype..."
    # #file_path=all_snps_on_gene()

    print "Machine Learning..."
    print "PERCENT_PT_DIFF_MIN, PERCENT_DIFF_MIN = ", PERCENT_PT_DIFF_MIN, PERCENT_DIFF_MIN
    ML_forest(file_path)
    ML_dtree(file_path)
