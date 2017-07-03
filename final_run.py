import pandas as pd
import re
import csv
from sklearn.preprocessing import Imputer
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.model_selection import train_test_split
import sklearn.metrics as skm
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn import tree

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

    data_bg.drop(labels=['total_user_per_snp'], axis=1, inplace=True)

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

    data_br.drop(labels=['total_user_per_snp'], axis=1, inplace=True)

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

        # if ((bg_fm_pct >= 75 and br_fm_pct <= 25))

    blue_recessive = set(bg_fm_gtupper['Rsid']).intersection(br_fm_ltlower['Rsid'])

    # if(br_fm_pct>=75 and bg_fm_pct<=25)

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

    # if (bg_nm_pct <= 25 and br_nm_pct >= 75)

    brown_dominant = set(bg_nm_ltlower['Rsid']).intersection(br_nm_gtupper['Rsid'])

    # if (br_nm_pct <= 25 and bg_nm_pct >= 75)

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

def data_pheno_merge():
        '''
        Pre-process the data to combine brown and blue-green dataframes, make SNPs columns instead of rows, and add the 
        phenotype column

        :return: 
        '''

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

        pheno = pd.read_csv('Data/final_eyecolor_phenotypes.csv')

        print pheno.head()

        blue_df = data_bg.loc[data_bg['Rsid'].isin(combined)]

        # print blue_df.head()
        # blue_df.to_csv('blue_df.csv')
        blue_df.drop(labels='Gene_info', axis=1, inplace=True)
        blue_df.drop_duplicates(subset='Rsid', inplace=True)

        blue_df_t = blue_df.transpose()

        blue_df_t = blue_df_t.rename(columns=blue_df_t.iloc[0])
        blue_df_t = blue_df_t[1:]
        blue_df_t['user_id'] = blue_df_t.index

        blue_df_t['user_id'] = blue_df_t['user_id'].apply(pd.to_numeric)

        print blue_df_t['user_id'].head()
        # print blue_df_t.iloc[:,0]

        # blue_df_t.to_csv('Data/output_temp/blue_df_t.csv')

        print "Reading final file Brown..."
        data_br = pd.read_csv('Data/output_intermediate/final_snp_br.csv')

        colnames = data_br.columns.values

        brown_df = data_br.loc[data_br['Rsid'].isin(combined)]

        # print brown_df.head()
        # brown_df.to_csv('brown_df.csv')
        brown_df.drop(labels='Gene_info', axis=1, inplace=True)

        brown_df.drop_duplicates(subset='Rsid', inplace=True)

        brown_df_t = brown_df.transpose()
        brown_df_t = brown_df_t.rename(columns=brown_df_t.iloc[0])
        brown_df_t = brown_df_t[1:]

        brown_df_t['user_id'] = brown_df_t.index

        brown_df_t['user_id'] = brown_df_t['user_id'].apply(pd.to_numeric)

        print brown_df_t.head()

        # brown_df_t.to_csv('Data/output_temp/brown_df_t.csv')


        final_dataset_blue = blue_df_t.merge(pheno, how="inner", left_on='user_id', right_on='user_id')

        final_dataset_blue.phenotype[final_dataset_blue['phenotype']=='Green']='Blue'

        print final_dataset_blue.head()

        final_dataset_brown = brown_df_t.merge(pheno, on="user_id")
        #
        print final_dataset_brown.head()
        #
        final_dataset = pd.concat([final_dataset_blue, final_dataset_brown])
        #
        final_dataset.to_csv("Data/output_intermediate/final_dataset.csv", index=False)
        #
        print final_dataset.head()

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

    pheno = pd.read_csv('Data/final_eyecolor_phenotypes.csv')

    #print pheno.head()

    blue_genes = data_bg.loc[data_bg['Rsid'].isin(combined),'Gene_info']
    #print blue_genes

    blue_df= data_bg.loc[data_bg['Gene_info'].isin(blue_genes)]

    #print blue_df.head()
    # blue_df.to_csv('blue_df.csv')
    blue_df.drop(labels='Gene_info', axis=1, inplace=True)
    blue_df.drop_duplicates(subset='Rsid', inplace=True)

    blue_df_t = blue_df.transpose()

    blue_df_t = blue_df_t.rename(columns=blue_df_t.iloc[0])
    blue_df_t = blue_df_t[1:]
    blue_df_t['user_id'] = blue_df_t.index

    blue_df_t['user_id'] = blue_df_t['user_id'].apply(pd.to_numeric)

    #print blue_df_t['user_id'].head()
    # print blue_df_t.iloc[:,0]

    # blue_df_t.to_csv('Data/output_temp/blue_df_t.csv')

    print "Reading final file Brown..."
    data_br = pd.read_csv('Data/output_intermediate/final_snp_br.csv')



    brown_genes = data_br.loc[data_br['Rsid'].isin(combined), 'Gene_info']
    print brown_genes

    brown_df = data_br.loc[data_br['Gene_info'].isin(brown_genes)]

    # print brown_df.head()
    # brown_df.to_csv('brown_df.csv')
    brown_df.drop(labels='Gene_info', axis=1, inplace=True)

    brown_df.drop_duplicates(subset='Rsid', inplace=True)

    brown_df_t = brown_df.transpose()
    brown_df_t = brown_df_t.rename(columns=brown_df_t.iloc[0])
    brown_df_t = brown_df_t[1:]

    brown_df_t['user_id'] = brown_df_t.index

    brown_df_t['user_id'] = brown_df_t['user_id'].apply(pd.to_numeric)

    #print brown_df_t.head()

    # brown_df_t.to_csv('Data/output_temp/brown_df_t.csv')


    final_dataset_blue = blue_df_t.merge(pheno, how="inner", left_on='user_id', right_on='user_id')

    final_dataset_blue.phenotype[final_dataset_blue['phenotype'] == 'Green'] = 'Blue'

    #print final_dataset_blue.head()

    final_dataset_brown = brown_df_t.merge(pheno, on="user_id")
    #
    #print final_dataset_brown.head()
    #
    final_dataset = pd.concat([final_dataset_blue, final_dataset_brown])
    #
    final_dataset.to_csv("Data/output_intermediate/final_dataset.csv", index=False)
    #
    #print final_dataset.head()


def ML_forest():
    '''

    :return: 
    '''
    final_dataset = pd.read_csv('Data/output_intermediate/final_dataset.csv')
    final_dataset.drop('user_id', axis=1, inplace=True)
    print "Shape before Drop: ",final_dataset.shape

    final_dataset.dropna(axis=1, thresh=int(0.1*len(final_dataset)), inplace=True)

    print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    y = final_dataset.loc[:, 'phenotype']

    X_col_names = X.columns.values
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X)

    X = imp.transform(X)
    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X, y)

    print(clf.feature_importances_)

    imp_df = pd.DataFrame({"col_name": X_col_names, "importance": np.round(clf.feature_importances_, 5)})
    print imp_df.sort_values(by='importance', ascending=False)

    imp_df.to_csv("feature_importances.csv")
    imp_df.to_csv("feature_importances_ft.csv")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    classify_metric(y_test, y_pred)

def ML_dtree():
    '''

    :return: 
    '''
    final_dataset = pd.read_csv('Data/output_intermediate/final_dataset.csv')
    final_dataset.drop('user_id', axis=1, inplace=True)
    print "Shape before Drop: ",final_dataset.shape

    final_dataset.dropna(axis=1, thresh=int(0.1*len(final_dataset)), inplace=True)

    print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    y = final_dataset.loc[:, 'phenotype']

    X_col_names = X.columns.values
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X)

    X = imp.transform(X)
    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X, y)

    print(clf.feature_importances_)

    imp_df = pd.DataFrame({"col_name": X_col_names, "importance": np.round(clf.feature_importances_, 5)})
    print imp_df.sort_values(by='importance', ascending=False)

    imp_df.to_csv("feature_importances.csv")
    imp_df.to_csv("feature_importances_ft.csv")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    classify_metric(y_test, y_pred)


def Ml_100():
    final_dataset = pd.read_csv('Data/output_intermediate/final_dataset.csv')
    final_dataset.drop('user_id', axis=1, inplace=True)

    final_dataset.dropna(axis=1, thresh=int(0.1 * len(final_dataset)), inplace=True)

    print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    y = final_dataset.loc[:, 'phenotype']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    X_col_names = X.columns.values
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X_train)

    X_train = imp.transform(X_train)
    clf = ExtraTreesClassifier()
    clf = clf.fit(X_train, y_train)

    print(clf.feature_importances_)

    imp_df = pd.DataFrame({"SNP": X_col_names, "importance": np.round(clf.feature_importances_, 5)})
    print imp_df.sort_values(by='importance', ascending=False)

    model = SelectFromModel(clf, prefit=True)



    X_new = model.transform(X)

    print X_new.shape

    X_test = imp.transform(X_test)

    X_test




    chosen_snps=imp_df['SNP'].head(n=100)



    #clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    #clf = clf.fit(X_new, y)

    X_train, X_test, y_train, y_test = train_test_split(X_new, y, test_size=0.33, random_state=42)

    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    classify_metric(y_test, y_pred)


def fine_tune():
    '''
    
    :return: 
    '''
    importance= pd.read_csv('feature_importances_ft.csv')
    snp = importance.loc[importance['importance']>0,'col_name']

    snp=snp.tolist()

    snp= snp + ['phenotype']

    final= pd.read_csv('Data/output_intermediate/final_dataset.csv')
    final_dataset=final.loc[:,final.columns.isin(snp)]

    rows,cols = final_dataset.shape



#    final_dataset.dropna(axis=1, thresh=0.8, inplace=True)

    print "Shape after drop: ", final_dataset.shape

    X = final_dataset.drop(labels=['phenotype'], axis=1)
    y = final_dataset.loc[:, 'phenotype']

    X_col_names = X.columns.values
    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
    imp = imp.fit(X)

    X = imp.transform(X)
    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X, y)

    print(clf.feature_importances_)

    imp_df = pd.DataFrame({"col_name": X_col_names, "importance": np.round(clf.feature_importances_, 5)})
    #print imp_df.sort_values(by='importance', ascending=False)

    imp_df.to_csv("feature_importances_ft.csv")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    clf = RandomForestClassifier(n_estimators=1000, max_features='auto', n_jobs=-1)
    clf = clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    sensitivity, specificity, accuracy = classify_metric(y_test, y_pred)

    return sensitivity,specificity,accuracy,cols


if __name__=='__main__':
    print "preprocessing..."
    #data_processing()
    print "Computing percentages..."
    #pct_computation()
    print "Checking Dominant Traits..."
    dominant_check(65,35)
    print "Checking recessive traits..."
    recessive_check(65,35)
    print "Merging Phenotype..."
    #data_pheno_merge()
    all_snps_on_gene()
    print "Machine Learning..."
    ML_forest()
    #Ml_100()
    # sensitivity,specificity,accuracy,cols=fine_tune()
    # count = 0
    # while accuracy<=0.78 and count < 100:
    #     count +=1
    #     print "accuracy: , num_features:"
    #     print accuracy, cols - 1
    #     sensitivity, specificity, accuracy, cols = fine_tune()
    #
    # print "final accuracy: %f, final number of features: %d"%(accuracy,cols-1)
