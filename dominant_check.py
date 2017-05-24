import pandas as pd
from termcolor import colored

#data_bg<-pd.read_csv('output_final/final_file_Blue_Green.csv')

#data_br<-pd.read_csv('output_final/final_file_Brown.csv')

def dominant_check():
    '''
    Compare between groups which genes are dominant in one and recessive in another
    :return: List of selected dominant genes
    '''
    print "Reading Blue Green Grouped File..."
    data_gbg = pd.read_csv('Data/output_final/final_grouped_Blue_Green.csv')

    print "Reading Brown Grouped File..."
    data_gbr = pd.read_csv('Data/output_final/final_grouped_Brown.csv')

    blue_genes=data_gbg.Gene_info

    brown_genes=data_gbr.Gene_info

    common_genes=set(blue_genes).intersection(brown_genes)

    print "The number of genes in common: ",len(common_genes)
    selected_genes=[]

    for gene in common_genes:
        print "Gene: ",gene

        bg_pct=data_gbg.loc[data_gbg['Gene_info'] == gene, 'percent_users_Blue_Green'].iloc[0]

        #print "Gene's percentage in BG: ",bg_pct

        br_pct=data_gbr.loc[data_gbr['Gene_info'] == gene, 'percent_users_Brown'].iloc[0]

        #print "Gene's percentage in Brown: ", br_pct

        if ((bg_pct>=75 and br_pct<=25) or(bg_pct<=25 and br_pct>=75)):
            print colored(color='red',text="We have a hit!")
            print colored(color='red',text=gene)
            selected_genes+=[gene]

    print "Total Number of selected Genes: ",len(selected_genes)

    print selected_genes

    return selected_genes

if __name__=="__main_":
    dominant_check()

