# Simple Way to Read TSV Files in Python using pandas
# importing pandas library
import pandas as pd


def filter_table(filename):
    features_file = pd.read_csv(filename, sep='\t')
    features_file = features_file[features_file["# feature"] == "CDS"]
    return features_file


def extract_information(paf_file, first_feature_table, second_feature_table):
    interest_feature_table = filter_table(first_feature_table)
    tobecompared_feature_table = filter_table(second_feature_table)
    paf_file = pd.read_csv(paf_file, sep='\t')
    paf_file = paf_file.sort_values(
        paf_file.columns[10], ascending=False).head(n=10)
    for index, row in paf_file.iterrows():
        if row.iloc[4] == "+":
            print("Conserved alignment")
        elif row.iloc[4] == "-":
            print("Inversion")
        print(row.iloc[10])
        print(row.iloc[2], row.iloc[3])
        print(interest_feature_table[(interest_feature_table["start"] >= row.iloc[2]) & (
            interest_feature_table["end"] < row.iloc[3])]["name"])
        print(row.iloc[7], row.iloc[8])
        print("These functions")
        print(tobecompared_feature_table[(tobecompared_feature_table["start"] >= row.iloc[7]) & (
            tobecompared_feature_table["end"] < row.iloc[8])]["name"])


extract_information('~/GECO/comparaison_dgenies/interestVSproche1.paf', '~/GECO/proche1_GCF_900199705.1/GCF_900199705.1_PRJEB21894_feature_table.txt',
                    '~/GECO/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_feature_table.tsv')
extract_information('~/GECO/comparaison_dgenies/interestVSproche2.paf', '~/GECO/proche2_GCF_902501455.1/GCF_902501455.1_DSM13632_feature_table.txt',
                    '~/GECO/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_feature_table.tsv')
