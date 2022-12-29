# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np

#"C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_feature_count.txt"
def extract_stats(filename):
    f = pd.read_csv(filename, sep ='\t')
    
    print(f[["# Feature" ,"Class", "Placements"]])

extract_stats("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_faecalis/GCF_004135645.1_ASM413564v1_feature_count.txt")
"C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_genomic.fna"
def length_and_gc(filename):
    f = open(filename ,"r")
    
    line = f.readline()
    total_length = 0
    total_GC = 0
    while line :
        if line.startswith(">"):
            line = f.readline()
        total_GC += line.count('G')
        total_GC += line.count('C')
        total_length += len(line)
        line = f.readline()
        #Eliminate the backspace after each line
        line = line.replace("\n","")
    print(total_GC)
    print("Total length is: " + str(total_length))
    print("% GC = " + str(total_GC * 100 /total_length))
    f.close()

length_and_gc("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_faecalis/GCF_004135645.1_ASM413564v1_genomic.fna")
def cog_categories_EGGNOG(filename):
    f = pd.read_excel(filename,header=2,engine='openpyxl')
    return(f.COG_category.value_counts())
    #When we want to do the calculations of the % of COG for a certain function
    #Divide by the total number of proteins and not total number of queries proteins (NOT len(f))

df3 = pd.concat([cog_categories_EGGNOG("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/collinsella_aerofaciens/results_eggnog/MM_e74keo15.emapper.annotations.xlsx")*100/1934, cog_categories_EGGNOG("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/collinsella_provencensis/results_eggnog/MM_ypkmqep6.emapper.annotations.xlsx")*100/1436,
              cog_categories_EGGNOG('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/collinsella_intestinalis/results_eggnog/MM_dmqw6ljh.emapper.annotations.xlsx')*100/1522, cog_categories_EGGNOG('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_anaerobia/results_eggnog/MM_z2_drto0.emapper.annotations.xlsx')*100/1904,
               cog_categories_EGGNOG('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_faecalis/results_eggnog\MM_onwcr_5u.emapper.annotations.xlsx')*100/2211], axis = 1)
print(df3)
df3 = df3.set_axis(['Collinsella Aerofaciens', 'Collinsella Provencensis', 'collinsella Intestinalis', 'Senegalimassilia Anaerobia', 'Senegalimassilia Faecalis'], axis=1, inplace=False)
df3.to_csv("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/resultats_newesT_comparaison_5.tsv", sep='\t')  
def cog_categories(filename):
    f = pd.read_csv(filename, sep = "\t")
    f.Cat = f.Cat.str.replace(' ', '')
    return(f.Cat.value_counts())
df4 = pd.concat([cog_categories_EGGNOG("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/collinsella_aerofaciens/results_eggnog/MM_e74keo15.emapper.annotations.xlsx")*100/1934,cog_categories("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/collinsella_aerofaciens/cog_result_table.tsv")*100/1934], axis=  1)

df4.to_csv("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/resultats_comparaison_eggNOG_VS_COG.tsv", sep='\t') 

#Coding density = Size of all the CDS / Size of the genome


def CDS_coding_protein_stats(filename):
    f = pd.read_csv(filename, sep ='\t')
    f = f[f["# feature"] == "CDS"]
    f = f[f["class"] == "with_protein"]
    print("Max size of CDS: "+str(max(f["feature_interval_length"])))
    print("AVG size of CDS: "+str(f["feature_interval_length"].mean()))
    print("Total size of CDS: "+str(sum(f["feature_interval_length"])))
    print("Number of CDS: " +str(len(f["feature_interval_length"])))

    return(f)
    
# t = CDS_coding_protein_stats("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_feature_table.tsv")
# CDS_coding_protein_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/proche1_GCF_900078545.1/GCF_900078545.1_PRJEB13207_feature_table.txt')
# CDS_coding_protein_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/proche2_GCF_902501455.1/GCF_902501455.1_DSM13632_feature_table.txt')
#CDS_coding_protein_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_faecalis/GCF_004135645.1_ASM413564v1_feature_table.txt')
#CDS_coding_protein_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_anaerobia/GCF_000236865.1_ASM23686v1_feature_table.txt')

def gene_stats(filename):
    f = pd.read_csv(filename, sep ='\t')
    f = f[f["# feature"] == "gene"]
    print("Max size of genes: "+str(max(f["feature_interval_length"])))
    print("AVG size of genes: "+str(f["feature_interval_length"].mean()))
    print("Total size of genes: "+str(sum(f["feature_interval_length"])))
    print("Number of genes: " +str(len(f["feature_interval_length"])))
    return(f)
    
# t = gene_stats("C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/asleGCA_003856815.1/GCF_003856815.1_ASM385681v1_feature_table.tsv")
# gene_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/proche1_GCF_900078545.1/GCF_900078545.1_PRJEB13207_feature_table.txt')
# gene_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/proche2_GCF_902501455.1/GCF_902501455.1_DSM13632_feature_table.txt')
#gene_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_faecalis/GCF_004135645.1_ASM413564v1_feature_table.txt')
#gene_stats('C:/Users/User/Desktop/COURS M2 DLAD/S1/Geco/tp1/Senegalimassilia_anaerobia/GCF_000236865.1_ASM23686v1_feature_table.txt')