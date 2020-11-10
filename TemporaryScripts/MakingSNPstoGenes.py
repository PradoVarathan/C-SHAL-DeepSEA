# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:13:53 2020

@author: ppugale
"""

import pandas as pd
import os

os.chdir("C:/Users/ppugale/OneDrive - Indiana University/Documents/DeepSEA_Lab")
temp = pd.read_csv("Log_Results.csv")
temp['0'].unique()
snps_top_20 = list(temp['Unnamed: 0'][0:20])
from Bio import Entrez, SeqIO
import pandas as pd

snps = pd.read_csv("SNPs_Thres_0.05",sep=" ")
Entrez.email = "pradluzog@gmail.com"

record = Entrez.read(Entrez.elink(dbfrom="snp", 
                                  id=",".join(snps_top_20).replace('rs', ''), 
                                  db="gene"))

genes = []
for gene_id in record[0]['LinkSetDb'][0]['Link']:
    handle = Entrez.esummary(db="gene", id=gene_id['Id'])
    uid_record = Entrez.read(handle)
    handle.close()

    uid_summary = uid_record["DocumentSummarySet"]['DocumentSummary'][0]
    genes.append(uid_summary['Name'])