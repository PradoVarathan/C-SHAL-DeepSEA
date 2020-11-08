"""
Created on Wed Nov  4 13:30:44 2020

@author: ppugale
"""

from Bio import Entrez, SeqIO
import os
import pandas as pd
import pickle 

igap_thresholded_snp_list = pd.read_csv("SNPs_Thres_0.05",sep=" ")
sequences = {}
Entrez.email  = "*@gmail.com"
Entrez.api_key = "*"

def Get_seq(start_pos,end_pos,chr):
    
    try:
        if int(chr) < 10:
            id_chr = "".join(["NC_00000",chr])
        else:
            id_chr = "".join(["NC_0000",chr])

        handle = Entrez.efetch(db="nucleotide",
                        id = id_chr,
                        rettype = "fasta",
                        strand = 1,
                        seq_start = start_pos,
                        seq_stop  = end_pos)
        record = SeqIO.read(handle,"fasta")
        
        return str(record.seq)
    except:
        print("No proper chromosome found ...")


for i in range(0,len(igap_thresholded_snp_list)):
    print("Running ...")
    start_pos = int(igap_thresholded_snp_list["BP"][i]) - 500
    end_pos = int(igap_thresholded_snp_list["BP"][i]) + 499
    chr = str(igap_thresholded_snp_list["CHR"][i])
    rsid = igap_thresholded_snp_list["SNP"][i]
    act_allele = igap_thresholded_snp_list["A1"][i]
    alt_allele = igap_thresholded_snp_list["A2"][i]
    seq = Get_seq(start_pos,end_pos,chr)
    sequences[rsid] = [seq,act_allele,alt_allele]
    print(i)


   
output = open('Sequences.pkl','wb')
pickle.dump(sequences,output)
output.close()
