
from Bio import Entrez, SeqIO
import numpy as np
import pandas as pd
import os
import pickle
import torch
import torch.optim as optim
import torch.nn as nn
import torch.utils.data as Data
import torch.nn.functional as F
import time
from sklearn.metrics import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import sys
from pyfiglet import Figlet
f = Figlet(font='slant')
print(f.renderText('C-SHAL'))
print(" --------  CHROMATIN-SEQUENCE HOT ENCODED ANALYSIS PIPLINE WITH DEEPSEA  --------")
import click


@click.command()
@click.option("--ss", prompt="Summary Statistics File",help = "A summary statistics tab seperated file with columns and headers as SNP,CHR,BP,A1,A2. The columns containing Single Nulceotide Polymorphisms, chormosome, base pair locations, primary allele and secondary allele.")
@click.option("--w",prompt = "BP top obtain from SNP position",default = 500, type =int, help="The upstream and downstream width from base pair")
@click.option("--email",prompt="email id",help = "Email for the Entrez ID to obtain sequences")
@click.option("--ak",prompt="API KEY",help="API key")
@click.option("--det",prompt="Detail Level of output (log,all,both)",help="Options for the detail in the output file. log only gives basic log terms;all provides all 919 labels and values; both provides both the files")
def opt_values(ss,w,email,ak,det):
    click.echo("The process will begin now")
    return ss,w,email,ak,det



def one_hot_encode(seq):
    seq = seq.lower()
    mat_list = [np.eye(4)[i] for i in range(4)]
    mapping = dict(zip(['a','c','g','t'],mat_list))
    seq_one_hot =  np.stack([mapping[i] for i in seq]).T   
    return(mapping,seq_one_hot)


# Define the DeepSEA CNN model
class DeepSEA(nn.Module):
    def __init__(self, ):
        super(DeepSEA, self).__init__()
        self.Conv1 = nn.Conv1d(in_channels=4, out_channels=320, kernel_size=8)
        self.Conv2 = nn.Conv1d(in_channels=320, out_channels=480, kernel_size=8)
        self.Conv3 = nn.Conv1d(in_channels=480, out_channels=960, kernel_size=8)
        self.Maxpool = nn.MaxPool1d(kernel_size=4, stride=4)
        self.Drop1 = nn.Dropout(p=0.2)
        self.Drop2 = nn.Dropout(p=0.5)
        self.Linear1 = nn.Linear(53*960, 925)
        self.Linear2 = nn.Linear(925, 919)
    def forward(self, input):
        x = self.Conv1(input)
        x = F.relu(x)
        x = self.Maxpool(x)
        x = self.Drop1(x)
        x = self.Conv2(x)
        x = F.relu(x)
        x = self.Maxpool(x)
        x = self.Drop1(x)
        x = self.Conv3(x)
        x = F.relu(x)
        x = self.Drop2(x)
        x = x.view(-1, 53*960)
        x = self.Linear1(x)
        x = F.relu(x)
        x = self.Linear2(x)
        return x

def log_change(P_ref,P_alt):
  term1 = np.log(P_ref/(1-P_ref))
  term2 = np.log(P_alt/(1-P_alt))
  return abs(term1-term2)

def Run_Deepsea(seq,a1,a2,w):
    seq_ref = seq[0:w] + a1 + seq[(w+1):]
    mapping,seq_one_hot_ref = one_hot_encode(seq_ref)
    x_ref = torch.tensor(seq_one_hot_ref.reshape(1,4,1000),dtype=torch.float).to(device)  #reshaping, converting to tensor and putting on GPU. We nee shape (1,4,1000) since 1st shape represents no of data points
    seq_alt = seq[0:w] + a2 +seq[(w+1):]
    mapping, seq_one_hot_alt = one_hot_encode(seq_alt)
    x_alt =  torch.tensor(seq_one_hot_alt.reshape(1,4,1000),dtype=torch.float).to(device)
    CNN.eval()
    with torch.no_grad(): 
        y_pred_ref =  CNN(x_ref)
        y_pred_alt = CNN(x_alt)
    P_ref =  torch.sigmoid(y_pred_ref.data)   
    P_alt = torch.sigmoid(y_pred_alt.data)
    P_ref = P_ref.detach().cpu().reshape(-1).numpy()  
    P_alt = P_alt.detach().cpu().reshape(-1).numpy()  
    LC = log_change(P_ref,P_alt)
    sorted_LC = np.sort(LC)[::-1]
    return sorted_LC,P_ref,P_alt

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




if __name__ == '__main__':
    ss,w,email,ak,det = opt_values()
    igap_thresholded_snp_list = pd.read_csv(ss,sep=" ")
    sequences = {}
    Entrez.email  = email
    Entrez.api_key = ak
    w = width

    for i in range(0,len(igap_thresholded_snp_list)):
        start_pos = int(igap_thresholded_snp_list["BP"][i]) - 500
        end_pos = int(igap_thresholded_snp_list["BP"][i]) + 499
        chr = str(igap_thresholded_snp_list["CHR"][i])
        rsid = igap_thresholded_snp_list["SNP"][i]
        act_allele = igap_thresholded_snp_list["A1"][i]
        alt_allele = igap_thresholded_snp_list["A2"][i]
        seq = Get_seq(start_pos,end_pos,chr)
        sequences[rsid] = [seq,act_allele,alt_allele]
        


    print("Completed obtaining sequences. Stored in Sequences.pkl")
    output = open('Sequences.pkl','wb')
    pickle.dump(sequences,output)
    output.close()

    CNN = DeepSEA().to(device)
    best_model = torch.load("deepsea_bestmodel.pkl")
    CNN.load_state_dict(best_model)

    cost_function = nn.BCEWithLogitsLoss().to(device)

    seq_list = pickle.load(open("Sequences.pkl","rb"))
    final_results = {}
    detailed_final_results = {}
    for rsid in seq_list.keys():
        print(f"Running ... {rsid}")
        seq = seq_list[rsid][0]
        a1 = seq_list[rsid][1]
        a2 = seq_list[rsid][2]
        if seq != None:
            log_changes,P_ref,P_alt = Run_Deepsea(seq,a1,a2,w)
            detailed_final_results[rsid] = [P_ref,P_alt]
            final_results[rsid] = log_changes

    if det == "both":
        print("Final output files are present in pkl format")
        output = open('Final_Output.pkl','wb')
        pickle.dump(final_results,output)
        output.close()

        output2 = open('Detailed_Final_Output.pkl','wb')
        pickle.dump(detailed_final_results,output2)
        output2.close()
    elif det == "all":
        print("Final output file are present in pkl format")
        output2 = open('Detailed_Final_Output.pkl','wb')
        pickle.dump(detailed_final_results,output2)
        output2.close()
    else:
        print("Final output file are present in pkl format")
        output = open('Final_Output.pkl','wb')
        pickle.dump(final_results,output)
        output.close()
    

