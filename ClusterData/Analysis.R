
# Converting GWAS --> BP --------------------------------------------------


# gwas = read.csv("mygwas",header = T)
# gwas = gwas %>% filter(p < 0.05) %>% arrange(p)
# gwas_2000 = gwas[1:2000,]
# colnames(gwas_2000) = c("SNP","A1","A2","FREQ","B","SE","P")
# library(biomaRt)
# snpmart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
# listAttributes(snpmart)
# listFilters(snpmart)
# test <- getBM(attributes = c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = list(gwas_2000$SNP), mart = snpmart,useCache = FALSE)
# colnames(test) = c("SNP","CHR","BP")
# gwas_2000 = merge(gwas_2000,test,by = "SNP")
# write.csv(gwas_2000,"Top_2000_IGAP_SNPs.csv",quote = F,row.names = F)
gwas_2000 = read.csv("Top_2000_IGAP_SNPs.csv",header = T)


# BiocManager::install("QUBIC")
library("QUBIC")
results = read.csv("Results_Actual.csv",header = T)
row.names(results) = results$X
results$X = NULL
# BiCluster - QUBIC -------------------------------------------------------
results_mt = as.matrix(results)
RSID_Names = rownames(results_mt)
CHMATIN_Names = colnames(results_mt)
#results_mt = qudiscretize(results_mt)
results_mt[1:5,1:5]
res <- biclust::biclust(results_mt, method = BCQUD(),verbose = TRUE)
summary(res)

quheatmap(results_mt, res, number = c(1,5), showlabel = TRUE)

library(RColorBrewer)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(results_mt, res, number = c(1,2), showlabel = TRUE, col = paleta)

hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
                                 "#E0F3F8", "#91BFDB", "#4575B4")))(100)
par(mar = c(4, 5, 3, 5) + 0.1)
quheatmap(results_mt, res,number = c(1,2), col = hmcols, showlabel = TRUE)

# Checking_Average_PValue -------------------------------------------------
library(QUBIC)
library(biclust)
check_average_sig = function(res_obj,gwas_obj,cluster_number,label_rsid){
  ids_cols = biclusternumber(res_obj,cluster_number)[[1]]$Cols
  ids_rows = biclusternumber(res_obj,cluster_number)[[1]]$Rows
  gwas_res = gwas_obj[which(gwas_obj$SNP %in% label_rsid[ids_rows]),]
  avg_p = mean(gwas_res$P)
  avg_se = mean(gwas_res$SE)
  return(c(avg_p,avg_se))
}

get_snp_list = function(res_obj,labels_crm,cluster_number,label_rsid){
  ids_cols = biclusternumber(res_obj,cluster_number)[[1]]$Cols
  ids_rows = biclusternumber(res_obj,cluster_number)[[1]]$Rows
  snps = label_rsid[ids_rows]
  chrms = labels_crm[ids_cols]
  return(snps)
}
get_crms_list = function(res_obj,labels_crm,cluster_number,label_rsid){
  ids_cols = biclusternumber(res_obj,cluster_number)[[1]]$Cols
  ids_rows = biclusternumber(res_obj,cluster_number)[[1]]$Rows
  snps = label_rsid[ids_rows]
  chrms = labels_crm[ids_cols]
  return(chrms)
}
p_vals = c()
se_vals = c()
for (i in 1:res@Number){
  print(i)
  temp = check_average_sig(res,gwas_2000,i,RSID_Names)
  p_vals = c(p_vals,temp[1])
  se_vals = c(se_vals,temp[2])
}

sig_res = data.frame("ClusterNum"=1:res@Number,"Avg_P" = p_vals,"Avg_SE"=se_vals)

snps_1 = get_snp_list(res,CHMATIN_Names,1,RSID_Names)
crms_1 = get_crms_list(res,CHMATIN_Names,1,RSID_Names)
snps_2 = get_snp_list(res,CHMATIN_Names,2,RSID_Names)
crms_2 = get_crms_list(res,CHMATIN_Names,2,RSID_Names)

library(ggplot2)
ggplot(sig_res) +
  geom_bar( aes(x=ClusterNum, y=Avg_P), stat="identity", fill="skyblue", alpha=0.7) 


# XGR Analysis  -----------------------------------------------------------

# BiocManager::install("XGR")
library(XGR)
