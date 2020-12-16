
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
results = read.csv("~/Documents/GitHub/C-SHAL-DeepSEA/Data/Results_Actual.csv",header = T)
row.names(results) = results$X
results$X = NULL

res_temp = abs(results)
res_temp_logi = res_temp > max(res_temp)*0.35
res_temp = results*res_temp_logi
for (i in 1:nrow(res_temp)){
  if (sum(res_temp[i,]) == 0){
    res_temp[i,] = NA
  }
}
res_temp_r = na.omit(res_temp)
for (i in 1:ncol(res_temp_r)){
  if (sum(res_temp_r[,i]) == 0){
    res_temp_r[,i] = NULL
  }
}
res_temp_r = na.omit(res_temp_r)

library(pheatmap)
dim(res_temp_r)
pheatmap(res_temp_r)
heatmap(as.matrix(res_temp_r))
# BiCluster - QUBIC -------------------------------------------------------
results_mt = as.matrix(results)
RSID_Names = rownames(results_mt)
CHMATIN_Names = colnames(results_mt)
#results_mt = qudiscretize(results_mt)
results_mt[1:5,1:5]
res <- biclust::biclust(results_mt, method = BCQU(),verbose = TRUE)
summary(res)

quheatmap(results_mt, res, number = 1, showlabel = TRUE)

library(RColorBrewer)
paleta <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11)
quheatmap(results_mt, res, number = c(90,84), showlabel = TRUE, col = paleta)

hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
                                 "#E0F3F8", "#91BFDB", "#4575B4")))(100)
par(mar = c(4, 5, 3, 5) + 0.1)
quheatmap(results_mt, res,number = c(3,6), col = hmcols, showlabel = TRUE)

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

get_snp_p_vals = function(snp_list,gwas){
  gwas_res = gwas[which(gwas$SNP %in% snp_list),c('SNP','P')]
  return(gwas_res)
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
library(ggplot2)
ggplot(sig_res) +
  geom_bar( aes(x=ClusterNum, y=Avg_P), stat="identity", fill="skyblue", alpha=0.7) 

snps_3 = get_snp_list(res,CHMATIN_Names,34,RSID_Names)
crms_3 = get_crms_list(res,CHMATIN_Names,34,RSID_Names)
snps_6 = get_snp_list(res,CHMATIN_Names,84,RSID_Names)
crms_6 = get_crms_list(res,CHMATIN_Names,84,RSID_Names)

sig_res = sig_res %>% arrange(Avg_P)
jpeg("All_Clusters.jpeg")
par(mfrow=c(2,4))
for (i in 1:8){
  quheatmap(results_mt, res,number = sig_res$ClusterNum[i], col = hmcols, showlabel = TRUE)
  
}
dev.off()
# XGR Analysis  -----------------------------------------------------------

# BiocManager::install("remotes", dependencies=T)
# BiocManager::install("hfang-bristol/XGR", dependencies=T)
library(XGR)
library(RCircos)

# Cluster3
XGR_3 = xEnricherSNPs(snps_3,include.LD = 'EUR')
genes_3 = xSNP2nGenes(snps_3)
view_3 = xEnrichViewer(XGR_3)
write.csv(view_3,"XGR_Clus34.csv",quote=F,row.names = F)
data = data.frame("SNP"=gwas_2000$SNP[which(gwas_2000$SNP %in% snps_3)],"P" = gwas_2000$P[which(gwas_2000$SNP %in% snps_3)])
ld_scores_3 = xSNPscores(data)
write.csv(ld_scores_3,"LD_Scores_Clus34.csv",quote = F,row.names = F)
social_snps_3 = xSocialiserSNPs(snps_3)
xCircos(social_snps_3)
xSocialiserDAGplot(social_snps_3)
xEnrichChord(XGR_3)
snps_genes_con_ld_3 = xSubneterSNPs(data)
plot(snps_genes_con_ld_3)

# Cluster6
XGR_6 = xEnricherSNPs(snps_6)
genes_6 = xEnrichViewer(snps_6)
view_6 = xEnrichViewer(XGR_6)
ld_scores_6 = xSNPscores(snps_6)
social_snps_6 = xSocialiserSNPs(snps_6)
xCircos(social_snps_6)
xSocialiserDAGplot(social_snps_6)
xEnrichChord(XGR_6)
xSubneterSNPs(snps_6)
clus6_p = get_snp_p_vals(snps_6,gwas_2000)
snps_genes_con_ld_6 = xSubneterSNPs(clus6_p)
plot(snps_genes_con_ld_6)
