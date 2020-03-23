library(ape)
library(phytools)
library(alluvial)
setwd("/Volumes/data/MPV/projects/SCCmec/CCR")
setwd("E:/Github/SCCmec/CCR")
setwd("/Users/thej/Documents/GitHub/SCCmec/CCR")

sp_count_table = read.table("/Volumes/data/DB/refseq/Staphylococcus_species_counts_191128.txt",sep = "\t")
sp_count_table = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/Staphylococcus_species_counts_191128.txt",sep = "\t")
sp_count_vec = sp_count_table$V2
names(sp_count_vec) = sp_count_table$V1

triplet_color_vec = c("#61d2ff","#468bfa","#030bfc","#63ff7d","#23db64","#009133","#ff9999","#f03232","#c90000",
                      "#ffe194","#ffcc47","#e8a800","#f2a8ff","#e44dff","#9e00ba","#f8fc6a","#b8bd00","#786000",
                      "#c4c4c4","#545454","#000000")



#### ccr all ####

tbl = read.table("protein_blast/ccr_all_with_IWG_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl = read.table("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_all_with_IWG_uniq_table.txt?raw=true",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")

tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))

aln = read.dna("protein_blast/ccr_all_uniq_mafft.fasta","fasta")
aln = read.dna("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_all_uniq_mafft.fasta?raw=true","fasta")

dist_mat = read.table("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_muscle_percent_distance.txt?raw=true",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)




# ccrA #

tbl = tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrA_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))


dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrA_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
dist_mat = as.matrix(read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrA_all_uniq_sim.txt",sep = "\t",header=T, row.names=1))

aln_dist = as.dist(dist_mat)


fit = hclust(aln_dist,method = "complete")

v = c()
for (i in 10:40) {
  g = cutree(fit,h=i)
  v = c(v,length(unique(g)))
}
v_A = v


plot(10:40,v)
ccr_groups = cutree(fit,h=22)
rect.hclust(fit,h=22)
tbl$ccrA_group = NA
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrA_group[which(tbl$uniq_fasta_ID==name)] = group
}
tbl$ccr_allotype = paste0('ccrAn',as.vector(tbl$ccrA_group))
tbl$ccr_type = "ccrA"


IWG_ccrA_tbl = tbl[grep('ccrA',tbl$fasta_ID),]
table(as.vector(IWG_ccrA_tbl$fasta_ID),as.vector(IWG_ccrA_tbl$ccrA_group))
sort(unique(ccr_groups))

ccrA_tbl = tbl[which(tbl$ccr_type=="ccrA"),]
#allotype_colors = c(RColorBrewer::brewer.pal(12,"Paired"),"#d9d9d9")
allotype_colors = triplet_color_vec[1:19]
ccr_allotype_vec = unique(ccrA_tbl$ccr_allotype)
ccrA_tbl$allotype_color = NA

for (i in 1:length(allotype_colors)) {
  col = allotype_colors[i]
  ccr_type = ccr_allotype_vec[i]
  idx = which(ccrA_tbl$ccr_allotype==ccr_type)
  names = as.vector(ccrA_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrA_tbl$ccr_allotype[idx] = t
  }
  ccrA_tbl$allotype_color[idx] = col
}

ccrA_uniq_tbl = ccrA_tbl[which(!duplicated(as.vector(ccrA_tbl$uniq_fasta_ID))),]

write.table(ccrA_tbl,"ccrA_allotype_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrA_uniq_tbl,"ccrA_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$allotype_color,' ',ccrA_uniq_tbl$ccr_allotype)
writeLines(to_print,con = "ccrA_allotype_colors_QC_22.txt")

ccrA_dist_phylo = as.phylo(fit)
write.tree(ccrA_dist_phylo,file = "ccrA_dist_tree_QC.nwk")



#### ccrB ####

tbl = tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrB_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))


dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrB_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
dist_mat = as.matrix(read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrB_all_uniq_sim.txt",sep = "\t",header=T, row.names=1))

aln_dist = as.dist(dist_mat)


fit = hclust(aln_dist,method = "complete")

v = c()
for (i in 10:40) {
  g = cutree(fit,h=i)
  v = c(v,length(unique(g)))
}
v_B = v


plot(10:40,v)
ccr_groups = cutree(fit,h=22)
rect.hclust(fit,h=22)
tbl$ccrB_group = NA
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrB_group[which(tbl$uniq_fasta_ID==name)] = group
}
tbl$ccr_allotype = paste0('ccrBn',as.vector(tbl$ccrB_group))
tbl$ccr_type = "ccrB"


IWG_ccrB_tbl = tbl[grep('ccrB',tbl$fasta_ID),]
table(as.vector(IWG_ccrB_tbl$fasta_ID),as.vector(IWG_ccrB_tbl$ccrB_group))
sort(unique(ccr_groups))

ccrB_tbl = tbl[which(tbl$ccr_type=="ccrB"),]
#allotype_colors = c(RColorBrewer::brewer.pal(12,"Paired"),"#d9d9d9")
allotype_colors = triplet_color_vec[1:16]
ccr_allotype_vec = unique(ccrB_tbl$ccr_allotype)
ccrB_tbl$allotype_color = NA

for (i in 1:length(allotype_colors)) {
  col = allotype_colors[i]
  ccr_type = ccr_allotype_vec[i]
  idx = which(ccrB_tbl$ccr_allotype==ccr_type)
  names = as.vector(ccrB_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrB_tbl$ccr_allotype[idx] = t
  }
  ccrB_tbl$allotype_color[idx] = col
}

ccrB_uniq_tbl = ccrB_tbl[which(!duplicated(as.vector(ccrB_tbl$uniq_fasta_ID))),]

write.table(ccrB_tbl,"ccrB_allotype_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrB_uniq_tbl,"ccrB_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrB_uniq_tbl$uniq_fasta_ID,' ',ccrB_uniq_tbl$allotype_color,' ',ccrB_uniq_tbl$ccr_allotype)
writeLines(to_print,con = "ccrB_allotype_colors_QC_22.txt")

ccrB_dist_phylo = as.phylo(fit)
write.tree(ccrB_dist_phylo,file = "ccrB_dist_tree_QC.nwk")




#### ccrC ####

tbl = tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrC_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))


dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrC_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
dist_mat = as.matrix(read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrC_all_uniq_sim.txt",sep = "\t",header=T, row.names=1))

aln_dist = as.dist(dist_mat)


fit = hclust(aln_dist,method = "complete")

v = c()
for (i in 10:40) {
  g = cutree(fit,h=i)
  v = c(v,length(unique(g)))
}
v_C = v


plot(10:40,v)
ccr_groups = cutree(fit,h=22)
rect.hclust(fit,h=22)
tbl$ccrC_group = NA
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrC_group[which(tbl$uniq_fasta_ID==name)] = group
}
tbl$ccr_allotype = paste0('ccrCn',as.vector(tbl$ccrC_group))
tbl$ccr_type = "ccrC"


IWG_ccrC_tbl = tbl[grep('ccrC',tbl$fasta_ID),]
table(as.vector(IWG_ccrC_tbl$fasta_ID),as.vector(IWG_ccrC_tbl$ccrC_group))
sort(unique(ccr_groups))

ccrC_tbl = tbl[which(tbl$ccr_type=="ccrC"),]
allotype_colors = RColorBrewer::brewer.pal(6,"Set1")
#allotype_colors = triplet_color_vec[1:16]
ccr_allotype_vec = unique(ccrC_tbl$ccr_allotype)
ccrC_tbl$allotype_color = NA

for (i in 1:length(allotype_colors)) {
  col = allotype_colors[i]
  ccr_type = ccr_allotype_vec[i]
  idx = which(ccrC_tbl$ccr_allotype==ccr_type)
  names = as.vector(ccrC_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrC_tbl$ccr_allotype[idx] = t
  }
  ccrC_tbl$allotype_color[idx] = col
}

ccrC_uniq_tbl = ccrC_tbl[which(!duplicated(as.vector(ccrC_tbl$uniq_fasta_ID))),]

write.table(ccrC_tbl,"ccrC_allotype_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrC_uniq_tbl,"ccrC_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrC_uniq_tbl$uniq_fasta_ID,' ',ccrC_uniq_tbl$allotype_color,' ',ccrC_uniq_tbl$ccr_allotype)
writeLines(to_print,con = "ccrC_allotype_colors_QC_22.txt")

ccrC_dist_phylo = as.phylo(fit)
write.tree(ccrC_dist_phylo,file = "ccrC_dist_tree_QC.nwk")

ccrA_tbl_2 = ccrA_tbl
colnames(ccrA_tbl_2)[8] = "ccrX_group"
ccrB_tbl_2 = ccrB_tbl
colnames(ccrB_tbl_2)[8] = "ccrX_group"
ccrC_tbl_2 = ccrC_tbl
colnames(ccrC_tbl_2)[8] = "ccrX_group"
ccr_uniq_all = rbind(ccrA_tbl_2,ccrB_tbl_2,ccrC_tbl_2)

ccr_uniq_all$source = "RefSeq"
ccr_uniq_all$source[which(ccr_uniq_all$GCF_ID=="NA_NA")] = "IWG_reference"


ccrA_tbl$source = "RefSeq"
ccrA_tbl$source[which(ccrA_tbl$GCF_ID=="NA_NA")] = "IWG_reference"
freq_tbl = as.data.frame(table(ccrA_tbl$uniq_ID,ccrA_tbl$source))

ccrA_uniq_tbl$IWG_references = unlist(lapply(as.vector(ccrA_uniq_tbl$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
ccrA_uniq_tbl$IWG_reference = 1
ccrA_uniq_tbl$IWG_reference[which(ccrA_uniq_tbl$IWG_references == 0)] = 0

write.table(ccrA_uniq_tbl,"ccrA_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)


ccrB_tbl$source = "RefSeq"
ccrB_tbl$source[which(ccrB_tbl$GCF_ID=="NA_NA")] = "IWG_reference"
freq_tbl = as.data.frame(table(ccrB_tbl$uniq_ID,ccrB_tbl$source))

ccrB_uniq_tbl$IWG_references = unlist(lapply(as.vector(ccrB_uniq_tbl$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
ccrB_uniq_tbl$IWG_reference = 1
ccrB_uniq_tbl$IWG_reference[which(ccrB_uniq_tbl$IWG_references == 0)] = 0

write.table(ccrB_uniq_tbl,"ccrB_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)


ccrC_tbl$source = "RefSeq"
ccrC_tbl$source[which(ccrC_tbl$GCF_ID=="NA_NA")] = "IWG_reference"
freq_tbl = as.data.frame(table(ccrC_tbl$uniq_ID,ccrC_tbl$source))

ccrC_uniq_tbl$IWG_references = unlist(lapply(as.vector(ccrC_uniq_tbl$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
ccrC_uniq_tbl$IWG_reference = 1
ccrC_uniq_tbl$IWG_reference[which(ccrC_uniq_tbl$IWG_references == 0)] = 0

write.table(ccrC_uniq_tbl,"ccrC_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)



