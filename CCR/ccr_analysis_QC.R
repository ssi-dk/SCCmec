library(ape)
library(phytools)
library(alluvial)
setwd("/Volumes/data/MPV/projects/SCCmec/CCR")
setwd("E:/Github/SCCmec/CCR")

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

tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))

aln = read.dna("protein_blast/ccr_all_uniq_mafft.fasta","fasta")
aln = read.dna("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_all_uniq_mafft.fasta?raw=true","fasta")

dist_mat = read.table("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_muscle_percent_distance.txt?raw=true",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)





# ccrA #

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
ccr_groups = cutree(fit,h=18)
rect.hclust(fit,h=18)
tbl$ccrA_group = NA
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrA_group[which(tbl$uniq_fasta_ID==name)] = group
}
tbl$ccr_allotype = paste0(as.vector(tbl$ccr_haplotype),'n',as.vector(tbl$ccrA_group))



IWG_ccrA_tbl = tbl[grep('ccrA',tbl$fasta_ID),]
table(as.vector(IWG_ccrA_tbl$fasta_ID),as.vector(IWG_ccrA_tbl$ccrA_group))
sort(unique(ccr_groups))

ccrA_tbl = tbl[which(tbl$ccr_haplotype=="ccrA"),]
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

write.table(ccrA_tbl,"ccrA_allotype_table_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrA_uniq_tbl,"ccrA_allotype_uniq_table_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$allotype_color,' ',ccrA_uniq_tbl$ccr_allotype)
writeLines(to_print,con = "ccrA_allotype_colors_QC_18.txt")

ccrA_dist_phylo = as.phylo(fit)
write.tree(ccrA_dist_phylo,file = "ccrA_dist_tree_QC.nwk")
