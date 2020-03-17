library(ape)
library(phytools)
library(alluvial)
setwd("/Volumes/data/MPV/projects/SCCmec/CCR")

sp_count_table = read.table("/Volumes/data/DB/refseq/Staphylococcus_species_counts_191128.txt",sep = "\t")
sp_count_vec = sp_count_table$V2
names(sp_count_vec) = sp_count_table$V1

#### ccr all ####

tbl = read.table("protein_blast/ccr_all_with_IWG_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")

aln = read.dna("protein_blast/ccr_all_uniq_mafft.fasta","fasta")

dist_mat = read.table("protein_blast/ccr_muscle_percent_distance.txt",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)
dist_obj = as.dist(dist_mat)

plot(sort(rowSums(dist_mat)))

fit = hclust(dist_obj,method = "complete")

#### Haplotyes ####

plot(fit)
ccr_groups = cutree(fit,h=49)
rect.hclust(fit,h=49)

tbl$ccr_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_group[which(tbl$uniq_fasta_ID==name)] = group
}

tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_haplotype = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_haplotype[which(tbl$ccr_group==5)] = 'ccrB'
tbl$ccr_haplotype[which(tbl$ccr_group==6)] = 'ccrA'
tbl$ccr_haplotype[which(tbl$ccr_group==7)] = 'ccrC'
table(tbl$ccr_haplotype)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_haplotype = factor(tbl$ccr_haplotype)


haplotype_colors = RColorBrewer::brewer.pal(7,"Set1")
ccr_haplotype_vec = levels(tbl$ccr_haplotype)
tbl$haplotype_color = NA

for (i in 1:length(haplotype_colors)) {
  col = haplotype_colors[i]
  ccr_type = ccr_haplotype_vec[i]
  tbl$haplotype_color[which(tbl$ccr_haplotype==ccr_type)] = col
}

uniq_tbl = tbl[which(!duplicated(as.vector(tbl$uniq_fasta_ID))),]

#write.table(tbl,"protein_blast/ccr_haplotype_table.txt",sep = "\t",quote = F,row.names=F)
#write.table(uniq_tbl,"protein_blast/ccr_haplotype_uniq_table.txt",sep = "\t",quote = F,row.names=F)
# to_print = paste0(uniq_tbl$uniq_fasta_ID,' ',uniq_tbl$haplotype_color,' ',uniq_tbl$ccr_haplotype)
# writeLines(to_print,con = "protein_blast/col_test.txt")


#### Separate into 7 haplotype and align/analyze each haplotype ####

tree = read.newick("protein_blast/ccr_haplotype_uniq_fastas/ccrA_muscle_fasttree.nwk")

as.hclust(tree)


# ccrA #

dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrA_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
aln_dist = as.dist(dist_mat)
# aln = read.dna("protein_blast/ccr_haplotype_uniq_fastas/ccrA_aln.fasta","fasta")
# aln_dist = dist.dna(aln)
# aln_dist_mat = as.matrix(aln_dist)
# aln_dist_mat[1:10,1:10]

fit = hclust(aln_dist)

plot(fit)
ccr_groups = cutree(fit,h=25)
rect.hclust(fit,h=25)
tbl$ccrA_group = NA
new_ccr_count = 1
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
allotype_colors = RColorBrewer::brewer.pal(11,"Paired")
#allotype_colors = colorRamps::primary.colors(16)
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

write.table(ccrA_tbl,"protein_blast/ccrA_allotype_table.txt",sep = "\t",quote = F,row.names=F)
write.table(ccrA_uniq_tbl,"protein_blast/ccrA_allotype_uniq_table.txt",sep = "\t",quote = F,row.names=F)
to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$allotype_color,' ',ccrA_uniq_tbl$ccr_allotype)
#writeLines(to_print,con = "protein_blast/ccrA_allotype_colors.txt")

ccrA_dist_phylo = as.phylo(fit)
write.tree(ccrA_dist_phylo,file = "ccrA_dist_tree_2.nwk")

# ccrB #

dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrB_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
aln_dist = as.dist(dist_mat)

# aln = read.dna("protein_blast/ccr_haplotype_uniq_fastas/ccrB_aln.fasta","fasta")
# aln_dist = dist.dna(aln)
# aln_dist_mat = as.matrix(aln_dist)
#aln_dist_mat[1:10,1:10]

fit = hclust(aln_dist)
plot(fit)
ccr_groups = cutree(fit,h=20)
rect.hclust(fit,h=20)
tbl$ccrB_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrB_group[which(tbl$uniq_fasta_ID==name)] = group
  
}

IWG_ccrB_tbl = tbl[grep('ccrB',tbl$fasta_ID),]
table(as.vector(IWG_ccrB_tbl$fasta_ID),as.vector(IWG_ccrB_tbl$ccrB_group))
sort(unique(ccr_groups))


dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrC_pairwise_sim.txt",sep = "\t",header=T, row.names=1))
aln_dist = as.dist(dist_mat)

# aln = read.dna("protein_blast/ccr_haplotype_uniq_fastas/ccrC_aln.fasta","fasta")
# aln_dist = dist.dna(aln)
# aln_dist_mat = as.matrix(aln_dist)
#aln_dist_mat[1:10,1:10]

fit = hclust(aln_dist)
plot(fit)
ccr_groups = cutree(fit,h=22)
rect.hclust(fit,h=22)
tbl$ccrC_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccrC_group[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrC_tbl = tbl[grep('ccrC',tbl$fasta_ID),]
table(as.vector(IWG_ccrC_tbl$fasta_ID),as.vector(IWG_ccrC_tbl$ccrC_group))
sort(unique(ccr_groups))


tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

IWG_tbl = tbl[grep('_IWG',tbl$fasta_ID),]
table(as.vector(IWG_tbl$fasta_ID),as.vector(IWG_tbl$ccr_group))
sort(unique(ccr_groups))

tbl$ccr_haplotype = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_haplotype[which(tbl$ccr_group==5)] = 'ccrB'
tbl$ccr_haplotype[which(tbl$ccr_group==6)] = 'ccrA'
tbl$ccr_haplotype[which(tbl$ccr_group==7)] = 'ccrC'
table(tbl$ccr_haplotype)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_haplotype = factor(tbl$ccr_haplotype)


haplotype_colors = RColorBrewer::brewer.pal(7,"Set1")
ccr_haplotype_vec = levels(tbl$ccr_haplotype)
tbl$haplotype_color = NA

for (i in 1:length(haplotype_colors)) {
  col = haplotype_colors[i]
  ccr_type = ccr_haplotype_vec[i]
  tbl$haplotype_color[which(tbl$ccr_haplotype==ccr_type)] = col
}

uniq_tbl = tbl[which(!duplicated(as.vector(tbl$uniq_fasta_ID))),]





tbl$ID[which(!tbl$ccr_group %in% c(5,6,7))]




test_idx = which(rownames(dist_mat) %in% c("seq_103|9","seq_158|6","seq_518|2"))

dist_mat[test_idx,test_idx]

ccrA_table = read.table("ccrA_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrA_uniq_table = read.table("ccrA_genes_uniq.txt",sep = "\t")
colnames(ccrA_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrA_aln = read.dna("ccrA_genes_uniq_with_IWG_muscle.fasta","fasta")

ccrA_dist = dist.dna(ccrA_aln,model = "raw")

ccrA_fit = hclust(ccrA_dist,"complete")
plot(ccrA_fit)
ccrA_groups = cutree(ccrA_fit,h=0.18)
rect.hclust(ccrA_fit,h=0.18)

#ccrA_groups = cutree(ccrA_fit,k=5)
ccrA_table$cluster = NA
ccrA_table$cluster_gene = NA
ccrA_uniq_table$cluster = NA
ccrA_uniq_table$cluster_gene = NA
ccrA_IWG_table = data.frame('ID' = labels(ccrA_aln)[grep('_IWG',labels(ccrA_aln))], 'cluster'=NA, 'cluster_gene'=NA)
for (i in unique(ccrA_groups)) {
  print(paste0(i))
  seqs = names(ccrA_groups)[which(ccrA_groups==i)]
  IWG_IDs = seqs[grep('IWG',seqs)]
  ccrA_IWG_table$cluster[which(ccrA_IWG_table$ID %in% IWG_IDs)] = i
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrA_uniq_table$fasta_ID[which(ccrA_uniq_table$seq_ID %in% seqs)])
  ccrA_types = as.vector(ccrA_table$ccr_gene[which(ccrA_table$fasta_ID %in% IDs)])
  ccrA_table$cluster[which(ccrA_table$fasta_ID %in% IDs)] = i
  ccrA_uniq_table$cluster[which(ccrA_uniq_table$seq_ID %in% seqs)] = i
  pidents = as.vector(ccrA_table$pident[which(ccrA_table$fasta_ID %in% IDs)])
  print(table(ccrA_types))
  print(mean(pidents))
}

sort(unique(ccrA_groups[ccrA_fit$order]))
ccrA_cluster_names = c("ccrA4","ccrA12","ccrA9","ccrA2","ccrA10","ccrA11","ccrA13","ccrA14","ccrA15","ccrA1","ccrA7","ccrA16","ccrA8","ccrA3","ccrA5","ccrA6")
ccrA_cluster_names = ccrA_cluster_names[order(unique(ccrA_groups[ccrA_fit$order]))]
names(ccrA_cluster_names) = 1:length(ccrA_cluster_names)

for (i in 1:length(ccrA_cluster_names)) {
  name = ccrA_cluster_names[i]
  ccrA_table$cluster_gene[which(ccrA_table$cluster==i)] = name
  ccrA_uniq_table$cluster_gene[which(ccrA_uniq_table$cluster==i)] = name
  ccrA_IWG_table$cluster_gene[which(ccrA_IWG_table$cluster==i)] = name
}

sort(unique(ccrA_groups[ccrA_fit$order]))

#### ccrB ####

ccrB_table = read.table("ccrB_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrB_uniq_table = read.table("ccrB_genes_uniq.txt",sep = "\t")
colnames(ccrB_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrB_aln = read.dna("ccrB_genes_uniq_with_IWG_muscle.fasta","fasta")

ccrB_dist = dist.dna(ccrB_aln,model = "raw")

ccrB_fit = hclust(ccrB_dist,method = "complete")
plot(ccrB_fit)
ccrB_groups = cutree(ccrB_fit,h=0.18)
rect.hclust(ccrB_fit,h=0.18)

ccrB_table$cluster = NA
ccrB_table$cluster_gene = NA
ccrB_uniq_table$cluster = NA
ccrB_uniq_table$cluster_gene = NA
ccrB_IWG_table = data.frame('ID' = labels(ccrB_aln)[grep('_IWG',labels(ccrB_aln))], 'cluster'=NA, 'cluster_gene'=NA)

for (i in unique(ccrB_groups)) {
  print(paste0(i))
  seqs = names(ccrB_groups)[which(ccrB_groups==i)]
  IWG_IDs = seqs[grep('IWG',seqs)]
  ccrB_IWG_table$cluster[which(ccrB_IWG_table$ID %in% IWG_IDs)] = i
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrB_uniq_table$fasta_ID[which(ccrB_uniq_table$seq_ID %in% seqs)])
  ccrB_types = as.vector(ccrB_table$ccr_gene[which(ccrB_table$fasta_ID %in% IDs)])
  ccrB_table$cluster[which(ccrB_table$fasta_ID %in% IDs)] = i
  ccrB_uniq_table$cluster[which(ccrB_uniq_table$seq_ID %in% seqs)] = i
  pidents = as.vector(ccrB_table$pident[which(ccrB_table$fasta_ID %in% IDs)])
  print(table(ccrB_types))
  print(mean(pidents))
}

unique(ccrB_groups[ccrB_fit$order])

ccrB_cluster_names = c("ccrB12","ccrB4","ccrB2","ccrB10","ccrB6","ccrB1","ccrB8","ccrB7","ccrB5","ccrB9","ccrB13","ccrB11","ccrB3")
ccrB_cluster_names = ccrB_cluster_names[order(unique(ccrB_groups[ccrB_fit$order]))]
names(ccrB_cluster_names) = 1:length(ccrB_cluster_names)

for (i in 1:length(ccrB_cluster_names)) {
  name = ccrB_cluster_names[i]
  ccrB_table$cluster_gene[which(ccrB_table$cluster==i)] = name
  ccrB_uniq_table$cluster_gene[which(ccrB_uniq_table$cluster==i)] = name
  ccrB_IWG_table$cluster_gene[which(ccrB_IWG_table$cluster==i)] = name
}

#### ccrC ####

ccrC_table = read.table("ccrC_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrC_uniq_table = read.table("ccrC_genes_uniq.txt",sep = "\t")
colnames(ccrC_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrC_aln = read.dna("ccrC_genes_uniq_with_IWG_muscle.fasta","fasta")

ccrC_dist = dist.dna(ccrC_aln,model = "raw")

ccrC_fit = hclust(ccrC_dist)
plot(ccrC_fit)
ccrC_groups = cutree(ccrC_fit,h=0.18)
rect.hclust(ccrC_fit,h=0.18)

#ccrC_groups = cutree(ccrC_fit,k=5)

ccrC_table$cluster = NA
ccrC_table$cluster_gene = NA
ccrC_uniq_table$cluster = NA
ccrC_uniq_table$cluster_gene = NA
ccrC_IWG_table = data.frame('ID' = labels(ccrC_aln)[grep('_IWG',labels(ccrC_aln))], 'cluster'=NA, 'cluster_gene'=NA)

for (i in unique(ccrC_groups)) {
  print(paste0(i))
  seqs = names(ccrC_groups)[which(ccrC_groups==i)]
  IWG_IDs = seqs[grep('IWG',seqs)]
  ccrC_IWG_table$cluster[which(ccrC_IWG_table$ID %in% IWG_IDs)] = i
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrC_uniq_table$fasta_ID[which(ccrC_uniq_table$seq_ID %in% seqs)])
  ccrC_types = as.vector(ccrC_table$ccr_gene[which(ccrC_table$fasta_ID %in% IDs)])
  pidents = as.vector(ccrC_table$pident[which(ccrC_table$fasta_ID %in% IDs)])
  ccrC_table$cluster[which(ccrC_table$fasta_ID %in% IDs)] = i
  ccrC_uniq_table$cluster[which(ccrC_uniq_table$seq_ID %in% seqs)] = i
  print(table(ccrC_types))
  print(mean(pidents))
}

unique(ccrC_groups[ccrC_fit$order])


ccrC_cluster_names = c("ccrC2","ccrC6","ccrC5","ccrC1","ccrC4","ccrC3")
ccrC_cluster_names = ccrC_cluster_names[order(unique(ccrC_groups[ccrC_fit$order]))]
names(ccrC_cluster_names) = 1:length(ccrC_cluster_names)

for (i in 1:length(ccrC_cluster_names)) {
  name = ccrC_cluster_names[i]
  ccrC_table$cluster_gene[which(ccrC_table$cluster==i)] = name
  ccrC_uniq_table$cluster_gene[which(ccrC_uniq_table$cluster==i)] = name
  ccrC_IWG_table$cluster_gene[which(ccrC_IWG_table$cluster==i)] = name
}


sort(unique(ccrC_groups[ccrC_fit$order]))

#### ####

ccrA_gene_count = table(ccrA_table$ID)
ccrB_gene_count = table(ccrB_table$ID)
ccrC_gene_count = table(ccrC_table$ID)

ccrX_table = rbind(ccrA_table,ccrB_table,ccrC_table)


ccrX_table$ccrA_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrA_table$ID)==x))))
ccrX_table$ccrB_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrB_table$ID)==x))))
ccrX_table$ccrC_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrC_table$ID)==x))))

ccrX_table$ccrA_cluster = NA


single_AB_isolates = unique(as.vector(ccrX_table$ID[which(ccrX_table$ccrA_gene_count==1 & ccrX_table$ccrB_gene_count==1 & ccrX_table$ccrC_gene_count<2)]))
single_AB_table = ccrX_table[which(ccrX_table$ID %in% single_AB_isolates),]

single_AB_table_ccrA = ccrA_table[which(ccrA_table$ID %in% single_AB_isolates),]
single_AB_table_ccrA$ID = as.vector(single_AB_table_ccrA$ID)
single_AB_table_ccrA = single_AB_table_ccrA[order(single_AB_table_ccrA$ID),]

single_AB_table_ccrB = ccrB_table[which(ccrB_table$ID %in% single_AB_isolates),]
single_AB_table_ccrB$ID = as.vector(single_AB_table_ccrB$ID)
single_AB_table_ccrB = single_AB_table_ccrB[order(single_AB_table_ccrB$ID),]

single_AB_table_ccrC = ccrC_table[which(ccrC_table$ID %in% single_AB_isolates),]
single_AB_table_ccrC$ID = as.vector(single_AB_table_ccrC$ID)
single_AB_table_ccrC = single_AB_table_ccrC[order(single_AB_table_ccrC$ID),]

single_AB_table_combined = data.frame('ID'=single_AB_table_ccrA$ID,'species'=single_AB_table_ccrA$species,'ccrA_cluster'=single_AB_table_ccrA$cluster_gene,'ccrB_cluster'=single_AB_table_ccrB$cluster_gene)
single_AB_table$ccrC_cluster = 0
get_ccrC <- function(ID,ccrC_table) {
  if (ID %in% ccrC_table$ID) {
    ccrC_gene = ccrC_table$cluster_gene[which(ccrC_table$ID==ID)]
  } else {
    ccrC_gene = 0
  }
  return(ccrC_gene)
}

ccrC_vec = unlist(lapply(as.vector(single_AB_table_combined$ID), function(x) get_ccrC(x,single_AB_table_ccrC)))
single_AB_table_combined$ccrC_cluster = ccrC_vec
single_AB_table_combined$test_string = paste0(single_AB_table_combined$species,' ',single_AB_table_combined$ccrA_cluster,' ',single_AB_table_combined$ccrB_cluster,' ',single_AB_table_combined$ccrC_cluster)

single_AB_table_combined$freq = unlist(lapply(single_AB_table_combined$test_string, function(x) length(which(as.vector(single_AB_table_combined$test_string)==x))))

alluvial_table = single_AB_table_combined[!duplicated(single_AB_table_combined$test_string),]

alluvial_table$ccrA_cluster = factor(as.vector(alluvial_table$ccrA_cluster), levels =c("ccrA1","ccrA2","ccrA3","ccrA4","ccrA5","ccrA6","ccrA7","ccrA8","ccrA9","ccrA10","ccrA11","ccrA12","ccrA13","ccrA14"))
alluvial_table$ccrB_cluster = factor(as.vector(alluvial_table$ccrB_cluster), levels =c("ccrB1","ccrB2","ccrB3","ccrB4","ccrB5","ccrB6","ccrB7","ccrB8","ccrB9","ccrB10","ccrB11","ccrB12","ccrB13"))

species_count_sorted = sort(table(ccrX_table$species),decreasing = T)
sp_color_vec = c(RColorBrewer::brewer.pal(9,"Set1"),"#cccccc")
names(sp_color_vec) = c(names(sort(sp_count_vec,decreasing = T))[1:9],"Other")

get_col = function(sp,color_vec) {
  if(sp %in% names(color_vec)) {
    col = color_vec[which(names(color_vec)==sp)]
  } else {
    col = "#cccccc"
  }
}

#write.table(single_AB_table,"single_AB_table.txt",sep = "\t")

alluvial_table$species_color = unlist(lapply(alluvial_table$species, function(x) get_col(x,sp_color_vec)))

alluvial_table$logfreq = log((alluvial_table$freq+1),10)

alluvial_table_sub = alluvial_table[which(alluvial_table$freq>1),]

unique(as.vector(alluvial_table_sub$species))

alluvial(alluvial_table_sub[,c(3,4,5)], freq=alluvial_table_sub$logfreq,
         col = alluvial_table_sub$species_color,
         border = alluvial_table_sub$species_color,
         alpha = 0.6
         #hide = tit$Freq == 0,
         #cex = 0.7
)



alluvial(alluvial_table[,c(3,4)], freq=alluvial_table$logfreq,
         col = ifelse(alluvial_table$species == "Staphylococcus_aureus", "orange", "grey"),
         border = ifelse(alluvial_table$species == "Staphylococcus_aureus", "orange", "grey"),
         alpha = 0.6
         #hide = tit$Freq == 0,
         #cex = 0.7
)


alluvial_aureus = alluvial_table[which(alluvial_table$species=="Staphylococcus_aureus"),]

alluvial(alluvial_aureus[,c(3,4)], freq=alluvial_aureus$logfreq,
         col = ifelse(alluvial_aureus$species == "Staphylococcus_aureus", "orange", "grey"),
         border = ifelse(alluvial_aureus$species == "Staphylococcus_aureus", "orange", "grey")
         #hide = tit$Freq == 0,
         #cex = 0.7
)



alluvial(alluvial_table[,c(3,4)], freq=alluvial_table$logfreq,
         col = ifelse(alluvial_table$species == "Staphylococcus_aureus", "orange", "grey"),
         border = ifelse(alluvial_table$species == "Staphylococcus_aureus", "orange", "grey")
         #hide = tit$Freq == 0,
         #cex = 0.7
)


plot_colors <- function(color_vec) {
  df = data.frame('ID'=factor(names(color_vec),levels=names(color_vec)),'color'=color_vec,'count'=1)
  p <- plot_ly(data = df, type='bar',x=~ID,y=~count,color=~ID,colors=color_vec)
  print(df)
  print(df$ID)
  return(p)
}


#ggplotdata = 
  
plot_colors(sp_color_vec)

ggplot(data = )

grep('IWG',ccrA_uniq_table$fasta_ID)

ccrA_table[1:2,]





ccrA_dist_tree = as.phylo(ccrA_fit)
ccrB_dist_tree = as.phylo(ccrB_fit)
ccrC_dist_tree = as.phylo(ccrC_fit)

rename_function = function(name) {
  if (substr(name,1,4)=="seq_") {
    newname = strsplit(name,'|',fixed = T)[[1]][1]
  } else {
    newname = name
  }
  return(newname)
}

ccrA_dist_tree$tip.label = unlist(lapply(ccrA_dist_tree$tip.label, rename_function))
ccrB_dist_tree$tip.label = unlist(lapply(ccrB_dist_tree$tip.label, rename_function))
ccrC_dist_tree$tip.label = unlist(lapply(ccrC_dist_tree$tip.label, rename_function))

write.tree(ccrA_dist_tree,file = "ccrA_dist_tree.nwk")
write.tree(ccrB_dist_tree,file = "ccrB_dist_tree.nwk")
write.tree(ccrC_dist_tree,file = "ccrC_dist_tree.nwk")

ccrA_color_vec = c(RColorBrewer::brewer.pal(12,"Paired"),"#d9d9d9","#636363")
ccrB_color_vec = RColorBrewer::brewer.pal(12,"Paired")
ccrC_color_vec = RColorBrewer::brewer.pal(6,"Set1")

names(ccrA_color_vec) = c("ccrA1","ccrA2","ccrA3","ccrA4","ccrA5","ccrA6","ccrA7","ccrA8","ccrA9","ccrA10","ccrA11","ccrA12","ccrA13","ccrA14")
names(ccrB_color_vec) = c("ccrB1","ccrB2","ccrB3","ccrB4","ccrB5","ccrB6","ccrB7","ccrB8","ccrB9","ccrB10","ccrB11","ccrB12")
names(ccrC_color_vec) = c("ccrC1","ccrC2","ccrC3","ccrC4","ccrC5","ccrC6")

ccrA_table = ccrA_table[order(ccrA_table$fasta_ID),]
ccrA_uniq_table = ccrA_uniq_table[order(ccrA_uniq_table$fasta_ID),]
ccrB_table = ccrB_table[order(ccrB_table$fasta_ID),]
ccrB_uniq_table = ccrB_uniq_table[order(ccrB_uniq_table$fasta_ID),]
ccrC_table = ccrC_table[order(ccrC_table$fasta_ID),]
ccrC_uniq_table = ccrC_uniq_table[order(ccrC_uniq_table$fasta_ID),]

ccrA_uniq_table$species = ccrA_table$species
ccrB_uniq_table$species = ccrB_table$species
ccrC_uniq_table$species = ccrC_table$species

ccrA_species_table = as.data.frame.matrix(table(ccrA_uniq_table[,c(2,6)]))
ccrA_species_table_collapsed = ccrA_species_table[,which(colnames(ccrA_species_table) %in% names(sp_color_vec))]
ccrA_species_table_collapsed$Other = rowSums(ccrA_species_table[,which(!colnames(ccrA_species_table) %in% names(sp_color_vec))])

ccrB_species_table = as.data.frame.matrix(table(ccrB_uniq_table[,c(2,6)]))
ccrB_species_table_collapsed = ccrB_species_table[,which(colnames(ccrB_species_table) %in% names(sp_color_vec))]
ccrB_species_table_collapsed$Other = rowSums(ccrB_species_table[,which(!colnames(ccrB_species_table) %in% names(sp_color_vec))])

ccrC_species_table = as.data.frame.matrix(table(ccrC_uniq_table[,c(2,6)]))
ccrC_species_table_collapsed = ccrC_species_table[,which(colnames(ccrC_species_table) %in% names(sp_color_vec))]
ccrC_species_table_collapsed$Other = rowSums(ccrC_species_table[,which(!colnames(ccrC_species_table) %in% names(sp_color_vec))])

seq_species_count_list = list('ccrA'=ccrA_species_table_collapsed,'ccrB'=ccrB_species_table_collapsed,'ccrC'=ccrC_species_table_collapsed)

write.xlsx(seq_species_count_list,file = "seq_species_counts.xlsx",row.names =T)


ccrA_uniq_nondupl = ccrA_uniq_table[!duplicated(ccrA_uniq_table$seq_ID),]
ccrA_uniq_nondupl$ccr_color = unlist(lapply(ccrA_uniq_nondupl$cluster_gene, function(x) get_col(x,ccrA_color_vec)))

ccrB_uniq_nondupl = ccrB_uniq_table[!duplicated(ccrB_uniq_table$seq_ID),]
ccrB_uniq_nondupl$ccr_color = unlist(lapply(ccrB_uniq_nondupl$cluster_gene, function(x) get_col(x,ccrB_color_vec)))

ccrC_uniq_nondupl = ccrC_uniq_table[!duplicated(ccrC_uniq_table$seq_ID),]
ccrC_uniq_nondupl$ccr_color = unlist(lapply(ccrC_uniq_nondupl$cluster_gene, function(x) get_col(x,ccrC_color_vec)))

ccrA_IWG_table$species = 'IWG_reference'
colnames(ccrA_IWG_table)[1] = 'seq_ID'
print_ccrA = rbind(ccrA_uniq_nondupl[,c(2,5,6)],ccrA_IWG_table[,c(1,3,4)])
print_ccrA$sp_color = unlist(lapply(print_ccrA$species, function(x) get_col(x,sp_color_vec)))
print_ccrA$ccrA_color = unlist(lapply(print_ccrA$cluster_gene, function(x) get_col(x,ccrA_color_vec)))


ccrB_IWG_table$species = 'IWG_reference'
colnames(ccrB_IWG_table)[1] = 'seq_ID'
print_ccrB = rbind(ccrB_uniq_nondupl[,c(2,5,6)],ccrB_IWG_table[,c(1,3,4)])
print_ccrB$sp_color = unlist(lapply(print_ccrB$species, function(x) get_col(x,sp_color_vec)))
print_ccrB$ccrB_color = unlist(lapply(print_ccrB$cluster_gene, function(x) get_col(x,ccrB_color_vec)))


ccrC_IWG_table$species = 'IWG_reference'
colnames(ccrC_IWG_table)[1] = 'seq_ID'
print_ccrC = rbind(ccrC_uniq_nondupl[,c(2,5,6)],ccrC_IWG_table[,c(1,3,4)])
print_ccrC$sp_color = unlist(lapply(print_ccrC$species, function(x) get_col(x,sp_color_vec)))
print_ccrC$ccrC_color = unlist(lapply(print_ccrC$cluster_gene, function(x) get_col(x,ccrC_color_vec)))


ccr_xl_list = list('ccrA'=print_ccrA,'ccrB'=print_ccrB,'ccrC'=print_ccrC)

write.xlsx(ccr_xl_list,file = "/Volumes/data/MPV/projects/SCCmec/CCR/ccr_phylogeny_tables.xlsx")




##### Find IWG matches #####

ccrA_IWG_test = read.table("ccrA_IWG_uniq.txt",sep = "\t")
colnames(ccrA_IWG_test) = c("fasta_ID","uniq_seq","seq","seq_name","seq_count")

IWG_rows = grep('IWG',ccrA_IWG_test$fasta_ID) 
ccrA_IWG_uniq_seqs = unique(as.vector(ccrA_IWG_test$uniq_seq[IWG_rows]))
ccrA_match_IDs = as.vector(ccrA_IWG_test$fasta_ID[which(ccrA_IWG_test$uniq_seq %in% ccrA_IWG_uniq_seqs)])


ccrB_IWG_test = read.table("ccrB_IWG_uniq.txt",sep = "\t")
colnames(ccrB_IWG_test) = c("fasta_ID","uniq_seq","seq","seq_name","seq_count")

IWG_rows = grep('IWG',ccrB_IWG_test$fasta_ID) 
ccrB_IWG_uniq_seqs = unique(as.vector(ccrB_IWG_test$uniq_seq[IWG_rows]))
ccrB_match_IDs = as.vector(ccrB_IWG_test$fasta_ID[which(ccrB_IWG_test$uniq_seq %in% ccrB_IWG_uniq_seqs)])


ccrC_IWG_test = read.table("ccrC_IWG_uniq.txt",sep = "\t")
colnames(ccrC_IWG_test) = c("fasta_ID","uniq_seq","seq","seq_name","seq_count")

IWG_rows = grep('IWG',ccrC_IWG_test$fasta_ID) 
ccrC_IWG_uniq_seqs = unique(as.vector(ccrC_IWG_test$uniq_seq[IWG_rows]))
ccrC_match_IDs = as.vector(ccrC_IWG_test$fasta_ID[which(ccrC_IWG_test$uniq_seq %in% ccrC_IWG_uniq_seqs)])



which(!ccrA_match_IDs %in% ccrA_dist_tree$tip.label)
which(!ccrB_match_IDs %in% ccrB_dist_tree$tip.label)
which(!ccrC_match_IDs %in% ccrC_dist_tree$tip.label)


ccrA_refs = ccrA_match_IDs[-grep('IWG',ccrA_match_IDs)]
ccrB_refs = ccrB_match_IDs[-grep('IWG',ccrB_match_IDs)]
ccrC_refs = ccrC_match_IDs[-grep('IWG',ccrC_match_IDs)]

ccrA_r = unlist(lapply(ccrA_refs,function(x) strsplit(x,'|',fixed=T)[[1]][1]))
ccrB_r = unlist(lapply(ccrB_refs,function(x) strsplit(x,'|',fixed=T)[[1]][1]))
ccrC_r = unlist(lapply(ccrC_refs,function(x) strsplit(x,'|',fixed=T)[[1]][1]))

print_ccrA$IWG_ref = 0
print_ccrA$IWG_ref[which(print_ccrA$seq_ID %in% ccrA_r)] = 1
print_ccrB$IWG_ref = 0
print_ccrB$IWG_ref[which(print_ccrB$seq_ID %in% ccrB_r)] = 1
print_ccrC$IWG_ref = 0
print_ccrC$IWG_ref[which(print_ccrC$seq_ID %in% ccrC_r)] = 1

ccrA_dist_tree_noIWG = ccrA_dist_tree
ccrA_dist_tree_noIWG = drop.tip(ccrA_dist_tree_noIWG,ccrA_dist_tree_noIWG$tip.label[grep('IWG',ccrA_dist_tree_noIWG$tip.label)])
ccrB_dist_tree_noIWG = ccrB_dist_tree
ccrB_dist_tree_noIWG = drop.tip(ccrB_dist_tree_noIWG,ccrB_dist_tree_noIWG$tip.label[grep('IWG',ccrB_dist_tree_noIWG$tip.label)])
ccrC_dist_tree_noIWG = ccrC_dist_tree
ccrC_dist_tree_noIWG = drop.tip(ccrC_dist_tree_noIWG,ccrC_dist_tree_noIWG$tip.label[grep('IWG',ccrC_dist_tree_noIWG$tip.label)])


print_ccrA_noIWG = print_ccrA
print_ccrA_noIWG = print_ccrA_noIWG[-grep('IWG',print_ccrA_noIWG$seq_ID),]
print_ccrB_noIWG = print_ccrB
print_ccrB_noIWG = print_ccrB_noIWG[-grep('IWG',print_ccrB_noIWG$seq_ID),]
print_ccrC_noIWG = print_ccrC
print_ccrC_noIWG = print_ccrC_noIWG[-grep('IWG',print_ccrC_noIWG$seq_ID),]

ccr_xl_list = list('ccrA'=print_ccrA_noIWG,'ccrB'=print_ccrB_noIWG,'ccrC'=print_ccrC_noIWG)

write.xlsx(ccr_xl_list,file = "/Volumes/data/MPV/projects/SCCmec/CCR/ccr_phylogeny_tables_IWG_removed.xlsx")

write.tree(ccrA_dist_tree_noIWG,"/Volumes/data/MPV/projects/SCCmec/CCR/ccrA_dist_tree_noIWG.nwk")
write.tree(ccrB_dist_tree_noIWG,"/Volumes/data/MPV/projects/SCCmec/CCR/ccrB_dist_tree_noIWG.nwk")
write.tree(ccrC_dist_tree_noIWG,"/Volumes/data/MPV/projects/SCCmec/CCR/ccrC_dist_tree_noIWG.nwk")




#### mafft distance ####

#### ccrA ####

ccrA_table = read.table("ccrA_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrA_uniq_table = read.table("ccrA_genes_uniq.txt",sep = "\t")
colnames(ccrA_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrA_aln = read.dna("ccrA_genes_uniq_with_IWG_aln.fasta","fasta")

ccrA_dist = dist.dna(ccrA_aln,model = "raw")

ccrA_fit = hclust(ccrA_dist,"complete")
plot(ccrA_fit)
ccrA_groups = cutree(ccrA_fit,h=0.18)
rect.hclust(ccrA_fit,h=0.18)

#ccrA_groups = cutree(ccrA_fit,k=5)
ccrA_table$cluster = NA
ccrA_table$cluster_gene = NA
ccrA_uniq_table$cluster = NA
ccrA_uniq_table$cluster_gene = NA
for (i in unique(ccrA_groups)) {
  print(paste0(i))
  seqs = names(ccrA_groups)[which(ccrA_groups==i)]
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrA_uniq_table$fasta_ID[which(ccrA_uniq_table$seq_ID %in% seqs)])
  IWG_IDs = IDs[which(!grep('IWG',IDs))]
  ccrA_types = as.vector(ccrA_table$ccr_gene[which(ccrA_table$fasta_ID %in% IDs)])
  ccrA_table$cluster[which(ccrA_table$fasta_ID %in% IDs)] = i
  ccrA_uniq_table$cluster[which(ccrA_uniq_table$seq_ID %in% seqs)] = i
  pidents = as.vector(ccrA_table$pident[which(ccrA_table$fasta_ID %in% IDs)])
  print(table(ccrA_types))
  print(mean(pidents))
}

unique(ccrA_groups[ccrA_fit$order])
ccrA_cluster_names = c("ccrA4a","ccrA4b","ccrA4b","ccrN1","ccrAN2","ccrAN3","ccrAN4","ccrA1","ccrAN5","ccrAN6","ccrAN7","ccrA3","ccrA5","ccrA2","ccrAN8","ccrAN9")
ccrA_cluster_names = ccrA_cluster_names[order(unique(ccrA_groups[ccrA_fit$order]))]
names(ccrA_cluster_names) = 1:length(ccrA_cluster_names)

for (i in 1:length(ccrA_cluster_names)) {
  name = ccrA_cluster_names[i]
  ccrA_table$cluster_gene[which(ccrA_table$cluster==i)] = name
  ccrA_uniq_table$cluster_gene[which(ccrA_uniq_table$cluster==i)] = name
}

sort(unique(ccrA_groups[ccrA_fit$order]))

#### ccrB ####

ccrB_table = read.table("ccrB_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrB_uniq_table = read.table("ccrB_genes_uniq.txt",sep = "\t")
colnames(ccrB_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrB_aln = read.dna("ccrB_genes_uniq_with_IWG_aln.fasta","fasta")

ccrB_dist = dist.dna(ccrB_aln,model = "raw")

ccrB_fit = hclust(ccrB_dist,method = "complete")
plot(ccrB_fit)
ccrB_groups = cutree(ccrB_fit,h=0.18)
rect.hclust(ccrB_fit,h=0.18)

ccrB_table$cluster = NA
ccrB_table$cluster_gene = NA
ccrB_uniq_table$cluster = NA
ccrB_uniq_table$cluster_gene = NA

for (i in unique(ccrB_groups)) {
  print(paste0(i))
  seqs = names(ccrB_groups)[which(ccrB_groups==i)]
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrB_uniq_table$fasta_ID[which(ccrB_uniq_table$seq_ID %in% seqs)])
  ccrB_types = as.vector(ccrB_table$ccr_gene[which(ccrB_table$fasta_ID %in% IDs)])
  ccrB_table$cluster[which(ccrB_table$fasta_ID %in% IDs)] = i
  ccrB_uniq_table$cluster[which(ccrB_uniq_table$seq_ID %in% seqs)] = i
  pidents = as.vector(ccrB_table$pident[which(ccrB_table$fasta_ID %in% IDs)])
  print(table(ccrB_types))
  print(mean(pidents))
}

unique(ccrB_groups[ccrB_fit$order])

ccrB_cluster_names = c("ccrB4b","ccrB4a","ccrB2","ccrBN1","ccrBN2","ccrBN3","ccrB3b","ccrB3a","ccrBN4","ccrB1","ccrB6","ccrB8","ccrB7")
ccrB_cluster_names = ccrB_cluster_names[order(unique(ccrB_groups[ccrB_fit$order]))]
names(ccrB_cluster_names) = 1:length(ccrB_cluster_names)

for (i in 1:length(ccrB_cluster_names)) {
  name = ccrB_cluster_names[i]
  ccrB_table$cluster_gene[which(ccrB_table$cluster==i)] = name
  ccrB_uniq_table$cluster_gene[which(ccrB_uniq_table$cluster==i)] = name
}

#### ccrC ####

ccrC_table = read.table("ccrC_gene_table.txt",sep = "\t",row.names=NULL,header = T)

ccrC_uniq_table = read.table("ccrC_genes_uniq.txt",sep = "\t")
colnames(ccrC_uniq_table) = c("fasta_ID","seq_ID","seq")

ccrC_aln = read.dna("ccrC_genes_uniq_with_IWG_aln.fasta","fasta")

ccrC_dist = dist.dna(ccrC_aln,model = "raw")

ccrC_fit = hclust(ccrC_dist)
plot(ccrC_fit)
ccrC_groups = cutree(ccrC_fit,h=0.18)
rect.hclust(ccrC_fit,h=0.18)

#ccrC_groups = cutree(ccrC_fit,k=5)

ccrC_table$cluster = NA
ccrC_table$cluster_gene = NA
ccrC_uniq_table$cluster = NA
ccrC_uniq_table$cluster_gene = NA
for (i in unique(ccrC_groups)) {
  print(paste0(i))
  seqs = names(ccrC_groups)[which(ccrC_groups==i)]
  seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
  IDs = as.vector(ccrC_uniq_table$fasta_ID[which(ccrC_uniq_table$seq_ID %in% seqs)])
  ccrC_types = as.vector(ccrC_table$ccr_gene[which(ccrC_table$fasta_ID %in% IDs)])
  pidents = as.vector(ccrC_table$pident[which(ccrC_table$fasta_ID %in% IDs)])
  ccrC_table$cluster[which(ccrC_table$fasta_ID %in% IDs)] = i
  ccrC_uniq_table$cluster[which(ccrC_uniq_table$seq_ID %in% seqs)] = i
  print(table(ccrC_types))
  print(mean(pidents))
}

unique(ccrC_groups[ccrC_fit$order])


ccrC_cluster_names = c("ccrC2","ccrCN1","ccrCN2","ccrC1","ccrCN3","ccrCN4")
ccrC_cluster_names = ccrC_cluster_names[order(unique(ccrC_groups[ccrC_fit$order]))]
names(ccrC_cluster_names) = 1:length(ccrC_cluster_names)

for (i in 1:length(ccrC_cluster_names)) {
  name = ccrC_cluster_names[i]
  ccrC_table$cluster_gene[which(ccrC_table$cluster==i)] = name
  ccrC_uniq_table$cluster_gene[which(ccrC_uniq_table$cluster==i)] = name
}


sort(unique(ccrC_groups[ccrC_fit$order]))

#### ####

ccrA_gene_count = table(ccrA_table$ID)
ccrB_gene_count = table(ccrB_table$ID)
ccrC_gene_count = table(ccrC_table$ID)

ccrX_table = rbind(ccrA_table,ccrB_table,ccrC_table)

ccrX_table$ccrA_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrA_table$ID)==x))))
ccrX_table$ccrB_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrB_table$ID)==x))))
ccrX_table$ccrC_gene_count = unlist(lapply(ccrX_table$ID, function(x) length(which(as.vector(ccrC_table$ID)==x))))

ccrX_table$ccrA_cluster = NA


single_AB_isolates = unique(as.vector(ccrX_table$ID[which(ccrX_table$ccrA_gene_count==1 & ccrX_table$ccrB_gene_count==1 & ccrX_table$ccrC_gene_count==0)]))
single_AB_table = ccrX_table[which(ccrX_table$ID %in% single_AB_isolates),]

single_AB_table_ccrA = ccrA_table[which(ccrA_table$ID %in% single_AB_isolates),]
single_AB_table_ccrA$ID = as.vector(single_AB_table_ccrA$ID)
single_AB_table_ccrA = single_AB_table_ccrA[order(single_AB_table_ccrA$ID),]
single_AB_table_ccrB = ccrB_table[which(ccrB_table$ID %in% single_AB_isolates),]
single_AB_table_ccrB$ID = as.vector(single_AB_table_ccrB$ID)
single_AB_table_ccrB = single_AB_table_ccrB[order(single_AB_table_ccrB$ID),]

single_AB_table_combined = data.frame('ID'=single_AB_table_ccrA$ID,'species'=single_AB_table_ccrA$species,'ccrA_cluster'=single_AB_table_ccrA$cluster_gene,'ccrB_cluster'=single_AB_table_ccrB$cluster_gene)

single_AB_table_combined$test_string = paste0(single_AB_table_combined$species,' ',single_AB_table_combined$ccrA_cluster,' ',single_AB_table_combined$ccrB_cluster)

single_AB_table_combined$freq = unlist(lapply(single_AB_table_combined$test_string, function(x) length(which(as.vector(single_AB_table_combined$test_string)==x))))

alluvial_table = single_AB_table_combined[!duplicated(single_AB_table_combined$test_string),]

alluvial_table_sub = alluvial_table[which(alluvial_table$freq>1),]
alluvial_table_sub$logfreq = log(alluvial_table_sub$freq,10)

alluvial(alluvial_table_sub[,c(3,4)], freq=alluvial_table_sub$logfreq
         # col = ifelse(tit$Survived == "Yes", "orange", "grey"),
         #border = ifelse(tit$Survived == "Yes", "orange", "grey"),
         #hide = tit$Freq == 0,
         #cex = 0.7
)


