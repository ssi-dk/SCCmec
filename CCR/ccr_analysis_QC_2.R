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
tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/ccr_all_QC_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")


tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))
all_tbl = tbl

#aln = read.dna("protein_blast/ccr_all_uniq_mafft.fasta","fasta")
#aln = read.dna("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_all_uniq_mafft.fasta?raw=true","fasta")

dist_mat = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/ccr_all_QC_uniq_sim.txt",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)


dist_obj = as.dist(dist_mat)

fit = hclust(dist_obj)

plot(fit)


ccr_groups = cutree(fit,h=50)
rect.hclust(fit,h=50)
sort(unique(ccr_groups))

tbl$ccr_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_group[which(tbl$uniq_fasta_ID==name)] = group
}

tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_type = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_type[which(tbl$ccr_group==9)] = 'ccrA'
tbl$ccr_type[which(tbl$ccr_group==11)] = 'ccrB'
tbl$ccr_type[which(tbl$ccr_group==8)] = 'ccrC'
table(tbl$ccr_type)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_type = factor(tbl$ccr_type)

uniq_tbl = tbl[which(!duplicated(as.vector(tbl$uniq_fasta_ID))),]

all_tbl = tbl
all_tbl_uniq = uniq_tbl

ccr_dist_phylo = as.phylo(fit)
write.tree(ccr_dist_phylo,file = "protein_blast/ccr_combined/ccr_all_dist_tree_QC.nwk")


ccr_aln = read.FASTA("protein_blast/ccr_combined/ccr_nucleotide_uniq_muscle.fasta")

ccrA_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrA"),]
ccrB_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrB"),]
ccrC_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrC"),]
ccrA_uniq_fasta_print = paste0('>',as.vector(ccrA_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrA_uniq_tbl$seq))
writeLines(ccrA_uniq_fasta_print,"protein_blast/ccr_combined/ccrA_uniq_seqs.fasta")
ccrB_uniq_fasta_print = paste0('>',as.vector(ccrB_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrB_uniq_tbl$seq))
writeLines(ccrB_uniq_fasta_print,"protein_blast/ccr_combined/ccrB_uniq_seqs.fasta")
ccrC_uniq_fasta_print = paste0('>',as.vector(ccrC_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrC_uniq_tbl$seq))
writeLines(ccrC_uniq_fasta_print,"protein_blast/ccr_combined/ccrC_uniq_seqs.fasta")

# ccrA #


dist_mat = as.matrix(read.table("protein_blast/ccr_combined/ccrA_muscle_dist.txt",sep = "\t",header=T, row.names=1))
aln_dist = as.dist(dist_mat)


fit = hclust(aln_dist,method = "complete")

v = c()
for (i in 10:40) {
  g = cutree(fit,h=i)
  v = c(v,length(unique(g)))
}
v_A = v


plot(10:40,v)
plot(fit)
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

IWG_ccrA_tbl = tbl[grep('ccrA',tbl$fasta_ID),]
table(as.vector(IWG_ccrA_tbl$fasta_ID),as.vector(IWG_ccrA_tbl$ccrA_group))
sort(unique(ccr_groups))

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

#write.table(ccrA_tbl,"ccrA_allotype_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
#write.table(ccrA_uniq_tbl,"ccrA_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
#to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$allotype_color,' ',ccrA_uniq_tbl$ccr_allotype)
#writeLines(to_print,con = "ccrA_allotype_colors_QC_22.txt")

#ccrA_dist_phylo = as.phylo(fit)
#write.tree(ccrA_dist_phylo,file = "ccrA_dist_tree_QC.nwk")



#### ccrB ####

tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrB_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))
#tbl = tbl[which(as.vector(tbl$ccr_type)=="ccrB"),]

#dist_mat = as.matrix(read.table("protein_blast/ccr_haplotype_uniq_fastas/ccrB_pairwise_sim_2.txt",sep = "\t",header=T, row.names=1))
dist_mat = as.matrix(read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrB_all_uniq_sim.txt",sep = "\t",header=T, row.names=1))
#dist_mat_old = as.matrix(read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrB_all_uniq_sim_old.txt",sep = "\t",header=T, row.names=1))
dist_mat = as.matrix(read.table("fastas/ccrB_all_uniq_sim.txt",sep = "\t",header=T, row.names=1))

aln_dist = as.dist(dist_mat)


fit = hclust(aln_dist,method = "complete")

v = c()
for (i in 10:40) {
  g = cutree(fit,h=i)
  v = c(v,length(unique(g)))
}
v_B = v


plot(10:40,v)
plot(fit)
h = 22
ccr_groups = cutree(fit,h=h)
rect.hclust(fit,h=h)
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
allotype_colors = triplet_color_vec[1:19]
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
table(ccrB_tbl$ccr_allotype)

write.table(ccrB_tbl,"ccrB_allotype_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrB_uniq_tbl,"ccrB_allotype_uniq_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrB_uniq_tbl$uniq_fasta_ID,' ',ccrB_uniq_tbl$allotype_color,' ',ccrB_uniq_tbl$ccr_allotype)
writeLines(to_print,con = "ccrB_allotype_colors_QC_22.txt")

ccrB_dist_phylo = as.phylo(fit)
write.tree(ccrB_dist_phylo,file = "ccrB_dist_tree_QC.nwk")




#### ccrC ####

tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/fastas/ccrC_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
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




######### data tables ##########

get_top_n_from_variable = function(variable_factor,n) {
  variable_vector = as.vector(variable_factor)
  sorted_table = sort(table(variable_factor),decreasing = T)
  top_n_groups = names(sorted_table)[1:n]
  variable_vector[which(!variable_vector %in% top_n_groups)] = "Other"
  return_factor = factor(variable_vector, levels = c(top_n_groups,"Other"))
  return(return_factor)
}


ccrA_tbl_2 = ccrA_tbl
colnames(ccrA_tbl_2)[8] = "ccrX_group"
ccrB_tbl_2 = ccrB_tbl
colnames(ccrB_tbl_2)[8] = "ccrX_group"
ccrC_tbl_2 = ccrC_tbl
colnames(ccrC_tbl_2)[8] = "ccrX_group"
ccr_all = rbind(ccrA_tbl_2,ccrB_tbl_2,ccrC_tbl_2)

ccr_all$source = "RefSeq"
ccr_all$source[which(ccr_all$GCF_ID=="NA_NA")] = "IWG_reference"

ccr_uniq_all = ccr_all[!duplicated(ccr_all$seq),]

match_idx = match(as.vector(all_tbl_uniq$seq),as.vector(ccr_uniq_all$seq))

all_tbl_uniq$ccr_allotype = as.vector(ccr_uniq_all$ccr_allotype)[match_idx]
all_tbl_uniq$uniq_fasta_ID_split = as.vector(ccr_uniq_all$uniq_fasta_ID)[match_idx]
all_tbl_uniq$source = as.vector(ccr_uniq_all$source)[match_idx]
all_tbl_uniq$source[which(is.na(all_tbl_uniq$source))] = "RefSeq"

head(all_tbl_uniq)
head(ccr_uniq_all)

all_tbl$source = "RefSeq"
all_tbl$source[which(all_tbl$GCF_ID=="NA_NA")] = "IWG_reference"

freq_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$source))
all_tbl_uniq$IWG_references = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
all_tbl_uniq$IWG_reference = 1
all_tbl_uniq$IWG_reference[which(all_tbl_uniq$IWG_references == 0)] = 0

all_tbl_uniq$RefSeq_count = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="RefSeq")]))


match_idx = match(as.vector(all_tbl$seq),as.vector(ccr_all$seq))
all_tbl$ccr_allotype = as.vector(ccr_all$ccr_allotype)[match_idx]
all_tbl$uniq_fasta_ID_split = as.vector(ccr_all$uniq_fasta_ID)[match_idx]

all_tbl$species = unlist(lapply(as.vector(all_tbl$fasta_ID), function(x) strsplit(x,'__')[[1]][1]))
all_tbl$species[which(all_tbl$GCF_ID=="NA_NA")] = "IWG_reference"

all_tbl$species_simplified = get_top_n_from_variable(all_tbl$species,13)

sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))
sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species_simplified))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species_simplified =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))

write.table(all_tbl,"ccr_table_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)

write.table(all_tbl_uniq,"ccr_table_uniq_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)




ccrA_tbl$source = "RefSeq"
ccrA_tbl$source[which(ccrA_tbl$GCF_ID=="NA_NA")] = "IWG_reference"
freq_tbl = as.data.frame(table(ccrA_tbl$uniq_ID,ccrA_tbl$source))

ccrA_uniq_tbl$IWG_references = unlist(lapply(as.vector(ccrA_uniq_tbl$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
ccrA_uniq_tbl$IWG_reference = 1
ccrA_uniq_tbl$IWG_reference[which(ccrA_uniq_tbl$IWG_references == 0)] = 0

ccrA_tbl$species = unlist(lapply(as.vector(ccrA_tbl$fasta_ID), function(x) strsplit(x,'__')[[1]][1]))
ccrA_tbl$species[which(ccrA_tbl$GCF_ID=="NA_NA")] = "IWG_reference"

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


setup_color_lines = function(ID_vector,variable_factor,colors,legend_title) {
  colors = colors[1:length(levels(variable_factor))]
  lvls = levels(variable_factor)
  variable_vector = as.vector(variable_factor)
  color_vector = variable_vector
  for (i in 1:length(lvls)) {
    color_vector[which(variable_vector==lvls[i])] = colors[i]
  }
  itol_lines = paste0(ID_vector,"\t",color_vector,"\t",variable_vector)
  title_line = paste0("LEGEND_TITLE\t",legend_title)
  label_color_vec = c("LEGEND_COLORS",colors)
  label_name_vec = c("LEGEND_LABELS",lvls)
  label_shape_vec = c("LEGEND_SHAPES",rep(1,length(lvls)))
  label_color_line = paste0(label_color_vec, collapse = "\t")
  label_name_line = paste0(label_name_vec, collapse = "\t")
  label_shape_line = paste0(label_shape_vec, collapse = "\t")
  return(list('legend_title'=title_line,'label_shapes'=label_shape_line,'label_colors'=label_color_line,'label_names'=label_name_line,'data_lines'=itol_lines))
}

print_template = function(template,itol_lines,output_file,legend_title) {
  template[16] = paste0("DATASET_LABEL\t",legend_title)
  toprint = c(template,itol_lines$legend_title,itol_lines$label_shapes,itol_lines$label_colors,itol_lines$label_names,"DATA",itol_lines$data_lines)
  writeLines(toprint,output_file)
}


get_top_n_from_variable = function(variable_factor,n) {
  variable_vector = as.vector(variable_factor)
  sorted_table = sort(table(variable_factor),decreasing = T)
  top_n_groups = names(sorted_table)[1:n]
  variable_vector[which(!variable_vector %in% top_n_groups)] = "Other"
  return_factor = factor(variable_vector, levels = c(top_n_groups,"Other"))
  return(return_factor)
}


d_tbl = read.table("ccr_table_uniq_QC_22.txt", sep = "\t", header=T,comment.char = "",quote = "")

template = readLines("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/itol_template.txt")

itol_lines = setup_color_lines(d_tbl$uniq_fasta_ID,d_tbl$ccr_type,RColorBrewer::brewer.pal(11,"Paired"),legend_title = "ccr_type")
print_template(template,itol_lines,"ccr_type_all_colorstrip.txt",legend_title = "ccr_type")

itol_lines = setup_color_lines(ccrA_uniq_tbl$uniq_fasta_ID,ccrA_uniq_tbl$IWG_reference,c("#FFFFFF","#000000"),legend_title = "IWG_reference")
print_template(template,itol_lines,"ccrA_IWG_ref_colorstrip.txt",legend_title = "IWG_reference")

ccrB_uniq_tbl$IWG_reference = factor(ccrB_uniq_tbl$IWG_reference)
itol_lines = setup_color_lines(ccrB_uniq_tbl$uniq_fasta_ID,ccrB_uniq_tbl$IWG_reference,c("#FFFFFF","#000000"),legend_title = "IWG_reference")
print_template(template,itol_lines,"ccrB_IWG_ref_colorstrip.txt",legend_title = "IWG_reference")

ccrC_uniq_tbl$IWG_reference = factor(ccrC_uniq_tbl$IWG_reference)
itol_lines = setup_color_lines(ccrC_uniq_tbl$uniq_fasta_ID,ccrC_uniq_tbl$IWG_reference,c("#FFFFFF","#000000"),legend_title = "IWG_reference")
print_template(template,itol_lines,"ccrC_IWG_ref_colorstrip.txt",legend_title = "IWG_reference")





tbl = read.table("protein_blast/ccr_all_with_IWG_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/ccr_all_QC_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")


tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))
all_tbl = tbl

#aln = read.dna("protein_blast/ccr_all_uniq_mafft.fasta","fasta")
#aln = read.dna("https://github.com/ssi-dk/SCCmec/blob/master/CCR/ccr_all_uniq_mafft.fasta?raw=true","fasta")

dist_mat = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/CCR/ccr_all_QC_uniq_sim.txt",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)


dist_obj = as.dist(dist_mat)

fit = hclust(dist_obj)

plot(fit)


ccr_groups = cutree(fit,h=50)
rect.hclust(fit,h=50)
sort(unique(ccr_groups))

tbl$ccr_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_group[which(tbl$uniq_fasta_ID==name)] = group
}

tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_type = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_type[which(tbl$ccr_group==9)] = 'ccrA'
tbl$ccr_type[which(tbl$ccr_group==11)] = 'ccrB'
tbl$ccr_type[which(tbl$ccr_group==8)] = 'ccrC'
table(tbl$ccr_type)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_type = factor(tbl$ccr_type)



ccr_groups = cutree(fit,h=22)
rect.hclust(fit,h=22)
sort(unique(ccr_groups))

tbl$ccr_subgroup = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subgroup[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccr_tbl = tbl[grep('_IWG',tbl$fasta_ID),]
table(as.vector(IWG_ccr_tbl$fasta_ID),as.vector(IWG_ccr_tbl$ccr_subgroup))
table(tbl$ccr_subgroup)
sort(unique(ccr_groups))

tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_type = paste0('ccrN',as.vector(tbl$ccr_group))


tbl$ccr_type[which(tbl$ccr_group==9)] = 'ccrA'
tbl$ccr_type[which(tbl$ccr_group==11)] = 'ccrB'
tbl$ccr_type[which(tbl$ccr_group==8)] = 'ccrC'
table(tbl$ccr_type)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_type = factor(tbl$ccr_type)


