library(ape)
library(phytools)
library(alluvial)

#### functions ####

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

#### data ####

setwd("/Volumes/data/MPV/projects/SCCmec/CCR")
setwd("E:/Github/SCCmec/CCR")
#setwd("/Users/thej/Documents/GitHub/SCCmec/CCR")



sp_count_table = read.table("/Volumes/data/DB/refseq/Staphylococcus_species_counts_191128.txt",sep = "\t")
sp_count_table = read.table("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/Staphylococcus_species_counts_191128.txt",sep = "\t")
sp_count_vec = sp_count_table$V2
names(sp_count_vec) = sp_count_table$V1

triplet_color_vec = c("#61d2ff","#468bfa","#030bfc","#63ff7d","#23db64","#009133","#ff9999","#f03232","#c90000",
                      "#ffe194","#ffcc47","#e8a800","#f2a8ff","#e44dff","#9e00ba","#f8fc6a","#b8bd00","#786000",
                      "#c4c4c4","#545454","#000000")


#### blast outputs ####
## perform blastp
## blastp -query ../../references_ccrA/ccrA1_protein.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/protein_SCCmec_QC -out ccrA1_blast.tab -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
## blastp -query ../../references_ccrB/ccrB1_protein.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/protein_SCCmec_QC -out ccrB1_blast.tab -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
## blastp -query ../../references_ccrC/ccrC1_protein.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/protein_SCCmec_QC -out ccrC1_blast.tab -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"

blast_headers = c("qseqid", "sseqid", "qlen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
ccrA_blast_table = read.table("/Volumes/data/MPV/projects/SCCmec/CCR/protein_blast/ccr_combined/ccrA1_blast.tab", sep = "\t")
colnames(ccrA_blast_table) = blast_headers
ccrA_blast_table$slen = ccrA_blast_table$send-ccrA_blast_table$sstart+1
ccrA_blast_table$percent_len = ccrA_blast_table$slen/ccrA_blast_table$qlen*100
plot(ccrA_blast_table$pident,ccrA_blast_table$percent_len)

ccrB_blast_table = read.table("/Volumes/data/MPV/projects/SCCmec/CCR/protein_blast/ccr_combined/ccrB2_blast.tab", sep = "\t")
colnames(ccrB_blast_table) = blast_headers
ccrB_blast_table$slen = ccrB_blast_table$send-ccrB_blast_table$sstart+1
ccrB_blast_table$percent_len = ccrB_blast_table$slen/ccrB_blast_table$qlen*100
plot(ccrB_blast_table$pident,ccrB_blast_table$percent_len)

ccrC_blast_table = read.table("/Volumes/data/MPV/projects/SCCmec/CCR/protein_blast/ccr_combined/ccrC1_blast.tab", sep = "\t")
colnames(ccrC_blast_table) = blast_headers
ccrC_blast_table$slen = ccrC_blast_table$send-ccrC_blast_table$sstart+1
ccrC_blast_table$percent_len = ccrC_blast_table$slen/ccrC_blast_table$qlen*100
plot(ccrC_blast_table$pident,ccrC_blast_table$percent_len)


## concatenate ccr blast outputs:
## cat ccrA1_blast.tab ccrB2_blast.tab ccrC1_blast.tab > ccr_combined_blast.tab
## extract cds fastas for hits over 50% length and 50 pident
## python3 /srv/data/MPV/THEJ/scripts/extract_nucleotide_refseq_genus.py ccr_combined_blast.tab /srv/data/DB/refseq/cds/Staphylococcus ccr_combined_nucleotide.fasta 50 50
## add IWG ref sequences:
## cat ccr_combined_nucleotide.fasta ../../references_ccr/IWG_ccr*_all.fasta > ccr_nucleotide_with_IWG.fasta
## uniquefy:
## python3 /srv/data/MPV/THEJ/scripts/uniquify_generic.py ccr_nucleotide_with_IWG.fasta ccr_nucleotide_uniq.fasta ccr_nucleotide_uniq_table.txt
## run muscle alignment
## muscle -in ccr_nucleotide_uniq.fasta -out ccr_nucleotide_uniq_muscle.fasta
## sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccr_nucleotide_uniq.fasta -out ccr_nucleotide_uniq_muscle.fasta"
## calculate distance matrix
## python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccr_nucleotide_uniq_muscle.fasta > ccr_nucleotide_muscle_dist.txt





### ccr All

tbl = read.table("protein_blast/ccr_combined/ccr_nucleotide_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))
all_tbl = tbl

dist_mat = read.table("protein_blast/ccr_combined/ccr_nucleotide_muscle_dist.txt",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)
dist_obj = as.dist(dist_mat)

fit = hclust(dist_obj)

plot(fit)


ccr_groups = cutree(fit,h=45)
rect.hclust(fit,h=45)
sort(unique(ccr_groups))

tbl$ccr_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_group[which(tbl$uniq_fasta_ID==name)] = group
}

ccr_dist_phylo = as.phylo(fit)
write.tree(ccr_dist_phylo,file = "protein_blast/ccr_combined/ccr_all_dist_tree_QC.nwk")


tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_type = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_type[which(tbl$ccr_group==2)] = 'ccrA'
tbl$ccr_type[which(tbl$ccr_group==4)] = 'ccrB'
tbl$ccr_type[which(tbl$ccr_group==1)] = 'ccrC'
table(tbl$ccr_type)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_type = factor(tbl$ccr_type)



all_tbl = tbl
all_tbl$source = "RefSeq"
all_tbl$source[which(all_tbl$GCF_ID=="NA_NA")] = "IWG_reference"

all_tbl$species = unlist(lapply(as.vector(all_tbl$fasta_ID), function(x) strsplit(x,'__')[[1]][1]))
all_tbl$species[which(all_tbl$GCF_ID=="NA_NA")] = "IWG_reference"
all_tbl$species_simplified = get_top_n_from_variable(all_tbl$species,13)


all_tbl_uniq = all_tbl[which(!duplicated(as.vector(tbl$uniq_fasta_ID))),]

freq_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$source))
all_tbl_uniq$IWG_references = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
all_tbl_uniq$IWG_reference = 1
all_tbl_uniq$IWG_reference[which(all_tbl_uniq$IWG_references == 0)] = 0

all_tbl_uniq$RefSeq_count = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="RefSeq")]))

sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))
sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species_simplified))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species_simplified =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))

#write.table(all_tbl,"protein_blast/ccr_combined/ccr_table_QC.txt",sep = "\t",quote = FALSE,row.names=FALSE)

#write.table(all_tbl_uniq,"protein_blast/ccr_combined/ccr_table_uniq_QC.txt",sep = "\t",quote = FALSE,row.names=FALSE)



# ccrA_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrA"),]
# ccrB_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrB"),]
# ccrC_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrC"),]
# ccrA_uniq_fasta_print = paste0('>',as.vector(ccrA_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrA_uniq_tbl$seq))
# writeLines(ccrA_uniq_fasta_print,"protein_blast/ccr_combined/ccrA_uniq_seqs.fasta")
# ccrB_uniq_fasta_print = paste0('>',as.vector(ccrB_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrB_uniq_tbl$seq))
# writeLines(ccrB_uniq_fasta_print,"protein_blast/ccr_combined/ccrB_uniq_seqs.fasta")
# ccrC_uniq_fasta_print = paste0('>',as.vector(ccrC_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrC_uniq_tbl$seq))
# writeLines(ccrC_uniq_fasta_print,"protein_blast/ccr_combined/ccrC_uniq_seqs.fasta")

### align ccrA, ccrB, ccrC individually using muscle
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrA_uniq_seqs.fasta -out ccrA_uniq_muscle.fasta"
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrB_uniq_seqs.fasta -out ccrB_uniq_muscle.fasta"
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrC_uniq_seqs.fasta -out ccrC_uniq_muscle.fasta"
#
### do pairwise sim matrices
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrA_uniq_muscle.fasta > ccrA_muscle_dist.txt
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrB_uniq_muscle.fasta > ccrB_muscle_dist.txt
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrC_uniq_muscle.fasta > ccrC_muscle_dist.txt


all_tbl_uniq$ccr_subtype = "-"
all_tbl$ccr_subtype = "-"

#### ccrA ####

tbl = all_tbl[which(all_tbl$ccr_type=="ccrA"),]
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
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrA_tbl = tbl[grep('ccrA',tbl$fasta_ID),]
table(as.vector(IWG_ccrA_tbl$fasta_ID),as.vector(IWG_ccrA_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrAn',as.vector(tbl$ccr_subtype))
ccrA_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrA_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrA_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrA_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrA_tbl$ccr_subtype[idx] = t
  }
  ccrA_tbl$subtype_color[idx] = col
}

ccrA_uniq_tbl = ccrA_tbl[which(!duplicated(as.vector(ccrA_tbl$uniq_fasta_ID))),]

ccrA_sub_df = as.data.frame.matrix(table(as.vector(ccrA_tbl$uniq_fasta_ID),ccrA_tbl$ccr_subtype))
fasta_vec = c()
for (i in 1:ncol(ccrA_sub_df)) {
  top_ID = rownames(ccrA_sub_df)[order(ccrA_sub_df[,i], decreasing = TRUE)[1]]
  top_seq = ccrA_uniq_tbl$seq[which(ccrA_uniq_tbl$uniq_fasta_ID==top_ID)]
  printline = paste0('>',colnames(ccrA_sub_df)[i],'\n',top_seq)
  fasta_vec = c(fasta_vec,printline)
}
writeLines(fasta_vec,"protein_blast/subtype_fastas/ccrA_reps_from_protein_clustering.fasta")

# write.table(ccrA_tbl,"protein_blast/ccr_combined/ccrA_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# write.table(ccrA_uniq_tbl,"protein_blast/ccr_combined/ccrA_uniq_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$subtype_color,' ',ccrA_uniq_tbl$ccr_subtype)
# writeLines(to_print,con = "protein_blast/ccr_combined/ccrA_subtype_colors.txt")
# 
# ccrA_dist_phylo = as.phylo(fit)
# write.tree(ccrA_dist_phylo,file = "protein_blast/ccr_combined/ccrA_dist_tree_QC_22.nwk")



#### ccrB ####

# dist_mat = as.matrix(read.table("references_ccr/IWG_ccrB_dist.txt",sep = "\t",header=T, row.names=1))
# aln_dist = as.dist(dist_mat)


tbl = all_tbl[which(all_tbl$ccr_type=="ccrB"),]
dist_mat = as.matrix(read.table("protein_blast/ccr_combined/ccrB_muscle_dist.txt",sep = "\t",header=T, row.names=1))
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
ccr_groups = cutree(fit,h=18)
rect.hclust(fit,h=18)
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrB_tbl = tbl[grep('ccrB',tbl$fasta_ID),]
table(as.vector(IWG_ccrB_tbl$fasta_ID),as.vector(IWG_ccrB_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrBn',as.vector(tbl$ccr_subtype))
ccrB_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrB_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrB_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrB_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrB_tbl$ccr_subtype[idx] = t
  }
  ccrB_tbl$subtype_color[idx] = col
}

ccrB_uniq_tbl = ccrB_tbl[which(!duplicated(as.vector(ccrB_tbl$uniq_fasta_ID))),]
table(ccrB_tbl$ccr_subtype)

ccrB_sub_df = as.data.frame.matrix(table(as.vector(ccrB_tbl$uniq_fasta_ID),ccrB_tbl$ccr_subtype))
fasta_vec = c()
for (i in 1:ncol(ccrB_sub_df)) {
  top_ID = rownames(ccrB_sub_df)[order(ccrB_sub_df[,i], decreasing = TRUE)[1]]
  top_seq = ccrB_uniq_tbl$seq[which(ccrB_uniq_tbl$uniq_fasta_ID==top_ID)]
  printline = paste0('>',colnames(ccrB_sub_df)[i],'\n',top_seq)
  fasta_vec = c(fasta_vec,printline)
}
writeLines(fasta_vec,"protein_blast/subtype_fastas/ccrB_reps_from_protein_clustering.fasta")


# write.table(ccrB_tbl,"protein_blast/ccr_combined/ccrB_tbl_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# write.table(ccrB_uniq_tbl,"protein_blast/ccr_combined/ccrB_uniq_tbl_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# to_print = paste0(ccrB_uniq_tbl$uniq_fasta_ID,' ',ccrB_uniq_tbl$subtype_color,' ',ccrB_uniq_tbl$ccr_subtype)
# writeLines(to_print,con = "protein_blast/ccr_combined/ccrB_subtype_colors.txt")
# 
# ccrB_dist_phylo = as.phylo(fit)
# write.tree(ccrB_dist_phylo,file = "protein_blast/ccr_combined/ccrB_dist_tree_QC_18.nwk")


#### ccrC ####

tbl = all_tbl[which(all_tbl$ccr_type=="ccrC"),]
dist_mat = as.matrix(read.table("protein_blast/ccr_combined/ccrC_muscle_dist.txt",sep = "\t",header=T, row.names=1))
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
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrC_tbl = tbl[grep('ccrC',tbl$fasta_ID),]
table(as.vector(IWG_ccrC_tbl$fasta_ID),as.vector(IWG_ccrC_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrCn',as.vector(tbl$ccr_subtype))
ccrC_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrC_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrC_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrC_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrC_tbl$ccr_subtype[idx] = t
  }
  ccrC_tbl$subtype_color[idx] = col
}

ccrC_uniq_tbl = ccrC_tbl[which(!duplicated(as.vector(ccrC_tbl$uniq_fasta_ID))),]
table(ccrC_tbl$ccr_subtype)

ccrC_sub_df = as.data.frame.matrix(table(as.vector(ccrC_tbl$uniq_fasta_ID),ccrC_tbl$ccr_subtype))
fasta_vec = c()
for (i in 1:ncol(ccrC_sub_df)) {
  top_ID = rownames(ccrC_sub_df)[order(ccrC_sub_df[,i], decreasing = TRUE)[1]]
  top_seq = ccrC_uniq_tbl$seq[which(ccrC_uniq_tbl$uniq_fasta_ID==top_ID)]
  printline = paste0('>',colnames(ccrC_sub_df)[i],'\n',top_seq)
  fasta_vec = c(fasta_vec,printline)
}
writeLines(fasta_vec,"protein_blast/subtype_fastas/ccrC_reps_from_protein_clustering.fasta")




# write.table(ccrC_tbl,"protein_blast/ccr_combined/ccrC_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# write.table(ccrC_uniq_tbl,"protein_blast/ccr_combined/ccrC_uniq_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
# to_print = paste0(ccrC_uniq_tbl$uniq_fasta_ID,' ',ccrC_uniq_tbl$subtype_color,' ',ccrC_uniq_tbl$ccr_subtype)
# writeLines(to_print,con = "protein_blast/ccr_combined/ccrC_subtype_colors.txt")
# 
# ccrC_dist_phylo = as.phylo(fit)
# write.tree(ccrC_dist_phylo,file = "protein_blast/ccr_combined/ccrC_dist_tree_QC_22.nwk")



##### Do followup NT blast #####

# blastn -query protein_blast/subtype_fastas/ccrA_reps_from_protein_clustering.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/cds_loc_SCCmec_QC -out protein_blast/followup_nt_blast/ccrA_p_nt_blast.txt -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# blastn -query protein_blast/subtype_fastas/ccrB_reps_from_protein_clustering.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/cds_loc_SCCmec_QC -out protein_blast/followup_nt_blast/ccrB_p_nt_blast.txt -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# blastn -query protein_blast/subtype_fastas/ccrC_reps_from_protein_clustering.fasta -db /srv/data/MPV/projects/SCCmec/blast_DB/cds_loc_SCCmec_QC -out protein_blast/followup_nt_blast/ccrC_p_nt_blast.txt -num_alignments 10000000 -outfmt "6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# combine: cat ccrA_p_nt_blast.txt ccrB_p_nt_blast.txt ccrC_p_nt_blast.txt > ccr_all_p_nt_blast.txt
## extract ccr cds seqs
# python3 ../../extract_cds_from_cds_blast_combined.py ccr_all_p_nt_blast.txt /srv/data/DB/refseq/cds/Staphylococcus ccr_all_seqs.fasta 80 80
### add IWG
# cat ../../references_ccr/IWG_ccrA_all.fasta ../../references_ccr/IWG_ccrB_all.fasta ../../references_ccr/IWG_ccrC_all.fasta ccr_all_seqs.fasta > ccr_all_with_IWG.fasta
# uniquify:
# python3 /srv/data/MPV/THEJ/scripts/uniquify_generic.py ccr_all_with_IWG.fasta ccr_all_uniq.fasta ccr_all_uniq_table.txt
# muscle align:
# sbatch -c 8 --mem=24G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccr_all_uniq.fasta -out ccr_all_uniq_muscle.fasta"
# dist matrix:
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccr_all_uniq_muscle.fasta > ccr_all_muscle_dist.txt




#### ccr all followup NT ####

tbl = read.table("protein_blast/followup_nt_blast/ccr_all_uniq_table.txt",sep = "\t",row.names=NULL,header = T,comment.char = "",check.names = F,quote = "")
tbl$uniq_fasta_ID = paste0(as.vector(tbl$uniq_ID),'|',as.vector(tbl$seq_count))
tbl$GCF_ID = unlist(lapply(as.vector(tbl$ID), function(x) paste0(strsplit(x,"_")[[1]][4:5],collapse="_")))
tbl$source = "RefSeq"
tbl$source[which(tbl$GCF_ID=="NA_NA")] = "IWG_reference"

tbl$species = unlist(lapply(as.vector(tbl$fasta_ID), function(x) strsplit(x,'__')[[1]][1]))
tbl$species[which(tbl$GCF_ID=="NA_NA")] = "IWG_reference"
tbl$species_simplified = get_top_n_from_variable(tbl$species,13)

exclude_idx = which(tbl$source=="RefSeq" & duplicated(tbl$ID))
tbl = tbl[-exclude_idx,]

all_tbl = tbl


dist_mat = read.table("protein_blast/followup_nt_blast/ccr_all_muscle_dist.txt",sep="\t",header=T,row.names=1)
dist_mat = as.matrix(dist_mat)
dist_obj = as.dist(dist_mat)

fit = hclust(dist_obj)

plot(fit)


ccr_groups = cutree(fit,h=45)
rect.hclust(fit,h=45)
sort(unique(ccr_groups))

tbl$ccr_group = NA
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_group[which(tbl$uniq_fasta_ID==name)] = group
}

ccr_dist_phylo = as.phylo(fit)
write.tree(ccr_dist_phylo,file = "protein_blast/followup_nt_blast/ccr_all_dist_tree_QC.nwk")


tbl$ccr_group[grep('ccrA',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrB',tbl$fasta_ID)]
tbl$ccr_group[grep('ccrC',tbl$fasta_ID)]

tbl$ccr_type = paste0('ccrN',as.vector(tbl$ccr_group))
tbl$ccr_type[which(tbl$ccr_group==2)] = 'ccrA'
tbl$ccr_type[which(tbl$ccr_group==1)] = 'ccrB'
tbl$ccr_type[which(tbl$ccr_group==3)] = 'ccrC'
table(tbl$ccr_type)
tbl[grep('IWG',tbl$fasta_ID),]
tbl$ccr_type = factor(tbl$ccr_type)



all_tbl = tbl


all_tbl_uniq = all_tbl[which(!duplicated(as.vector(tbl$uniq_fasta_ID))),]

freq_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$source))
all_tbl_uniq$IWG_references = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="IWG_reference")]))
all_tbl_uniq$IWG_reference = 1
all_tbl_uniq$IWG_reference[which(all_tbl_uniq$IWG_references == 0)] = 0

all_tbl_uniq$RefSeq_count = unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) freq_tbl$Freq[which(freq_tbl$Var1 == x & freq_tbl$Var2=="RefSeq")]))

sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))
sp_tbl = as.data.frame(table(all_tbl$uniq_ID,all_tbl$species_simplified))
sp_tbl = sp_tbl[which(sp_tbl$Freq>0),]
all_tbl_uniq$species_simplified =unlist(lapply(as.vector(all_tbl_uniq$uniq_ID), function(x) paste0(as.vector(sp_tbl$Var2)[which(sp_tbl$Var1==x)],collapse = ',')))

table(all_tbl$source,all_tbl$ccr_type)

write.table(all_tbl,"protein_blast/followup_nt_blast/ccr_table_QC.txt",sep = "\t",quote = FALSE,row.names=FALSE)

write.table(all_tbl_uniq,"protein_blast/followup_nt_blast/ccr_table_uniq_QC.txt",sep = "\t",quote = FALSE,row.names=FALSE)



ccrA_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrA"),]
ccrB_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrB"),]
ccrC_uniq_tbl = all_tbl_uniq[which(all_tbl_uniq$ccr_type=="ccrC"),]
ccrA_uniq_fasta_print = paste0('>',as.vector(ccrA_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrA_uniq_tbl$seq))
writeLines(ccrA_uniq_fasta_print,"protein_blast/followup_nt_blast/ccrA_uniq_seqs.fasta")
ccrB_uniq_fasta_print = paste0('>',as.vector(ccrB_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrB_uniq_tbl$seq))
writeLines(ccrB_uniq_fasta_print,"protein_blast/followup_nt_blast/ccrB_uniq_seqs.fasta")
ccrC_uniq_fasta_print = paste0('>',as.vector(ccrC_uniq_tbl$uniq_fasta_ID),'\n',as.vector(ccrC_uniq_tbl$seq))
writeLines(ccrC_uniq_fasta_print,"protein_blast/followup_nt_blast/ccrC_uniq_seqs.fasta")

### align ccrA, ccrB, ccrC individually using muscle
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrA_uniq_seqs.fasta -out ccrA_uniq_muscle.fasta"
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrB_uniq_seqs.fasta -out ccrB_uniq_muscle.fasta"
# sbatch -c 4 --mem=12G --time=24:00:00 -J "muscle" -p daytime --wrap="muscle -in ccrC_uniq_seqs.fasta -out ccrC_uniq_muscle.fasta"
#
### do pairwise sim matrices
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrA_uniq_muscle.fasta > ccrA_muscle_dist.txt
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrB_uniq_muscle.fasta > ccrB_muscle_dist.txt
# python3 /srv/data/MPV/THEJ/scripts/pairwise_distance_percent_ignore_gaps.py ccrC_uniq_muscle.fasta > ccrC_muscle_dist.txt


all_tbl_uniq$ccr_subtype = "-"
all_tbl$ccr_subtype = "-"





#### ccrA ####

tbl = all_tbl[which(all_tbl$ccr_type=="ccrA"),]
dist_mat = as.matrix(read.table("protein_blast/followup_nt_blast/ccrA_muscle_dist.txt",sep = "\t",header=T, row.names=1))
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
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrA_tbl = tbl[grep('ccrA',tbl$fasta_ID),]
table(as.vector(IWG_ccrA_tbl$fasta_ID),as.vector(IWG_ccrA_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrAn',as.vector(tbl$ccr_subtype))
ccrA_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrA_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrA_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrA_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrA_tbl$ccr_subtype[idx] = t
  }
  ccrA_tbl$subtype_color[idx] = col
}

ccrA_uniq_tbl = ccrA_tbl[which(!duplicated(as.vector(ccrA_tbl$uniq_fasta_ID))),]


write.table(ccrA_tbl,"protein_blast/followup_nt_blast/ccrA_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrA_uniq_tbl,"protein_blast/followup_nt_blast/ccrA_uniq_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrA_uniq_tbl$uniq_fasta_ID,' ',ccrA_uniq_tbl$subtype_color,' ',ccrA_uniq_tbl$ccr_subtype)
writeLines(to_print,con = "protein_blast/followup_nt_blast/ccrA_subtype_colors.txt")

ccrA_dist_phylo = as.phylo(fit)
write.tree(ccrA_dist_phylo,file = "protein_blast/followup_nt_blast/ccrA_dist_tree_QC_22.nwk")





#### ccrB ####

tbl = all_tbl[which(all_tbl$ccr_type=="ccrB"),]
dist_mat = as.matrix(read.table("protein_blast/followup_nt_blast/ccrB_muscle_dist.txt",sep = "\t",header=T, row.names=1))
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
ccr_groups = cutree(fit,h=18)
rect.hclust(fit,h=18)
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrB_tbl = tbl[grep('ccrB',tbl$fasta_ID),]
table(as.vector(IWG_ccrB_tbl$fasta_ID),as.vector(IWG_ccrB_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrBn',as.vector(tbl$ccr_subtype))
ccrB_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrB_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrB_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrB_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrB_tbl$ccr_subtype[idx] = t
  }
  ccrB_tbl$subtype_color[idx] = col
}
table(ccrB_tbl$ccr_subtype)
ccrB_tbl$ccr_subtype[which(ccrB_tbl$uniq_fasta_ID %in% names(ccr_groups)[which(ccr_groups==10)])] = "ccrBn3"
table(ccrB_tbl$ccr_subtype)


ccrB_uniq_tbl = ccrB_tbl[which(!duplicated(as.vector(ccrB_tbl$uniq_fasta_ID))),]
table(ccrB_tbl$ccr_subtype)

write.table(ccrB_tbl,"protein_blast/followup_nt_blast/ccrB_tbl_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrB_uniq_tbl,"protein_blast/followup_nt_blast/ccrB_uniq_tbl_QC_18.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrB_uniq_tbl$uniq_fasta_ID,' ',ccrB_uniq_tbl$subtype_color,' ',ccrB_uniq_tbl$ccr_subtype)
writeLines(to_print,con = "protein_blast/followup_nt_blast/ccrB_subtype_colors.txt")

ccrB_dist_phylo = as.phylo(fit)
write.tree(ccrB_dist_phylo,file = "protein_blast/followup_nt_blast/ccrB_dist_tree_QC_18.nwk")



#### ccrC ####

tbl = all_tbl[which(all_tbl$ccr_type=="ccrC"),]
dist_mat = as.matrix(read.table("protein_blast/followup_nt_blast//ccrC_muscle_dist.txt",sep = "\t",header=T, row.names=1))
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
new_ccr_count = 1
unique(ccr_groups)
for (i in 1:length(ccr_groups)) {
  name = names(ccr_groups)[i]
  group = ccr_groups[i]
  tbl$ccr_subtype[which(tbl$uniq_fasta_ID==name)] = group
}

IWG_ccrC_tbl = tbl[grep('ccrC',tbl$fasta_ID),]
table(as.vector(IWG_ccrC_tbl$fasta_ID),as.vector(IWG_ccrC_tbl$ccr_subtype))
sort(unique(ccr_groups))

tbl$ccr_subtype = paste0('ccrCn',as.vector(tbl$ccr_subtype))
ccrC_tbl = tbl
subtype_colors = triplet_color_vec[1:12]
ccr_subtype_vec = unique(tbl$ccr_subtype)
ccrC_tbl$subtype_color = NA

for (i in 1:length(subtype_colors)) {
  col = subtype_colors[i]
  ccr_type = ccr_subtype_vec[i]
  idx = which(ccrC_tbl$ccr_subtype==ccr_type)
  names = as.vector(ccrC_tbl$fasta_ID[idx])
  test = names[grep('IWG',names)]
  if (length(test)>0) {
    print(test)
    t = strsplit(test[1],'_')[[1]][1]
    ccrC_tbl$ccr_subtype[idx] = t
  }
  ccrC_tbl$subtype_color[idx] = col
}

ccrC_uniq_tbl = ccrC_tbl[which(!duplicated(as.vector(ccrC_tbl$uniq_fasta_ID))),]
table(ccrC_tbl$ccr_subtype)




write.table(ccrC_tbl,"protein_blast/followup_nt_blast/ccrC_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
write.table(ccrC_uniq_tbl,"protein_blast/followup_nt_blast/ccrC_uniq_tbl_QC_22.txt",sep = "\t",quote = FALSE,row.names=FALSE)
to_print = paste0(ccrC_uniq_tbl$uniq_fasta_ID,' ',ccrC_uniq_tbl$subtype_color,' ',ccrC_uniq_tbl$ccr_subtype)
writeLines(to_print,con = "protein_blast/followup_nt_blast/ccrC_subtype_colors.txt")

ccrC_dist_phylo = as.phylo(fit)
write.tree(ccrC_dist_phylo,file = "protein_blast/followup_nt_blast/ccrC_dist_tree_QC_22.nwk")

ccr_combined_tbl = rbind(ccrA_tbl,ccrB_tbl,ccrC_tbl)
write.table(ccr_combined_tbl,"protein_blast/followup_nt_blast/ccr_combined_table.txt",sep = "\t",quote = FALSE,row.names=FALSE)
ccr_combined_uniq_tbl = rbind(ccrA_uniq_tbl,ccrB_uniq_tbl,ccrC_uniq_tbl)
write.table(ccr_combined_uniq_tbl,"protein_blast/followup_nt_blast/ccr_combined_uniq_table.txt",sep = "\t",quote = FALSE,row.names=FALSE)

table(ccr_combined_tbl$ccr_subtype)
table(ccr_combined_uniq_tbl$ccr_subtype)


#### combine ccr types and add mec complex ####

ccr_all_tbl = ccr_combined_tbl
GCF_IDs = unique(ccr_all_tbl$GCF_ID)
GCF_IDs = GCF_IDs[which(!GCF_IDs=="NA_NA")]

A_vec = c()
A_count = c()
B_vec = c()
B_count = c()
C_vec = c()
C_count = c()
for (ID in GCF_IDs) {
  sub_A = ccrA_tbl[which(ccrA_tbl$GCF_ID==ID),]
  sub_B = ccrB_tbl[which(ccrB_tbl$GCF_ID==ID),]
  sub_C = ccrC_tbl[which(ccrC_tbl$GCF_ID==ID),]
  A_count = c(A_count,nrow(sub_A))
  B_count = c(B_count,nrow(sub_B))
  C_count = c(C_count,nrow(sub_C))
  A_vec = c(A_vec,paste0(as.vector(sub_A$ccr_subtype),collapse=','))
  B_vec = c(B_vec,paste0(as.vector(sub_B$ccr_subtype),collapse=','))
  C_vec = c(C_vec,paste0(as.vector(sub_C$ccr_subtype),collapse=','))
}

ccr_simple_tbl = data.frame("GCF_ID"=GCF_IDs,"ccrA_count"=A_count,"ccrB_count"=B_count,"ccrC_count"=C_count,'ccrA_types'=A_vec,'ccrB_types'=B_vec,'ccrC_types'=C_vec)
single_AB_tbl = ccr_simple_tbl[which(ccr_simple_tbl$ccrA_count==1 & ccr_simple_tbl$ccrB_count==1),]
table(single_AB_tbl$ccrA_types,single_AB_tbl$ccrB_types)

test_IDs = ccr_simple_tbl$GCF_ID[which(ccr_simple_tbl$ccrA_types=="ccrA1" & ccr_simple_tbl$ccrB_types=="ccrB1")]

#### mec ####

setwd("/Volumes/data/MPV/projects/SCCmec/mecA_complex/")
mec_table = read.table("blast/mec_all_table.txt",sep = "\t",row.names=NULL,header = T)
mec_table$GCF_ID = unlist(lapply(as.vector(mec_table$ID), function(x) paste0(strsplit(x,'_')[[1]][1:2],collapse = "_")))

mec_uniq_table = read.table("blast/mec_all_refs_uniq_table.txt",sep = "\t",row.names=1,header = F)
colnames(mec_uniq_table) = c("seq_ID","seq","fasta_ID","count")
mec_uniq_table$species = unlist(lapply(as.vector(mec_uniq_table$fasta_ID), function(x) strsplit(x,'__')[[1]][2]))
mec_uniq_table$ID = unlist(lapply(as.vector(mec_uniq_table$fasta_ID), function(x) strsplit(x,'__')[[1]][3]))
mec_uniq_table$contig = unlist(lapply(as.vector(mec_uniq_table$fasta_ID), function(x) strsplit(x,'__')[[1]][4]))
mec_uniq_table$position = unlist(lapply(as.vector(mec_uniq_table$fasta_ID), function(x) strsplit(x,'__')[[1]][5]))
mec_uniq_table$GCF_ID = unlist(lapply(as.vector(mec_uniq_table$ID), function(x) paste0(strsplit(x,'_')[[1]][1:2],collapse = "_")))

mec_uniq_aln = read.dna("blast/mec_all_refs_uniq_muscle_cleaned.fasta","fasta")

mec_ref_nwk = read.tree("references_mec/mec_reps/mec_all_mafft_fasttree.nwk")
mec_ref_aln = read.dna("references_mec/mec_reps/mec_all_mafft.fasta","fasta")

mec_uniq_table[which(mec_uniq_table$fasta_ID %in% mec_ref_nwk$tip.label),]
mec_uniq_table[which(is.na(mec_uniq_table$contig)),]



dupl_IDs = as.vector(mec_table$ID[which(duplicated(mec_table$ID))])

mec_table[which(mec_table$ID %in% dupl_IDs),]


mec_dist = dist.dna(mec_uniq_aln,"raw")

fit = hclust(mec_dist,method="complete")
plot(fit)
rect.hclust(fit,h=0.05)
mec_groups = cutree(fit,h=0.05)

mec_table$cluster = NA
mec_table$cluster_gene = NA
mec_uniq_table$cluster = NA
mec_uniq_table$cluster_gene = NA
mec_IWG_table = data.frame('ID' = labels(mec_ref_aln)[grep('_IWG',labels(mec_ref_aln))], 'cluster'=NA, 'cluster_gene'=NA)
# for (i in unique(mec_groups)) {
#   print(paste0(i))
#   seqs = names(mec_groups)[which(mec_groups==i)]
#   #IWG_IDs = seqs[grep('IWG',seqs)]
#   #mec_IWG_table$cluster[which(mec_IWG_table$ID %in% IWG_IDs)] = i
#   seqs = unlist(lapply(seqs, function(x) strsplit(x,'|',fixed=T)[[1]][1]))
#   IDs = as.vector(mec_uniq_table$fasta_ID[which(mec_uniq_table$seq_ID %in% seqs)])
#   mec_types = as.vector(mec_table$ccr_gene[which(mec_table$fasta_ID %in% IDs)])
#   mec_table$cluster[which(mec_table$fasta_ID %in% IDs)] = i
#   mec_uniq_table$cluster[which(mec_uniq_table$seq_ID %in% seqs)] = i
#   pidents = as.vector(mec_table$pident[which(mec_table$fasta_ID %in% IDs)])
#   print(table(mec_types))
#   print(mean(pidents))
# }

sort(unique(mec_groups[fit$order]))
mec_cluster_names = c("mecB","mecD","mecC3","mecC","mecA4","mecA3","mecA1","mecA2","mecA")
mec_cluster_names = mec_cluster_names[order(unique(mec_groups[fit$order]))]
names(mec_cluster_names) = 1:length(mec_cluster_names)

for (i in 1:length(mec_cluster_names)) {
  name = mec_cluster_names[i]
  mec_table$cluster_gene[which(mec_table$cluster==i)] = name
  mec_uniq_table$cluster_gene[which(mec_uniq_table$cluster==i)] = name
  #mec_IWG_table$cluster_gene[which(mec_IWG_table$cluster==i)] = name
}

sort(unique(mec_groups[fit$order]))

mec_count_table = table(mec_uniq_table$ID)
mec_dupl_IDs = names(which(mec_count_table>1))


#single_mec_table = mec_uniq_table[which(!mec_uniq_table$ID %in% mec_dupl_IDs & !is.na(mec_uniq_table$contig)),]

single_mec_table = mec_table[which(!duplicated(mec_table$GCF_ID)),]


mec_single_AB_table = single_mec_table[which(as.vector(single_mec_table$GCF_ID) %in% as.vector(single_AB_tbl$GCF_ID)),]
single_AB_mec_table = single_AB_tbl[which(as.vector(single_AB_tbl$GCF_ID) %in% as.vector(mec_single_AB_table$GCF_ID)),]

mec_single_AB_table = mec_single_AB_table[order(mec_single_AB_table$GCF_ID),]
single_AB_mec_table = single_AB_mec_table[order(single_AB_mec_table$GCF_ID),]

single_AB_mec_table$mec_gene = mec_single_AB_table$mec_gene

single_AB_mec_table$species = mec_table$species[match(single_AB_mec_table$GCF_ID,mec_table$GCF_ID)]

ccr_mec_freq_table = melt(table(single_AB_mec_table$mec_gene,single_AB_mec_table$ccrA_types,single_AB_mec_table$ccrB_types,single_AB_mec_table$ccrC_types,single_AB_mec_table$species))

colnames(ccr_mec_freq_table) = c("mec_gene","ccrA_gene","ccrB_gene","ccrC_gene","species","freq")
ccr_mec_freq_table=ccr_mec_freq_table[which(ccr_mec_freq_table$freq>0),]
ccr_mec_freq_table$logfreq = log10(ccr_mec_freq_table$freq+1)
ccr_mec_freq_table$species_color = unlist(lapply(ccr_mec_freq_table$species, function(x) get_col(x,sp_color_vec)))

mec_levels = c("mecA","mecA1","mecA2","mecA3","mecA4","mecB","mecC","mecC3","mecD")
ccrA_levels = paste0('ccrA',1:14)
ccrB_levels = paste0('ccrB',1:14)
ccrC_levels = c(0,paste0('ccrC',1:6))

ccr_mec_freq_table$ccrA_gene = factor(as.vector(ccr_mec_freq_table$ccrA_gene),levels=ccrA_levels)
ccr_mec_freq_table$ccrB_gene = factor(as.vector(ccr_mec_freq_table$ccrB_gene),levels=ccrB_levels)
ccr_mec_freq_table$ccrC_gene = factor(as.vector(ccr_mec_freq_table$ccrC_gene),levels=ccrC_levels)



alluvial(ccr_mec_freq_table[,1:4], freq=ccr_mec_freq_table$logfreq,
         #col = ccr_mec_freq_table$species_color,
         #border = ccr_mec_freq_table$species_color,
         alpha = 0.6
         #hide = tit$Freq == 0,
         #cex = 0.7
)


