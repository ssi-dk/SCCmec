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



template = readLines("https://raw.githubusercontent.com/ssi-dk/SCCmec/master/itol_template.txt")

#### all ###

uniq_tbl = read.table("ccr_table_uniq_QC_22.txt", sep = "\t", header=T,comment.char = "",quote = "")
tbl = read.table("ccr_table_QC_22.txt", sep = "\t", header=T,comment.char = "",quote = "")

itol_lines = setup_color_lines(uniq_tbl$uniq_fasta_ID,uniq_tbl$ccr_type,RColorBrewer::brewer.pal(11,"Paired"),legend_title = "ccr_type")
print_template(template,itol_lines,"ccr_type_all_colorstrip.txt",legend_title = "ccr_type")



#### investigate ccrB outliers ####

outlier_IDs = unique(as.vector(tbl$fasta_ID[which(tbl$ccr_allotype %in% c("ccrBn1","ccrBn2"))]))
outlier_ccrA_tbl = tbl[which(tbl$fasta_ID %in% outlier_IDs & tbl$ccr_type == "ccrA"),]

table(as.vector(outlier_ccrA_tbl$fasta_ID),as.vector(outlier_ccrA_tbl$ccr_allotype))

outlier_ccrB_tbl = tbl[which(tbl$fasta_ID %in% outlier_IDs & tbl$ccr_type == "ccrB"),]
table(as.vector(outlier_ccrB_tbl$fasta_ID),as.vector(outlier_ccrB_tbl$ccr_allotype))

outlier_ccrC_tbl = tbl[which(tbl$fasta_ID %in% outlier_IDs & tbl$ccr_type == "ccrC"),]
table(as.vector(outlier_ccrC_tbl$fasta_ID),as.vector(outlier_ccrC_tbl$ccr_allotype))

outlier_tbl = tbl[which(tbl$fasta_ID %in% outlier_IDs),]
table(as.vector(outlier_tbl$fasta_ID),as.vector(outlier_tbl$ccr_allotype))


#### ccrA ####

ccrA_tbl = d_tbl[which(d_tbl$ccr_type=="ccrA"),]

plot_factor = factor(as.vector(ccrA_tbl$ccr_allotype))
levels(plot_factor)
length(levels(plot_factor))

itol_lines = setup_color_lines(ccrA_tbl$uniq_fasta_ID_split,plot_factor,c(RColorBrewer::brewer.pal(12,"Paired"),"#A9A9A9"),legend_title = "ccrA allotype")
print_template(template,itol_lines,"colorstrip_ccrA_allotype.txt",legend_title = "ccrA allotype")

itol_lines = setup_color_lines(ccrA_tbl$uniq_fasta_ID_split,ccrA_tbl$source,c("#000000","#FFFFFF"),legend_title = "IWG reference")
print_template(template,itol_lines,"colorstrip_ccrA_IWGref.txt",legend_title = "IWG_reference")


#### ccrB ####

ccrB_tbl = d_tbl[which(d_tbl$ccr_type=="ccrB"),]

plot_factor = factor(as.vector(ccrB_tbl$ccr_allotype))
levels(plot_factor)
length(levels(plot_factor))

itol_lines = setup_color_lines(ccrB_tbl$uniq_fasta_ID_split,plot_factor,c(RColorBrewer::brewer.pal(11,"Paired")),legend_title = "ccrB allotype")
print_template(template,itol_lines,"colorstrip_ccrB_allotype.txt",legend_title = "ccrB allotype")

itol_lines = setup_color_lines(ccrB_tbl$uniq_fasta_ID_split,ccrB_tbl$source,c("#000000","#FFFFFF"),legend_title = "IWG reference")
print_template(template,itol_lines,"colorstrip_ccrB_IWGref.txt",legend_title = "IWG_reference")
