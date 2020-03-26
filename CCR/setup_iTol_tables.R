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

d_tbl = read.table("ccr_table_uniq_QC_22.txt", sep = "\t", header=T,comment.char = "",quote = "")


itol_lines = setup_color_lines(d_tbl$uniq_fasta_ID,d_tbl$ccr_type,RColorBrewer::brewer.pal(11,"Paired"),legend_title = "ccr_type")
print_template(template,itol_lines,"ccr_type_all_colorstrip.txt",legend_title = "ccr_type")


#### ccrA ####

ccrA_tbl = read.table("ccrA_allotype_uniq_table_QC_22.txt", sep = "\t", header=T,comment.char = "",quote = "")

plot_factor = ccrA_tbl$ccr_allotype
levels(plot_factor)
length(levels(plot_factor))

itol_lines = setup_color_lines(ccrA_tbl$uniq_fasta_ID,plot_factor,c(RColorBrewer::brewer.pal(12,"Paired"),"#A9A9A9"),legend_title = "ccrA allotype")
print_template(template,itol_lines,"colorstrip_ccrA_allotype.txt",legend_title = "ccrA allotype")
