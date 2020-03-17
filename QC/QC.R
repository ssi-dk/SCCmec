library(ggplot2)
setwd("/Volumes/data/MPV/projects/SCCmec/QC/")
assembly_stats = read.table("/Volumes/data/DB/refseq/Staphylococcus_contig_length_stats_reclassified.txt",sep="\t")

setwd("E:/Github/SCCmec/QC")
assembly_stats = read.table("Staphylococcus_contig_length_stats_reclassified.txt",sep="\t")


sub_stats = assembly_stats[which(assembly_stats$species=="Staphylococcus_epidermidis"),]

plot_dist_table = data.frame("Contig_count"=sub_stats$contigs)
## Histogram
p <- ggplot(plot_dist_table,aes(x=Contig_count)) + geom_histogram(position="dodge",binwidth = 2) + scale_fill_manual(values = RColorBrewer::brewer.pal(3,"Set1"))
p

## Boxplot - number of contigs
p <- ggplot(assembly_stats,aes(x=species,y=contigs,fill=species)) + geom_boxplot() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) +
  ggtitle("Distribution of contig counts, no filtering")
p

## Boxplot - genome size
p <- ggplot(assembly_stats,aes(x=species,y=genome_size,fill=species)) + geom_boxplot() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) +
  ggtitle("Distribution of genome sizes, no filtering")
p

## Filter 
filtered_set = assembly_stats[which(assembly_stats$contigs>300),]  ## 366 genomer fjernes hvis vi filterer => 500 contigs
dim(filtered_set)
# Alternativ
length(filtered_set$contigs) ## lige meget hvad jeg vælger, altså contigs, species, eller size da det vil være det samme

x <- table(assembly_stats$species)  # number of isolates for each speacies
## skriver tabbellen ud
write.table(x, "species_table.txt", quote = FALSE, sep="\t", row.names = FALSE)
y <- table(filtered_set$species)

table(filtered_set$species)/table(assembly_stats$species)*100



#### IQR evaluation of species genome size ####

filter_size_IQR = function(assembly_stats,species,IQR_threshold_factor) {
  sub_df = assembly_stats[which(assembly_stats$species==species),]
  size_vec = sub_df$genome_size
  iqr = IQR(size_vec)
  median = median(size_vec)
  include_idx = which(size_vec>(median-iqr*IQR_threshold_factor) & size_vec<(median+iqr*IQR_threshold_factor))
  include_IDs = rownames(sub_df)[include_idx]
  return(include_IDs)
}

IQR_species_apply = function(assembly_stats, IQR_threshold_factor) {
  species_list = levels(assembly_stats$species)
  include_IDs = c()
  for (species in species_list) {
    sp_IDs = filter_size_IQR(assembly_stats,species,IQR_threshold_factor)
    include_IDs = c(include_IDs,sp_IDs)
  }
  return_df = assembly_stats[which(rownames(assembly_stats) %in% include_IDs),]
  return(return_df)
}



filtered_set = assembly_stats[which(assembly_stats$contigs<=300),]



### updated plot
filtered_stats = IQR_species_apply(filtered_set,IQR_threshold_factor = 2.0)
p <- ggplot(filtered_stats,aes(x=sp,y=genome_size,fill=sp)) + geom_boxplot() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) +
  ggtitle("Genome size distribution after filtering contig counts > 300 and genome sizes outside 2.0 IQR of the median for each species")
p

### whats been removed
table(assembly_stats$species)-table(filtered_stats$species)



contig_filtered_set = assembly_stats[which(assembly_stats$contigs<=300),]

size_filtered_set = IQR_species_apply(assembly_stats,IQR_threshold_factor = 2.5)

both_filtered_set = IQR_species_apply(contig_filtered_set,IQR_threshold_factor = 2.5)


species_removed_summary = data.frame('Species'=levels(assembly_stats$sp),'Original'=as.vector(table(assembly_stats$sp)),
                                     'Failing_contig_requirement'=as.vector(table(assembly_stats$sp))-as.vector(table(contig_filtered_set$sp)),
                                     'Failing_genome_size_requirement'=as.vector(table(assembly_stats$sp))-as.vector(table(size_filtered_set$sp)),
                                     'Failing_any_requirement'=as.vector(table(assembly_stats$sp))-as.vector(table(both_filtered_set$sp)),
                                     'Final_remaining'=as.vector(table(both_filtered_set$sp)))



write.table(species_removed_summary,"excluded_isolates_species_table.txt",sep="\t",quote = FALSE)

original_IDs = as.vector(assembly_stats$sp_ID)
final_IDs = as.vector(both_filtered_set$sp_ID)

excluded_isolates = original_IDs[which(!original_IDs %in% final_IDs)]
excluded_isolates = unlist(lapply(excluded_isolates, function(x) strsplit(x,'_genomic')[[1]][1]))

writeLines(excluded_isolates,"excluded_isolates.txt")

write.table(both_filtered_set,"assembly_stats_post_filtering.txt",sep = "\t")







aureus_start_IDs = assembly_stats$sp_ID[which(assembly_stats$sp=="aureus")]
aureus_filter_IDs = filtered_stats$sp_ID[which(filtered_stats$sp=="aureus")]

filtered_IDs = aureus_start_IDs[which(!aureus_start_IDs %in% aureus_filter_IDs)]

sub_mlst = mlst[which(mlst$V1 %in% filtered_IDs),]

sort(table(sub_mlst$V3))

test = table(sub_mlst$V3)/table(mlst$V3)*100
