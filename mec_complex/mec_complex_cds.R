

setwd("/Users/thej/Documents/GitHub/SCCmec/mec_complex/cds_blast_loc/")

blast_header = c("query","subject","pident","len","qlen","slen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sstrand")

#### functions ####
get_positions = function(pos_string) {
  gene_complete = "complete"
  if (substr(pos_string,1,10) == "complement") {
    pos = substr(pos_string,12,(nchar(pos_string)-1))
    if (substr(pos,1,4) == "join") {
      pos = substr(pos,6,(nchar(pos_string)-1))
      split_pos = strsplit(pos,',',fixed=TRUE)[[1]]
      spos1 = strsplit(split_pos[1],'..',fixed=TRUE)[[1]][1]
      epos1 = strsplit(split_pos[1],'..',fixed=TRUE)[[1]][2]
      spos2 = strsplit(split_pos[2],'..',fixed=TRUE)[[1]][1]
      epos2 = strsplit(split_pos[2],'..',fixed=TRUE)[[1]][2]
      gene_complete = "fragmented"
    } else {
      spos1 = strsplit(pos,'..',fixed=TRUE)[[1]][1]
      epos1 = strsplit(pos,'..',fixed=TRUE)[[1]][2]
      spos2 = NA
      epos2 = NA
    }
    
  } else if (substr(pos_string,1,4) == "join") {
    pos = substr(pos_string,6,(nchar(pos_string)-1))
    split_pos = strsplit(pos,',',fixed=TRUE)[[1]]
    spos1 = strsplit(split_pos[1],'..',fixed=TRUE)[[1]][1]
    epos1 = strsplit(split_pos[1],'..',fixed=TRUE)[[1]][2]
    spos2 = strsplit(split_pos[2],'..',fixed=TRUE)[[1]][1]
    epos2 = strsplit(split_pos[2],'..',fixed=TRUE)[[1]][2]
    gene_complete = "fragmented"
  }
  else {
    pos = pos_string
    spos1 = strsplit(pos,'..',fixed=TRUE)[[1]][1]
    epos1 = strsplit(pos,'..',fixed=TRUE)[[1]][2]
    spos2 = NA
    epos2 = NA
  }
  return(c(pos,spos1,epos1,spos2,epos2))
}

split_positions

get_start_pos = function(pos_string) {
  if (substr(pos_string,1,1)=='<') {
    pos = substr(pos_string,2,nchar(pos_string))
  } else {
    pos = pos_string
  }
  return(pos)
}

get_end_pos = function(pos_string) {
  if (substr(pos_string,1,1)=='>') {
    pos = substr(pos_string,2,nchar(pos_string))
  } else {
    pos = pos_string
  }
  return(pos)
}

split_join_pos = function(pos_string) {
  pos_split = strsplit(pos_string,',')[[1]]
  if (length(pos_split)>1) {
    return_obj = pos_string
  } else {
    return_obj = c(pos_split[1],pos_split[0])
  }
  return(pos_split)
}
tbl = mecA_tbl

add_data_columns <- function(tbl) {
  tbl$species = unlist(lapply(as.vector(tbl$subject),function(x) strsplit(x,"__GCF")[[1]][1]))
  tbl$assembly_ID = unlist(lapply(as.vector(tbl$subject),function(x) strsplit(x,"__")[[1]][2]))
  tbl$GCF_ID = unlist(lapply(as.vector(tbl$assembly_ID),function(x) strsplit(x,"_")[[1]][2]))
  contig_temp = unlist(lapply(as.vector(tbl$subject),function(x) strsplit(x,"|",fixed = TRUE)[[1]][2]))
  tbl$contig = unlist(lapply(as.vector(contig_temp),function(x) strsplit(x,"_cds")[[1]][1]))
  loc_temp = unlist(lapply(as.vector(tbl$subject),function(x) strsplit(x,"loc=")[[1]][2]))
  tbl$loc_string = loc_temp
  tbl$assembly_strand = "plus"
  tbl$assembly_strand[grep('complement',loc_temp)] = "minus"
  pos_df = t(data.frame(lapply(loc_temp,get_positions)))
  colnames(pos_df) = c("gene_pos","gene_start_pos","gene_end_pos","gene_start_pos_2","gene_end_pos_2")
  rownames(pos_df) = rownames(tbl)
  tbl = data.frame(tbl,pos_df)
  tbl$gene_start_num = as.numeric(unlist(lapply(as.vector(tbl$gene_start_pos),get_start_pos)))
  tbl$gene_end_num = as.numeric(unlist(lapply(as.vector(tbl$gene_end_pos),get_end_pos)))
  tbl$gene_complete = "complete"
  partial_index = c(grep('>',tbl$gene_pos),grep('<',tbl$gene_pos))
  tbl$gene_complete[partial_index] = "partial"
  tbl$gene_split = "complete"
  tbl$gene_split[grep('join',loc_temp)] = "fragmented"
  #tbl$gene_length = tbl$gene_end_num-tbl$gene_start_num+1
  return(tbl)
}


#### load data ####
mecA_tbl = read.table("mecA_cds_blast.txt",sep = "\t", header = F)
colnames(mecA_tbl) = blast_header
mecA_tbl = add_data_columns(mecA_tbl)

mecI_tbl = read.table("mecI_cds_blast.txt",sep = "\t", header = F)
colnames(mecI_tbl) = blast_header
mecI_tbl = add_data_columns(mecI_tbl)

mecR2_tbl = read.table("mecRtype2_cds_blast.txt",sep = "\t", header = F)
colnames(mecR2_tbl) = blast_header
mecR2_tbl = add_data_columns(mecR2_tbl)

mecR4c_tbl = read.table("mecRtype4c_cds_blast.txt",sep = "\t", header = F)
colnames(mecR4c_tbl) = blast_header
mecR4c_tbl = add_data_columns(mecR4c_tbl)

mecR4e_tbl = read.table("mecRtype4e_cds_blast.txt",sep = "\t", header = F)
colnames(mecR4e_tbl) = blast_header
mecR4e_tbl = add_data_columns(mecR4e_tbl)

mecR5_tbl = read.table("mecRtype5_cds_blast.txt",sep = "\t", header = F)
colnames(mecR5_tbl) = blast_header
mecR5_tbl = add_data_columns(mecR5_tbl)

IS1272_tbl = read.table("IS1272_cds_blast.txt",sep = "\t", header = F)
colnames(IS1272_tbl) = blast_header
IS1272_tbl = add_data_columns(IS1272_tbl)

IS431_tbl = read.table("IS431_cds_blast.txt",sep = "\t", header = F)
colnames(IS431_tbl) = blast_header
IS431_tbl = add_data_columns(IS431_tbl)

mecA_assembly_tbl = read.table("../assembly_blast/mecA_assembly_blast.txt",sep = "\t", header = F)
colnames(mecA_assembly_tbl) = blast_header
mecA_assembly_tbl$species = unlist(lapply(as.vector(mecA_assembly_tbl$subject), function(x) strsplit(x,'__')[[1]][1]))
mecA_assembly_tbl$assembly_ID = unlist(lapply(as.vector(mecA_assembly_tbl$subject), function(x) strsplit(x,'__')[[1]][2]))
mecA_assembly_tbl$GCF_ID = unlist(lapply(as.vector(mecA_assembly_tbl$assembly_ID), function(x) strsplit(x,'_')[[1]][2]))
mecA_assembly_tbl$contig = unlist(lapply(as.vector(mecA_assembly_tbl$subject), function(x) strsplit(x,'__')[[1]][3]))
mecA_assembly_sub = mecA_assembly_tbl[which(mecA_assembly_tbl$pident>=98 & mecA_assembly_tbl$len>=1800),]
table(table(mecA_assembly_sub$subject))

IS431_sub = IS431_tbl[which(IS431_tbl$pident>=95),]
plot(sort(IS431_sub$slen))


mecA_sub = mecA_tbl[which(mecA_tbl$pident>=98 & mecA_tbl$gene_complete=="complete" & mecA_tbl$slen>=1800),]
table(mecA_sub$assembly_ID)
table(table(mecA_sub$assembly_ID))
# 8955 isolates with 1 mecA gene, 33 isolates with 2, 1 isolate with 4

n = 1
n = 1000
ass_to_cds_idx = c()
for (n in 1:nrow(mecA_sub)) {
  idx = which(mecA_assembly_sub$GCF_ID==mecA_sub$GCF_ID[n] & mecA_assembly_sub$contig == mecA_sub$contig[n] & 
                mecA_assembly_sub$sstart >= (mecA_sub$gene_start_num[n]-3000) & mecA_assembly_sub$sstart <= (mecA_sub$gene_start_num[n]+3000))
  
  if (length(idx)==1) {
    ass_to_cds_idx = c(ass_to_cds_idx,idx)
  } else if (length(idx>1)) {
    ass_to_cds_idx = c(ass_to_cds_idx,idx[1])
    print(idx)
  } else {
    ass_to_cds_idx = c(ass_to_cds_idx,NA)
  }
}
length(ass_to_cds_idx)
dim(mecA_sub)
mecA_sub$contig_length = mecA_assembly_sub$slen[ass_to_cds_idx]


mat = matrix(nrow=0,ncol=24)
for (n in 1:nrow(mecA_sub)) {
  assembly = as.character(as.vector(mecA_sub$assembly_ID)[n])
  contig = as.character(as.vector(mecA_sub$contig)[n])
  mecA_sub[n,]
  mecA_direction = mecA_sub$assembly_strand[n]
  mecA_start = mecA_sub$gene_start_num[n]
  mecA_end = mecA_sub$gene_end_num[n]
  mecA_len = mecA_sub$slen[n]
  if (mecA_direction == 'plus') {
    IS431_start_ID_tbl = IS431_tbl[which(IS431_tbl$assembly_ID==assembly & IS431_tbl$contig==contig & IS431_tbl$gene_end_num>=(mecA_start-5000) & IS431_tbl$gene_end_num <= mecA_start ),]
    IS431_end_ID_tbl = IS431_tbl[which(IS431_tbl$assembly_ID==assembly & IS431_tbl$contig==contig & IS431_tbl$gene_end_num>=mecA_end & IS431_tbl$gene_end_num < (mecA_end+5000) ),]
    IS1272_ID_tbl = IS1272_tbl[which(IS1272_tbl$assembly_ID==assembly & IS1272_tbl$contig==contig & IS1272_tbl$gene_end_num>=(mecA_start-5000) & IS1272_tbl$gene_end_num <= mecA_start ),]
    mecI_ID_tbl = mecI_tbl[which(mecI_tbl$assembly_ID==assembly & mecI_tbl$contig==contig & mecI_tbl$gene_end_num>=(mecA_start-5000) & mecI_tbl$gene_end_num <= mecA_start ),]
    mecR2_ID_tbl = mecR2_tbl[which(mecR2_tbl$assembly_ID==assembly & mecR2_tbl$contig==contig & mecR2_tbl$gene_end_num>=(mecA_start-500) & mecR2_tbl$gene_end_num <= mecA_start ),]
    mecR4c_ID_tbl = mecR4c_tbl[which(mecR4c_tbl$assembly_ID==assembly & mecR4c_tbl$contig==contig & mecR4c_tbl$gene_end_num>=(mecA_start-500) & mecR4c_tbl$gene_end_num <= mecA_start ),]
    #mecR4e_ID_tbl = mecR4e_tbl[which(mecR4e_tbl$assembly_ID==assembly & mecR4e_tbl$contig==contig),]
    #mecR5_ID_tbl = mecR5_tbl[which(mecR5_tbl$assembly_ID==assembly & mecR5_tbl$contig==contig),]
    if (nrow(IS431_start_ID_tbl)>0) {
      if (IS431_start_ID_tbl$assembly_strand=='minus') {
        IS431_start_vec = c(as.vector(IS431_start_ID_tbl$gene_start_pos[1]),as.vector(IS431_start_ID_tbl$slen[1]),as.vector(IS431_start_ID_tbl$gene_end_pos[1]),'<')
      } else {
        IS431_start_vec = c(as.vector(IS431_start_ID_tbl$gene_start_pos[1]),as.vector(IS431_start_ID_tbl$slen[1]),as.vector(IS431_start_ID_tbl$gene_end_pos[1]),'>')
      }
    } else {
      IS431_start_vec = c("-","-","-","-")
    }
    if (nrow(IS1272_ID_tbl)>0) {
      if (IS1272_ID_tbl$assembly_strand=='minus') {
        IS1272_vec = c(as.vector(IS1272_ID_tbl$gene_start_pos[1]),as.vector(IS1272_ID_tbl$slen[1]),as.vector(IS1272_ID_tbl$gene_end_pos[1]),'<')
      } else {
        IS127_vec = c(as.vector(IS1272_ID_tbl$gene_start_pos[1]),as.vector(IS1272_ID_tbl$slen[1]),as.vector(IS1272_ID_tbl$gene_end_pos[1]),'>')
      }
    } else {
      IS1272_vec = c("-","-","-","-")
    }
    if (nrow(mecI_ID_tbl)>0) {
      if (mecI_ID_tbl$assembly_strand=='minus') {
        mecI_vec = c(as.vector(mecI_ID_tbl$gene_start_pos[1]),as.vector(mecI_ID_tbl$slen[1]),as.vector(mecI_ID_tbl$gene_end_pos[1]),'<')
      } else {
        mecI_vec = c(as.vector(mecI_ID_tbl$gene_start_pos[1]),as.vector(mecI_ID_tbl$slen[1]),as.vector(mecI_ID_tbl$gene_end_pos[1]),'>')
      }
    } else {
      mecI_vec = c("-","-","-","-")
    }
    if (nrow(mecR2_ID_tbl)>0) {
      if (mecR2_ID_tbl$assembly_strand=='minus') {
        mecR2_vec = c(as.vector(mecR2_ID_tbl$gene_start_pos[1]),as.vector(mecR2_ID_tbl$slen[1]),as.vector(mecR2_ID_tbl$gene_end_pos[1]),'<')
      } else {
        mecR2_vec = c(as.vector(mecR2_ID_tbl$gene_start_pos[1]),as.vector(mecR2_ID_tbl$slen[1]),as.vector(mecR2_ID_tbl$gene_end_pos[1]),'>')
      }
    } else {
      mecR2_vec = c("-","-","-","-")
    }
    if (nrow(IS431_end_ID_tbl)>0) {
      if (IS431_end_ID_tbl$assembly_strand=='minus') {
        IS431_end_vec = c(as.vector(IS431_end_ID_tbl$gene_start_pos[1]),as.vector(IS431_end_ID_tbl$slen[1]),as.vector(IS431_end_ID_tbl$gene_end_pos[1]),'<')
      } else {
        IS431_end_vec = c(as.vector(IS431_end_ID_tbl$gene_start_pos[1]),as.vector(IS431_end_ID_tbl$slen[1]),as.vector(IS431_end_ID_tbl$gene_end_pos[1]),'>')
      }
    } else {
      IS431_end_vec = c("-","-","-","-")
    }
    
    
  } else {
    IS431_start_ID_tbl = IS431_tbl[which(IS431_tbl$assembly_ID==assembly & IS431_tbl$contig==contig & IS431_tbl$gene_start_num >= mecA_end & IS431_tbl$gene_start_num<(mecA_end+5000) ) ,]
    IS431_end_ID_tbl = IS431_tbl[which(IS431_tbl$assembly_ID==assembly & IS431_tbl$contig==contig & IS431_tbl$gene_end_num<=mecA_start & IS431_tbl$gene_end_num > (mecA_start-5000) ),]
    IS1272_ID_tbl = IS1272_tbl[which(IS1272_tbl$assembly_ID==assembly & IS1272_tbl$contig==contig & IS1272_tbl$gene_start_num>=mecA_end & IS1272_tbl$gene_start_num < (mecA_end+5000) ),]
    mecI_ID_tbl = mecI_tbl[which(mecI_tbl$assembly_ID==assembly & mecI_tbl$contig==contig & mecI_tbl$gene_start_num>=mecA_end & mecI_tbl$gene_start_num < (mecA_end+5000) ),]
    mecR2_ID_tbl = mecR2_tbl[which(mecR2_tbl$assembly_ID==assembly & mecR2_tbl$contig==contig & mecR2_tbl$gene_start_num>=mecA_end & mecR2_tbl$gene_start_num < (mecA_end+500) ),]
    mecR4c_ID_tbl = mecR4c_tbl[which(mecR4c_tbl$assembly_ID==assembly & mecR4c_tbl$contig==contig & mecR4c_tbl$gene_start_num>=mecA_end & mecR4c_tbl$gene_start_num < (mecA_end+500) ),]
    #mecR4e_ID_tbl = mecR4e_tbl[which(mecR4e_tbl$assembly_ID==assembly & mecR4e_tbl$contig==contig),]
    #mecR5_ID_tbl = mecR5_tbl[which(mecR5_tbl$assembly_ID==assembly & mecR5_tbl$contig==contig),]
    if (nrow(IS431_start_ID_tbl)>0) {
      if (IS431_start_ID_tbl$assembly_strand=='minus') {
        IS431_start_vec = c(as.vector(IS431_start_ID_tbl$gene_start_pos[1]),as.vector(IS431_start_ID_tbl$slen[1]),as.vector(IS431_start_ID_tbl$gene_end_pos[1]),'>')
      } else {
        IS431_start_vec = c(as.vector(IS431_start_ID_tbl$gene_start_pos[1]),as.vector(IS431_start_ID_tbl$slen[1]),as.vector(IS431_start_ID_tbl$gene_end_pos[1]),'<')
      }
    } else {
      IS431_start_vec = c("-","-","-","-")
    }
    if (nrow(IS1272_ID_tbl)>0) {
      if (IS1272_ID_tbl$assembly_strand=='minus') {
        IS1272_vec = c(as.vector(IS1272_ID_tbl$gene_start_pos[1]),as.vector(IS1272_ID_tbl$slen[1]),as.vector(IS1272_ID_tbl$gene_end_pos[1]),'>')
      } else {
        IS127_vec = c(as.vector(IS1272_ID_tbl$gene_start_pos[1]),as.vector(IS1272_ID_tbl$slen[1]),as.vector(IS1272_ID_tbl$gene_end_pos[1]),'<')
      }
    } else {
      IS1272_vec = c("-","-","-","-")
    }
    if (nrow(mecI_ID_tbl)>0) {
      if (mecI_ID_tbl$assembly_strand=='minus') {
        mecI_vec = c(as.vector(mecI_ID_tbl$gene_start_pos[1]),as.vector(mecI_ID_tbl$slen[1]),as.vector(mecI_ID_tbl$gene_end_pos[1]),'>')
      } else {
        mecI_vec = c(as.vector(mecI_ID_tbl$gene_start_pos[1]),as.vector(mecI_ID_tbl$slen[1]),as.vector(mecI_ID_tbl$gene_end_pos[1]),'<')
      }
    } else {
      mecI_vec = c("-","-","-","-")
    }
    if (nrow(mecR2_ID_tbl)>0) {
      if (mecR2_ID_tbl$assembly_strand=='minus') {
        mecR2_vec = c(as.vector(mecR2_ID_tbl$gene_start_pos[1]),as.vector(mecR2_ID_tbl$slen[1]),as.vector(mecR2_ID_tbl$gene_end_pos[1]),'>')
      } else {
        mecR2_vec = c(as.vector(mecR2_ID_tbl$gene_start_pos[1]),as.vector(mecR2_ID_tbl$slen[1]),as.vector(mecR2_ID_tbl$gene_end_pos[1]),'<')
      }
    } else {
      mecR2_vec = c("-","-","-","-")
    }
    if (nrow(IS431_end_ID_tbl)>0) {
      if (IS431_end_ID_tbl$assembly_strand=='minus') {
        IS431_end_vec = c(as.vector(IS431_end_ID_tbl$gene_start_pos[1]),as.vector(IS431_end_ID_tbl$slen[1]),as.vector(IS431_end_ID_tbl$gene_end_pos[1]),'>')
      } else {
        IS431_end_vec = c(as.vector(IS431_end_ID_tbl$gene_start_pos[1]),as.vector(IS431_end_ID_tbl$slen[1]),as.vector(IS431_end_ID_tbl$gene_end_pos[1]),'<')
      }
    } else {
      IS431_end_vec = c("-","-","-","-")
    }
  }
  add_vec = c(IS431_start_vec,IS1272_vec,mecI_vec,mecR2_vec,mecA_start,mecA_len,mecA_end,'>',IS431_end_vec)
  mat = rbind(mat,add_vec)
  
  #IS431_ID_tbl = IS431_tbl[which(IS431_tbl$assembly_ID==assembly & IS431_tbl$contig==contig),]
  #IS1272_ID_tbl = IS1272_tbl[which(IS1272_tbl$assembly_ID==assembly & IS1272_tbl$contig==contig),]
  #mecI_ID_tbl = mecI_tbl[which(mecI_tbl$assembly_ID==assembly & mecI_tbl$contig==contig),]
  #mecR2_ID_tbl = mecR2_tbl[which(mecR2_tbl$assembly_ID==assembly & mecR2_tbl$contig==contig),]
  #mecR4c_ID_tbl = mecR4c_tbl[which(mecR4c_tbl$assembly_ID==assembly & mecR4c_tbl$contig==contig),]
  #mecR4e_ID_tbl = mecR4e_tbl[which(mecR4e_tbl$assembly_ID==assembly & mecR4e_tbl$contig==contig),]
  #mecR5_ID_tbl = mecR5_tbl[which(mecR5_tbl$assembly_ID==assembly & mecR5_tbl$contig==contig),]
}

dim(mat)
dim(mecA_sub)
rownames(mat) = rownames(mecA_sub)
complex_df = as.data.frame(mat)
colnames(complex_df) = c("IS431s_start","IS431s_len","IS431s_end","IS431s_dir","IS1272s_start","IS1272s_len","IS1272s_end","IS1272s_dir","mecI_start","mecI_len","mecI_end","mecI_dir",
                         "mecR_start","mecR_len","mecR_end","mecR_dir","mecA_start","mecA_len","mecA_end","mecA_dir","IS431e_start","IS431e_len","IS431e_end","IS431e_dir")
#rownames(complex_df) = rownames(mecA_sub)

complex_df$assembly = mecA_sub$assembly_ID
complex_df$GCF_ID = mecA_sub$GCF_ID
complex_df$contig = mecA_sub$contig
complex_df$species = mecA_sub$species
complex_df$contig_length = mecA_sub$contig_length
complex_df$assembly_strand = mecA_sub$assembly_strand
plus_vec = as.numeric(as.vector(complex_df$mecA_start[which(complex_df$assembly_strand=="plus")]))
minus_vec = as.numeric(as.vector(complex_df$contig_length[which(complex_df$assembly_strand=="minus")]))-as.numeric(as.vector(complex_df$mecA_end[which(complex_df$assembly_strand=="minus")]))
complex_df$mecA_to_end_length = NA
complex_df$mecA_to_end_length[which(complex_df$assembly_strand=="plus")] = plus_vec
complex_df$mecA_to_end_length[which(complex_df$assembly_strand=="minus")] = minus_vec
complex_df$end_length = complex_df$mecA_to_end_length
complex_df$end_length[which(complex_df$end_length>=5000)] = 5000

mecR_start_num = as.numeric(unlist(lapply(as.character(complex_df$mecR_start), get_start_pos)))
mecR_end_num = as.numeric(unlist(lapply(as.character(complex_df$mecR_end), get_end_pos)))
mecR_min = unlist(lapply(1:length(mecR_start_num), function(x) min(mecR_start_num[x],mecR_end_num[x])))
mecR_max = unlist(lapply(1:length(mecR_start_num), function(x) max(mecR_start_num[x],mecR_end_num[x])))
plus_vec = as.numeric(as.vector(mecR_min[which(complex_df$assembly_strand=="plus")]))
minus_vec = as.numeric(as.vector(complex_df$contig_length[which(complex_df$assembly_strand=="minus")]))-as.numeric(as.vector(mecR_max[which(complex_df$assembly_strand=="minus")]))
complex_df$mecR_to_end_length = NA
complex_df$mecR_to_end_length[which(complex_df$assembly_strand=="plus")] = plus_vec
complex_df$mecR_to_end_length[which(complex_df$assembly_strand=="minus")] = minus_vec
complex_df$mecR_end_length = complex_df$mecR_to_end_length
complex_df$mecR_end_length[which(complex_df$mecR_to_end_length>=5000)] = 5000

plot(sort(complex_df$end_length))

mecR_len_num = as.numeric(as.vector(complex_df$mecR_len))
table(mecR_len_num)
complex_df$mecR_len_simple = "-"
complex_df$mecR_len_simple[which(!is.na(mecR_len_num))] = "NT"
complex_df$mecR_len_simple[which(mecR_len_num>=2000)] = "NT_3918"
complex_df$mecR_len_simple[which(mecR_len_num>=1740 & mecR_len_num<=1770)] = "type2"
complex_df$mecR_len_simple[which(mecR_len_num>1000 & mecR_len_num<1740)] = "NT_1001-1739"
complex_df$mecR_len_simple[which(mecR_len_num>=900 & mecR_len_num<=1000)] = "type4"
complex_df$mecR_len_simple[which(mecR_len_num>110 & mecR_len_num<900)] = "NT_111-984"
complex_df$mecR_len_simple[which(mecR_len_num>=100 & mecR_len_num<=110)] = "type5"
complex_df$mecR_len_simple[which(mecR_len_num<=100)] = "NT_0-99"
complex_df$mecI_present = 1
complex_df$mecI_present[which(complex_df$mecI_len=="-")] = 0
complex_df$IS1272_present = 1
complex_df$IS1272_present[which(complex_df$IS1272s_len=="-")] = 0
complex_df$IS431s_present = 1
complex_df$IS431s_present[which(complex_df$IS431s_len=="-")] = 0
complex_df$mecR_type = complex_df$mecR_len_simple
complex_df$ISs_status = "Absent"
complex_df$ISs_status[which((complex_df$IS431s_present==1 | complex_df$IS1272_present==1) & complex_df$end_length < (as.numeric(complex_df$mecR_len)+1000))]
complex_df$ISs_status[which(complex_df$IS431s_present==1 | complex_df$IS1272_present==1)] = "Present"

complex_df$mec_class = "NT"
complex_df$mec_class[which(complex_df$mecI_present==1 & complex_df$mecR_type=="type2")] = "A"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="type4" & complex_df$IS1272_present==1)] = "B"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="type4" & complex_df$IS1272_present==0)] = "C1/D"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="type4" & complex_df$IS1272_present==0 & complex_df$IS431s_present==1)] = "C1"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="type4" & complex_df$IS1272_present==0 & complex_df$IS431s_present==0 & complex_df$end_length>2000)] = "D"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="type5" & complex_df$IS1272_present==0)] = "C2"
complex_df$mec_class[which(complex_df$mecI_present==0 & complex_df$mecR_type=="-" & complex_df$IS431s_present==1 & complex_df$IS1272_present==0)] = "mecA_only"

table(complex_df$mec_class)

complex_df$mec_class_inferred = "NT"

plot(sort(complex_df$end_length[which(complex_df$mec_class=="A")]))
plot(sort(complex_df$end_length[which(complex_df$mec_class=="B")]))
plot(sort(complex_df$end_length[which(complex_df$mec_class=="C1/D")]))
plot(sort(complex_df$end_length[which(complex_df$mec_class=="C2")]))
plot(sort(complex_df$end_length[which(complex_df$mec_class=="NT")]))


complex_df[which(complex_df$mec_class=="A"),]
complex_df[which(complex_df$mec_class=="B"),]
complex_df[which(complex_df$mec_class=="C1/D"),]
complex_df[which(complex_df$mec_class=="C2"),]
complex_df[which(complex_df$mec_class=="mecA_only"),]
complex_df[which(complex_df$mec_class=="NT"),]


p <- ggplot(data.frame("x"=order(order(complex_df$end_length[which(complex_df$mecR_len=='-')])),"y"=complex_df$end_length[which(complex_df$mecR_len=='-')],"col"=complex_df$ISs_status[which(complex_df$mecR_len=='-')]),
            aes(x=x,y=y,color=col)) + geom_point() + labs(color="Upstream IS431 presence") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("Number of base pairs upstream of mecA in contig in isolates missing mecR")
p

complex_sub = complex_df[which(complex_df$mec_class=="B"),]
p_df = data.frame("x"=order(order(complex_sub$end_length)),"y"=complex_sub$end_length,"col"=complex_sub$ISs_status)
p <- ggplot(p_df,aes(x=x,y=y,color=col)) + geom_point() + labs(color="Upstream IS element presence") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("Number of base pairs upstream of mecA in contig in mec class B isolates")
p

complex_sub = complex_df[which(complex_df$end_length<5000),]
p_df = data.frame("x"=order(order(complex_sub$end_length)),"y"=complex_sub$end_length,"col"=complex_sub$mec_class)
p <- ggplot(p_df,aes(x=x,y=y,color=col)) + geom_point(aes(shape=complex_sub$ISs_status)) + labs(color="mec gene complex") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("Number of base pairs upstream of mecA in contig")
p

p_df = data.frame("x"=complex_sub$mec_class,"y"=complex_sub$end_length,"col"=complex_sub$mec_class)
p <- ggplot(p_df,aes(x=x,y=y,fill=col)) + geom_boxplot() + labs(fill="inferred mec class") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("Number of base pairs upstream of mecA in contig in mec class NT isolates")
p

p_df = data.frame("x"=complex_sub$mecR_type,"y"=complex_sub$end_length,"col"=complex_sub$mecR_type)
p <- ggplot(p_df,aes(x=x,y=y,fill=col)) + geom_boxplot() + labs(fill="mecR type") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("Number of base pairs upstream of mecA in contig in mec class NT isolates")
p

complex_sub = complex_df[which(!complex_df$mecR_len %in% c("-","3918")),]
p_df = data.frame("x"=order(order(as.numeric(as.vector(complex_sub$mecR_len)))),"y"=as.numeric(as.vector(complex_sub$mecR_len)),"col"=complex_sub$mecR_type)
p <- ggplot(p_df,aes(x=x,y=y,color=col)) + geom_point() + labs(color="Upstream IS431 presence") + ylab("Length from start of mecA to end of contig") + xlab("") + 
  ggtitle("mecR length")
p

p <- ggplot(complex_df[which(complex_df$mecR_end_length<5000),],aes(x=ISs_status,y=mecR_end_length)) + geom_boxplot()
p

b = seq(from = 0, to = 5000, by = 100)
p_df = complex_df
p_df$bin_type = .bincode(complex_df$end_length,b)

plot_df = as.data.frame(table(p_df$mec_class,p_df$bin_type))
colnames(plot_df) = c("mec_class","bin","Freq")
plot_df$end_length = as.numeric(as.vector(plot_df$bin))*100

p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="fill", stat="identity")
p
p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")
p

plot_sub = plot_df[which(plot_df$end_length<5000),]
p <- ggplot(plot_sub,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")
p

plot_df = as.data.frame(table(p_df$mecR_type,p_df$bin_type))
colnames(plot_df) = c("mec_class","bin","Freq")
plot_df$end_length = as.numeric(as.vector(plot_df$bin))*100

p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p
p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p

plot_sub = plot_df[which(plot_df$end_length<5000),]
p <- ggplot(plot_sub,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p

b = seq(from = 0, to = 5000, by = 100)
p_df = complex_df
p_df$bin_type = .bincode(complex_df$mecR_end_length,b)

plot_df = as.data.frame(table(p_df$mecR_type,p_df$bin_type))
colnames(plot_df) = c("mec_class","bin","Freq")
plot_df$end_length = as.numeric(as.vector(plot_df$bin))*100

p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p
p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p

plot_sub = plot_df[which(plot_df$end_length<5000),]
p <- ggplot(plot_sub,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p



plot_df = as.data.frame(table(p_df$IS431s_present,p_df$bin_type))
colnames(plot_df) = c("mec_class","bin","Freq")
plot_df$end_length = as.numeric(as.vector(plot_df$bin))*100

p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p
p <- ggplot(plot_df,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p

plot_sub = plot_df[which(plot_df$end_length<5000),]
p <- ggplot(plot_sub,aes(x=end_length,y=Freq,fill=mec_class)) + geom_bar(position="stack", stat="identity")  + scale_fill_manual(values=RColorBrewer::brewer.pal(8,"Set1"))
p

plot(sort(complex_df$mecR_end_length[which(complex_df$IS431s_present==0)]))
plot(sort(complex_df$mecR_end_length[which(complex_df$IS431s_present==1)]))


complex_df[which(complex_df$mec_class=="NT" & complex_df$mecR_type %in% c("type2","type4","type5")),]
complex_df[which(!complex_df$mecR_type %in% c("type2","type4","type5","-")),]

table(complex_df$mecI_present,complex_df$mecR_len_simple)
table(complex_df$IS1272_present,complex_df$mecR_len_simple)
table(complex_df$IS1272_present,complex_df$mecR_len_simple,complex_df$mecI_present)

x_df =complex_df[which(complex_df$mecR_len_simple=="type2" & complex_df$mecI_present==0),]

table(x_df$species)




#### combine with ccr ####

ccr_tbl = read.table("../../CCR/ccr_table_final.txt",sep = "\t",header=TRUE,comment.char = "",check.names = F,quote = "")
head(ccr_tbl)
dim(ccr_tbl)


refseq_tbl = ccr_tbl[which(!ccr_tbl$species=="IWG_reference"),]
dim(refseq_tbl)

isolate_df = as.data.frame(table(refseq_tbl$fasta_ID,refseq_tbl$ccr_type,refseq_tbl$ccr_allotype))
isolate_df = isolate_df[which(isolate_df$Freq>0),]
head(isolate_df)
dim(isolate_df)


GCF_IDs = unique(complex_df$GCF_ID)
complex_df$ccrA_count = 0
complex_df$ccrB_count = 0
complex_df$ccrC_count = 0
complex_df$ccrA_types = '-'
complex_df$ccrB_types = '-'
complex_df$ccrC_types = '-'
complex_df$mecA_count = 0

n = 1
for (n in 1:nrow(complex_df)) {
  mecA_count = length(which(complex_df$GCF_ID==complex_df$GCF_ID[n]))
  GCF_ID = paste0('GCF_',complex_df$GCF_ID[n])
  ccr_sub = ccr_tbl[which(ccr_tbl$GCF_ID==GCF_ID),]
  ccrA_sub = ccr_sub[which(ccr_sub$ccr_type=="ccrA"),]
  ccrB_sub = ccr_sub[which(ccr_sub$ccr_type=="ccrB"),]
  ccrC_sub = ccr_sub[which(ccr_sub$ccr_type=="ccrC"),]
  complex_df$ccrA_count[n] = nrow(ccrA_sub)
  complex_df$ccrB_count[n] = nrow(ccrB_sub)
  complex_df$ccrC_count[n] = nrow(ccrC_sub)
  if (nrow(ccrA_sub)>0) {
    complex_df$ccrA_types[n] = paste0(as.vector(ccrA_sub$ccr_allotype),collapse = ",")
  }
  if (nrow(ccrB_sub)>0) {
    complex_df$ccrB_types[n] = paste0(as.vector(ccrB_sub$ccr_allotype),collapse = ",")
  }
  if (nrow(ccrC_sub)>0) {
    complex_df$ccrC_types[n] = paste0(as.vector(ccrC_sub$ccr_allotype),collapse = ",")
  }
  complex_df$mecA_count[n] = mecA_count
}

head(complex_df)

table(complex_df$mecA_count)

write.table(complex_df,"SCCmec_master_table.txt",sep = "\t",quote = FALSE,row.names=FALSE)

SCCmec_single_mecA = complex_df[which(complex_df$mecA_count==1),]

ccr_single_idx = which((SCCmec_single_mecA$ccrA_count==1 & SCCmec_single_mecA$ccrB_count==1 & SCCmec_single_mecA$ccrC_count %in% c(0,1)) | (SCCmec_single_mecA$ccrA_count==0 & SCCmec_single_mecA$ccrB_count==0 & SCCmec_single_mecA$ccrC_count==1))
SCCmec_confirmed = SCCmec_single_mecA[ccr_single_idx,]
dim(SCCmec_confirmed)

SCCmec_confirmed$SCCmec_type = "NT"
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA1" & SCCmec_confirmed$ccrB_types=="ccrB1" & SCCmec_confirmed$mec_class=="B")] = 1
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA2" & SCCmec_confirmed$ccrB_types=="ccrB2" & SCCmec_confirmed$mec_class=="A")] = 2
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA3" & SCCmec_confirmed$ccrB_types=="ccrB3" & SCCmec_confirmed$mec_class=="A")] = 3
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA2" & SCCmec_confirmed$ccrB_types=="ccrB2" & SCCmec_confirmed$mec_class=="B")] = 4
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="-" & SCCmec_confirmed$ccrB_types=="-" & SCCmec_confirmed$ccrC_types=="ccrC1" & SCCmec_confirmed$mec_class=="C2")] = 5
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA4" & SCCmec_confirmed$ccrB_types=="ccrB4" & SCCmec_confirmed$mec_class=="B")] = 6
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="-" & SCCmec_confirmed$ccrB_types=="-" & SCCmec_confirmed$ccrC_types=="ccrC1" & SCCmec_confirmed$mec_class=="C1")] = 7
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA4" & SCCmec_confirmed$ccrB_types=="ccrB4" & SCCmec_confirmed$mec_class=="A")] = 8
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA1" & SCCmec_confirmed$ccrB_types=="ccrB1" & SCCmec_confirmed$mec_class=="C2")] = 9
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA1" & SCCmec_confirmed$ccrB_types=="ccrB6" & SCCmec_confirmed$mec_class=="C1")] = 10
#SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="ccrA1" & SCCmec_confirmed$ccrB_types=="ccrB6" & SCCmec_confirmed$mec_class=="C1")] = 11
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="-" & SCCmec_confirmed$ccrB_types=="-" & SCCmec_confirmed$ccrC_types=="ccrC2" & SCCmec_confirmed$mec_class=="C2")] = 12
SCCmec_confirmed$SCCmec_type[which(SCCmec_confirmed$ccrA_types=="-" & SCCmec_confirmed$ccrB_types=="-" & SCCmec_confirmed$ccrC_types=="ccrC2" & SCCmec_confirmed$mec_class=="A")] = 13
table(SCCmec_confirmed$SCCmec_type)
table(SCCmec_confirmed$species,SCCmec_confirmed$SCCmec_type)

write.table(SCCmec_confirmed,"SCCmec_master_table_singles.txt",sep = "\t",quote = FALSE,row.names=FALSE)

NT_tbl = SCCmec_confirmed[which(SCCmec_confirmed$SCCmec_type=="NT"),]

NT_df = as.data.frame(table(NT_tbl$ccrA_types,NT_tbl$ccrB_types,NT_tbl$ccrC_types,NT_tbl$mec_class))
NT_df = NT_df[which(NT_df$Freq>0),]
colnames(NT_df) = c("ccrA","ccrB","ccrC","mec_class","Freq")
NT_df[order(NT_df$Freq),]


t_df = as.data.frame(table(SCCmec_confirmed$ccrA_types,SCCmec_confirmed$ccrB_types,SCCmec_confirmed$ccrC_types,SCCmec_confirmed$mec_class,SCCmec_confirmed$SCCmec_type))
t_df = t_df[t_df$Freq>0,] 
t_df[order(t_df$Freq),]


D_tbl = SCCmec_confirmed[which(SCCmec_confirmed$ccrA_types=="ccrA2" & SCCmec_confirmed$ccrB_types=="ccrB2" & SCCmec_confirmed$mec_class=="D"),]

D_df = as.data.frame(table(D_tbl$ccrA_types,D_tbl$ccrB_types,D_tbl$ccrC_types,D_tbl$mec_class,D_tbl$SCCmec_type))


ccrA1_df = SCCmec_single_mecA[which(SCCmec_single_mecA$ccrA_types=="ccrA1"),]
sort(table(ccrA1_df$ccrB_types))

ccrB1_df = SCCmec_single_mecA[which(SCCmec_single_mecA$ccrB_types=="ccrB1"),]
sort(table(ccrA1_df$ccrB_types))
