## fingerprinting

library(stringr)

process_fp_variant_evidence <- function(fp) {
  # Process samtools mpileup output into a fingerprint dataframe
  
  # FS = uppercase calls, RS = lowercall calls
  fp$ref_FS <- str_count(fp$mpileup, pattern = "[.]")
  fp$ref_RS <- str_count(fp$mpileup, pattern = ",")
  fp$A_FS <- str_count(fp$mpileup, pattern = "A")
  fp$A_RS <- str_count(fp$mpileup, pattern = "a")
  fp$C_FS <- str_count(fp$mpileup, pattern = "C")
  fp$C_RS <- str_count(fp$mpileup, pattern = "c")
  fp$G_FS <- str_count(fp$mpileup, pattern = "G")
  fp$G_RS <- str_count(fp$mpileup, pattern = "g")
  fp$T_FS <- str_count(fp$mpileup, pattern = "T")
  fp$T_RS <- str_count(fp$mpileup, pattern = "t")
  fp$DEL <- str_count(fp$mpileup, pattern = "[*]")
  fp$INS <- str_count(fp$mpileup, pattern = "[+]")
  
  # [variant]_evidence: whether or not there was any support for that variant
  fp$A_evidence <- fp$A_FS>0 & fp$A_RS>0
  fp$C_evidence <- fp$C_FS>0 & fp$C_RS>0
  fp$T_evidence <- fp$T_FS>0 & fp$T_RS>0
  fp$G_evidence <- fp$G_FS>0 & fp$G_RS>0
  fp$INS_evidence <- fp$INS>0
  fp$DEL_evidence <- fp$DEL>0
  
  return(fp)
}

get_fp <- function(bam, ref, ref_fai, sampID, temp_dir = "/data",
                   no_coverage_depth = 2, low_coverage_depth = 30) {
  # Arguments:
  # bam: file path to a sample BAM file
  # ref: file path to a reference genome FASTA file
  # ref_fai: file path to the FASTA index (.fai) for reference
  # sampID: sample ID
  # temp_dir: optional argument for a temporary working directory
  # no_coverage_depth: if site has fewer than this many reads, it has "no coverage"
  # low_coverage_depth: if site has fewer than this many reads, it has "low coverage"
  
  load("fingerprint_snps.RData")
  
  # Run mpileup
  command <- paste0('samtools mpileup -a -l fingerprinting_locations.txt --fasta-ref ', ref,
                 ' ', bam,  '> ', temp_dir,'/fingerprint.txt')
  system(command)
  
  # Read and parse mpileup output into fingerprint
  fp <- read.csv(file = paste0(temp_dir, "/fingerprint.txt"), sep = "\t", 
                 stringsAsFactors = FALSE, header = FALSE, quote = "", fill = FALSE)
  
  colnames(fp) <- c("Chr", "Pos", "Ref", "Depth", "mpileup", "Qual")
  fp$CUR_POS <- paste(fp$Chr, fp$Pos, sep = "-")

  fp <- process_fp_variant_evidence(fp)
    
  fingerprint_snps$CHROM_POS_REF_ALT <- paste(
    fingerprint_snps$Chr, fingerprint_snps$Start, fingerprint_snps$Ref,
    fingerprint_snps$Alt, sep = "-")
  
  # Go through the fingerprint SNPs and compare them to the fingerprint
  nt <- c("A", "C", "T", "G")
  fp_count <- NULL
  
  for (i in 1:nrow(fingerprint_snps)) {
    CHROM_POS_REF_ALT <- fingerprint_snps$CHROM_POS_REF_ALT[i]
    curr.pos <- paste(fingerprint_snps$Chr[i], fingerprint_snps$Start[i], sep = "-")
    curr.ref <- fingerprint_snps$Ref[i]
    curr.alt <- fingerprint_snps$Alt[i]
    tmp <- fp[fp$CUR_POS == curr.pos,]
    
    depth_total <- tmp$Depth
    if (depth_total >= low_coverage_depth) {
      evidence = "Covered"
    } else if (depth_total < no_coverage_depth) {
      evidence = "No coverage"
    } else if (depth_total < low_coverage_depth) {
      evidence <- "Low coverage"
    }
    
    # Processing for SNVs
    if (curr.alt %in% nt & curr.ref %in% nt) {
      if (sum(tmp$A_evidence, tmp$T_evidence, tmp$C_evidence, tmp$G_evidence,
             tmp$INS_evidence, tmp$DEL_evidence) > 1) {
        evidence <- paste(evidence, "Multiallelic Locus", sep = ";")
      }
      
      if (curr.alt == "A" & tmp$A_evidence) {
        depth_alt <- tmp$A_FS + tmp$A_RS
        AF <- depth_alt/depth_total
      } else if (curr.alt == "C" & tmp$C_evidence) {
        depth_alt <- tmp$C_FS+tmp$C_RS
        AF <- depth_alt/depth_total
      } else if (curr.alt == "G" & tmp$G_evidence) {
        depth_alt <- tmp$G_FS+tmp$G_RS
        AF <- depth_alt/depth_total
      } else if (curr.alt == "T" & tmp$T_evidence) {
        depth_alt <- tmp$T_FS+tmp$T_RS
        AF <- depth_alt/depth_total
      } else {
        depth_alt = 0
        AF <- 0
      }
    }
    
    fp_count <- rbind(fp_count, 
                      cbind(CHROM_POS_REF_ALT, depth_total, depth_alt, AF, evidence))
  }
  
  # Clean up the fingerprint count vector
  fp_count <- data.frame(fp_count)
  
  for(i in 2:4){fp_count[,i] <- as.numeric(paste(fp_count[,i]))}

  
  fp_vector <- rep(NA, nrow(fp_count))
  fp_vector[fp_count$AF<0.15] <- "0"
  fp_vector[fp_count$AF>0.85] <- "2"
  fp_vector[fp_count$AF>0.25 & fp_count$AF<0.75] <- "1"
  fp_vector[fp_count$depth_total<2] <- "N"
  fp_vector[is.na(fp_vector)] <- "U"
  
  
  
  fp_vector <- data.frame(fp_vector)
  
  colnames(fp_vector) <- sampID
  colnames(fp_count) <- paste0(sampID, "_", colnames(fp_count))
  
  return(list(fp_count,fp_vector))
  
}

# Parse command line arguments and run the above function to get
# fingerprint vector and count
args <- commandArgs(trailingOnly = TRUE)

bam <- args[1]
ref <- args[2]
ref_fai <- args[3]
samp_id <- args[4]
fingerprint_path <- args[5]
fingerprint_vector_path <- args[6]
fingerprint_count_path <- args[7]

fp <- get_fp(bam, ref, ref_fai, samp_id)

save(fp, file = fingerprint_path)
fp_count <- fp[[1]]
fp_vector <- fp[[2]]

write.table(fp_vector, file = fingerprint_vector_path, quote = FALSE, row.names = FALSE)

write.table(fp_count, file = fingerprint_count_path, quote = FALSE, row.names = FALSE)
