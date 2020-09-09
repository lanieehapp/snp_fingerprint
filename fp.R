##fingerprinting

library(stringr)

load("fingerprint_snps.RData")

ref="GRCh38.p12.genome.plus.ERCC.fa"
ref_fai="GRCh38.p12.genome.plus.ERCC.fa.fai"
bam="gatk_filtersamreads__nonspliced__dna_GATK_FilterSamReads.21198_T_1.bam"

#comm<-paste0('parallel --colsep "\t" samtools mpileup -a -l fingerprinting.bed --fasta-ref ',ref, ' ', bam, ' -r {1} :::: ', ref_fai ,' > fingerprint.txt' )

comm<-paste0('samtools mpileup -a -l fingerprinting.bed --fasta-ref ',ref, ' ', bam,  '> fingerprint.txt' )

system(comm)

fp<-read.csv(file="fingerprint.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE, quote = "", fill=FALSE)

colnames(fp)<-c("Chr", "Pos", "Ref", "Depth", "mpileup", "Qual")
fp$CUR_POS<-paste(fp$Chr, fp$Pos, sep="-")

fingerprint_snps$CHROM_POS_REF_ALT<-paste(fingerprint_snps$Chr, fingerprint_snps$Start, fingerprint_snps$Ref, fingerprint_snps$Alt, sep="-")

fp$ref_FS<-str_count(fp$mpileup, pattern="[.]")
fp$ref_RS<-str_count(fp$mpileup, pattern=",")
fp$A_FS<-str_count(fp$mpileup, pattern="A")
fp$A_RS<-str_count(fp$mpileup, pattern="a")
fp$C_FS<-str_count(fp$mpileup, pattern="C")
fp$C_RS<-str_count(fp$mpileup, pattern="c")
fp$G_FS<-str_count(fp$mpileup, pattern="G")
fp$G_RS<-str_count(fp$mpileup, pattern="g")
fp$T_FS<-str_count(fp$mpileup, pattern="T")
fp$T_RS<-str_count(fp$mpileup, pattern="t")
fp$DEL<-str_count(fp$mpileup, pattern="[*]")
fp$INS<-str_count(fp$mpileup, pattern="[+]")

fp$A_evidence<-fp$A_FS>0 & fp$A_RS>0
fp$C_evidence<-fp$C_FS>0 & fp$C_RS>0
fp$T_evidence<-fp$T_FS>0 & fp$T_RS>0
fp$G_evidence<-fp$G_FS>0 & fp$G_RS>0
fp$INS_evidence<-fp$INS>0
fp$DEL_evidence<-fp$DEL>0

nt<-c("A", "C", "T", "G")
all_rna<-NULL
for(i in 1:nrow(fingerprint_snps)){
  curr.pos<-paste(fingerprint_snps$Chr[i], fingerprint_snps$Start[i], sep="-")
  curr.ref<-fingerprint_snps$Ref[i]
  curr.alt<-fingerprint_snps$Alt[i]
  tmp<-fp[fp$CUR_POS == curr.pos,]
  
  RNA_depth_total<-tmp$Depth
  if(RNA_depth_total>=30){
    RNA_evidence="Covered"
  } else if(RNA_depth_total < 2){
    RNA_evidence="No coverage"
  } else if(RNA_depth_total<30){
    RNA_evidence<-"Low coverage"
  }
  
  #processing for SNVs
  if(curr.alt %in% nt & curr.ref %in% nt){
    if(sum(tmp$A_evidence, tmp$T_evidence, tmp$C_evidence, tmp$G_evidence, tmp$INS_evidence, tmp$DEL_evidence)>1){
      RNA_evidence<-paste(RNA_evidence, "Multiallelic Locus", sep=";")
    }
    
    if(curr.alt == "A" & tmp$A_evidence){
      RNA_depth_alt<-tmp$A_FS+tmp$A_RS
      RNA_AF<-RNA_depth_alt/RNA_depth_total
      RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
    } else if(curr.alt == "C" & tmp$C_evidence){
      RNA_depth_alt<-tmp$C_FS+tmp$C_RS
      RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
      RNA_AF<-RNA_depth_alt/RNA_depth_total
    } else if(curr.alt == "G" & tmp$G_evidence){
      RNA_depth_alt<-tmp$G_FS+tmp$G_RS
      RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
      RNA_AF<-RNA_depth_alt/RNA_depth_total
    } else if(curr.alt == "T" & tmp$T_evidence){
      RNA_depth_alt<-tmp$T_FS+tmp$T_RS
      RNA_evidence<-paste(RNA_evidence, "TRUE", sep=";")
      RNA_AF<-RNA_depth_alt/RNA_depth_total
    } else {
      RNA_depth_alt=0
      RNA_evidence<-paste(RNA_evidence, "FALSE", sep=";")
      RNA_AF<-0
    }
    
  }
  
  
  
  
  all_rna<-rbind(all_rna, cbind(RNA_depth_total, RNA_depth_alt, RNA_AF, RNA_evidence))
  
}

colnames(all_rna)<-paste0(samp, ".", colnames(all_rna))
return(all_rna)



