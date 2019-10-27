# HHV6 script: This script imports bam files and makes a consensus sequence
# Pavitra Roychoudhury
# Aug 2017

# Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(ShortRead);
library(Biostrings);
library(RCurl);

#Get latest stable version of wgs_functions.R from github
# source('./wgs_functions.R'); #or locally
script<-getURL('https://raw.githubusercontent.com/proychou/ViralWGS/master/wgs_functions.R',
							 ssl.verifypeer=FALSE)
eval(parse(text=script));

#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
    print(args[[i]])
  }
}

#For testing (these args should come from command line)
# s1<-'/fh/fast/jerome_k/HHV6_PR/fastq_files/2017_08_10//ABI-HHV6A_S385_L001_R1_001.fastq.gz'

#Files, directories, target site
merged_bam_folder<-'./remapped_reads/'; 
sampname<-strsplit(basename(s1),'_R1_001.fastq*')[[1]][1];
mapped_reads_folder<-'./mapped_reads/';
con_seqs_dir<-'./consensus_seqs_all';

#Make consensus sequences against HHV6A and B--returns TRUE if this worked
conseq<-clean_consensus_hhv6(sampname,merged_bam_folder,mapped_reads_folder);

#Prepare seqs for annotation -- will make separate folders for A and B
if(conseq==TRUE){
	if(!dir.exists('./annotations_prokka_6A')) dir.create('./annotations_prokka_6A');
	if(!dir.exists('./annotations_prokka_6B')) dir.create('./annotations_prokka_6B');
	
	
	#Remove all Ns at the beginning and end of the seq, write to folder
	for(ref in c('_hhv6A_ref_U1102','_hhv6B_ref_z29')){
		fname<-grep(paste0('\\/',sampname,'_'),list.files(con_seqs_dir,ref,full.names=T),value=T);
		con_seq<-readDNAStringSet(fname);
		con_seq_trimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(con_seq))));
		names(con_seq_trimmed)<-substring(names(con_seq),1,20); #prokka needs contig name to be <=20 chars long
		if(ref=='_hhv6A_ref_U1102'){
			sampdir<-paste('./annotations_prokka_6A/',sampname,sep='');
			if(!dir.exists(sampdir)) dir.create(sampdir); #create folder for the sample
			writeXStringSet(con_seq_trimmed,file=paste(sampdir,'/',sampname,'.fa',sep=''),format='fasta');
		}else if(ref=='_hhv6B_ref_z29'){
			sampdir<-paste('./annotations_prokka_6B/',sampname,sep='');
			if(!dir.exists(sampdir)) dir.create(sampdir); #create folder for the sample
			writeXStringSet(con_seq_trimmed,file=paste(sampdir,'/',sampname,'.fa',sep=''),format='fasta');
		}
	}
	
}else{
	print('Failed to generate consensus sequences.')
}

