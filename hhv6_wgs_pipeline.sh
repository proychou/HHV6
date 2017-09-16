#!/bin/bash
#Pipeline for whole genome HHV6
#August 2017
#Pavitra Roychoudhury

#Usage: 
#First build reference for bowtie and make a copy of the ref seqs:
# 	module load bowtie2
# 	bowtie2-build './NC_000898.1.fasta' hhv6B_ref_z29
#   bowtie2-build './NC_001664.2.fasta' hhv6A_ref_U1102
# 	cp './NC_001664.2.fasta' hhv6A_ref_U1102.fasta
#  	cp './NC_000898.1.fasta' hhv6B_ref_z29.fasta
#For paired-end library
#		hhv6_wgs_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz
#For single-end library
#		hhv6_wgs_pipeline.sh -s yourreads.fastq.gz
#This is meant to be run on the cluster (typically through sbatch) so if run locally,
#first set the environment variable manually, e.g.
#		SLURM_CPUS_PER_TASK=8
#or whatever is the number of processors available for the process

#Load required tools
#Note that samtools, mugsy, spades, bbmap and last are all locally installed and need to be updated manually as required
# module load bowtie2
# module load FastQC/0.11.5-Java-1.8.0_92
# module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0-fh1
# module load prokka/1.11-foss-2016b-BioPerl-1.7.0

#To do: replace local version of samtools with module load since these have now caught up

PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.9.0-Linux/bin:$HOME/mugsy_x86-64-v1r2.2:$HOME/last759/:$HOME/bbmap/:$HOME/samtools-1.3.1/:
export MUGSY_INSTALL=$HOME/mugsy_x86-64-v1r2.2
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK
# echo "Path: "$PATH

while getopts ":1:2:s:f" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		s) in_fastq="$OPTARG"
			paired="false"
		;;
		f) filter="true"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@

#For testing single-end
# in_fastq='/fh/fast/jerome_k/HHV6_PR/fastq_files/2017_08_10//ABI-HHV6A_S385_L001_R1_001.fastq.gz' 
# paired="false"
# filter="true"
#For testing paired-end
# in_fastq_r1='/fh/fast/jerome_k/HHV6_PR/fastq_files/2017_08_10/HHV6-MT4-PFAR1_S137_L001_R1_001.fastq.gz'
# in_fastq_r2='/fh/fast/jerome_k/HHV6_PR/fastq_files/2017_08_10/HHV6-MT4-PFAR1_S137_L001_R2_001.fastq.gz'
# paired="true"


##  PAIRED-END  ##
if [[ $paired == "true" ]]
then
if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq_r1%%_R1_001.fastq*})

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2  out1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'  out1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' './trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' out1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' out2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20


#Use bbduk to filter reads that match HHV6 genomes
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hhv6_refs.fasta ... \n\n\n"
mkdir -p ./filtered_fastq/

bbduk.sh in1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' in2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' out1='./filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' out2='./filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' outm1='./filtered_fastq/'$sampname'_matched_hhv6_r1.fastq.gz' outm2='./filtered_fastq/'$sampname'_matched_hhv6_r2.fastq.gz' ref='./hhv6_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hhv6.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK

rm './filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' './filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' 
mv './filtered_fastq/'$sampname'_matched_hhv6_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz'
mv './filtered_fastq/'$sampname'_matched_hhv6_r2.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz'
fi

#FastQC report on processed reads
mkdir -p ./fastqc_reports_preprocessed
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
fastqc './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o ./fastqc_reports_preprocessed

#Map reads to reference
printf "\n\nMapping reads to reference seqs hhv6A_ref_U1102 and hhv6B_ref_z29 ... \n\n\n"
mkdir -p ./mapped_reads
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bowtie2 -x $ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes 
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

#Delete some spades folders to free up space
rm -r './contigs/'$sampname'/corrected' 



##  SINGLE-END  ##
else 
if [[ $paired == "false" ]]
then
if [ -z $in_fastq ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq%%_R1_001.fastq*})

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq -o ./fastqc_reports_raw

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in=$in_fastq out='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'  out='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' 

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' out='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#Use bbduk to filter reads that match HHV6 genomes
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hhv6_refs.fasta ... \n\n\n"
mkdir -p ./filtered_fastq/
bbduk.sh in='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' out='./filtered_fastq/'$sampname'_unmatched.fastq.gz' outm='./filtered_fastq/'$sampname'_matched_hhv6.fastq.gz' ref='./hhv6_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hhv6.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK
rm './filtered_fastq/'$sampname'_unmatched.fastq.gz' 
mv './filtered_fastq/'$sampname'_matched_hhv6.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz'
fi

#FastQC report on processed reads
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
mkdir -p ./fastqc_reports_preprocessed
fastqc './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o ./fastqc_reports_preprocessed

#Map reads to reference
printf "\n\nMapping reads to reference seqs hhv6A_ref_U1102 and hhv6B_ref_z29 ... \n\n\n"
mkdir -p ./mapped_reads
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bowtie2 -x $ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -s './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

fi
fi


#Generate sorted bams for mapped reads
printf "\n\nMaking and sorting bam files ... \n\n\n"
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
if [ -f './mapped_reads/'$sampname'_'$ref'.sam' ]
then
~/samtools-1.3.1/samtools view -bh -o './mapped_reads/'$sampname'_'$ref'.bam' './mapped_reads/'$sampname'_'$ref'.sam' -T $ref'.fasta'  
rm './mapped_reads/'$sampname'_'$ref'.sam'
~/samtools-1.3.1/samtools sort -o './mapped_reads/'$sampname'_'$ref'.sorted.bam' './mapped_reads/'$sampname'_'$ref'.bam' 
rm './mapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'Mapping to '$ref 'failed. No sam file found'
fi
done


#Map contigs to refs
printf "\n\nMapping scaffolds to reference seqs hhv6A_ref_U1102 hhv6B_ref_z29 ... \n\n\n"
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
mugsy --directory `readlink -f './contigs/'$sampname` --prefix 'aligned_scaffolds_'$ref $ref'.fasta' `readlink -f './contigs/'$sampname'/scaffolds.fasta'`
sed '/^a score=0/,$d' './contigs/'$sampname'/aligned_scaffolds_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf'
python ~/last-759/scripts/maf-convert sam -d './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
~/samtools-1.3.1/samtools view -bS -T $ref'.fasta' './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam' | ~/samtools-1.3.1/samtools sort > './contigs/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
rm './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
done
rm *.mugsy.log


#To do (maybe): replace mugsy step with scaffold_builder if that is better

#Make new reference sequence using scaffolds
printf "\n\nMaking a reference sequence for remapping ... \n\n\n"
mkdir -p ./ref_for_remapping
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bamfname='./contigs/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
reffname=$ref'.fasta'
Rscript --vanilla hhv6_make_reference.R bamfname=\"$bamfname\" reffname=\"$reffname\" 
done


#Remap reads to "new" reference
printf "\n\nRe-mapping reads to assembled sequence ... \n\n\n"
mkdir -p ./remapped_reads

if [[ $paired == "false" ]]
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bowtie2-build './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './remapped_reads/'$sampname'_'$ref'.sam'
done
fi
 
if [[ $paired == "true" ]]
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bowtie2-build './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './remapped_reads/'$sampname'_'$ref'.sam'
done
fi

#Make sorted bam
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
if [ -f './remapped_reads/'$sampname'_'$ref'.sam' ]
then
~/samtools-1.3.1/samtools view -bh -o './remapped_reads/'$sampname'_'$ref'.bam' './remapped_reads/'$sampname'_'$ref'.sam' -T './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta'
rm './remapped_reads/'$sampname'_'$ref'.sam'
~/samtools-1.3.1/samtools sort -o './remapped_reads/'$sampname'_'$ref'.sorted.bam' './remapped_reads/'$sampname'_'$ref'.bam'
rm './remapped_reads/'$sampname'_'$ref'.bam'
else
echo 'No sam file found'
fi
done

#Call R script to merge bams and generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs_all
mkdir -p ./stats
if [[ $paired == "true" ]]
then
Rscript --vanilla hhv6_generate_consensus.R s1=\"$in_fastq_r1\"
else
if [[ $paired == "false" ]]
then
Rscript --vanilla hhv6_generate_consensus.R s1=\"$in_fastq\"
fi
fi

#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka_6A
prokka --outdir './annotations_prokka_6A/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 6A' --species '' --proteins HHV6_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_6A/'$sampname/*.fa
mkdir -p ./annotations_prokka_6B
prokka --outdir './annotations_prokka_6B/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 6B' --species '' --proteins HHV6_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_6B/'$sampname/*.fa
