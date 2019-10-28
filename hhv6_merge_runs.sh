#!/bin/bash
#Rerunning Spades on samples that need to be merged
cd /fh/fast/jerome_k/HHV6_WGS/merging_runs/
module load bowtie2/2.2.5
module load R/3.6.0-foss-2016b-fh1
module load prokka/1.13-foss-2016b-BioPerl-1.7.0
module load parallel/20170222-foss-2016b
module load Perl/5.28.0-foss-2016b
cp -R /fh/fast/jerome_k/HHV6_WGS/refs/ ./


PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.11.1-Linux/bin:$HOME/mugsy_x86-64-v1r2.2:$HOME/last759/:$HOME/bbmap/:$HOME/samtools-1.3.1/:
export MUGSY_INSTALL=$HOME/mugsy_x86-64-v1r2.2
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK

while getopts "s:c:b:" opt; do
	case $opt in
		s) sampname="$OPTARG"
		;;
		c) spades_command="$OPTARG"
		;;
		b) bowtie_command="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@

#Assemble with SPAdes 
mkdir -p '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname
eval $spades_command' -t '$SLURM_CPUS_PER_TASK

#Delete some spades folders to free up space
rm -r '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/corrected' 


#Map contigs to refs
printf "\n\nMapping scaffolds to reference seqs hhv6A_ref_U1102 hhv6B_ref_z29 ... \n\n\n"
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
mugsy --directory `readlink -f '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname` --prefix 'aligned_scaffolds_'$ref ./refs/$ref'.fasta' `readlink -f '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/scaffolds.fasta'`
sed '/^a score=0/,$d' '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_'$ref'.maf' > '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf'
python ~/last-759/scripts/maf-convert sam -d '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf' > '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_'$ref'.sam'
~/samtools-1.3.1/samtools view -bS -T ./refs/$ref'.fasta' '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_'$ref'.sam' | ~/samtools-1.3.1/samtools sort > '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
rm '/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/aligned_scaffolds_'$ref'.sam'
done
rm *.mugsy.log



#Make new reference sequence using scaffolds---need to figure out why part of the seq is missing
printf "\n\nMaking a reference sequence for remapping ... \n\n\n"
mkdir -p ./ref_for_remapping
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bamfname='/fh/fast/jerome_k/HHV6_WGS/merging_runs/assembly/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
reffname=./refs/$ref'.fasta'
Rscript --vanilla /fh/fast/jerome_k/HHV6_WGS/hhv6_make_reference.R bamfname=\"$bamfname\" reffname=\"$reffname\" 
done

#Remap reads to "new" reference
printf "\n\nRe-mapping reads to assembled sequence ... \n\n\n"
mkdir -p ./remapped_reads
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
bowtie2-build './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
eval $bowtie_command' -p '$SLURM_CPUS_PER_TASK' -S ./remapped_reads/'$sampname'_'$ref'.sam -x ./ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
done
rm ./ref_for_remapping/$sampname*.bt2


#Make sorted bam
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
if [ -f './remapped_reads/'$sampname'_'$ref'.sam' ]
then
~/samtools-1.3.1/samtools view -bh -o './remapped_reads/'$sampname'_'$ref'.bam' './remapped_reads/'$sampname'_'$ref'.sam' -T './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta'
rm './remapped_reads/'$sampname'_'$ref'.sam'
~/samtools-1.3.1/samtools sort -o './remapped_reads/'$sampname'_'$ref'.sorted.bam' './remapped_reads/'$sampname'_'$ref'.bam'
rm './remapped_reads/'$sampname'_'$ref'.bam'
mv './remapped_reads/'$sampname'_'$ref'.sorted.bam'  './remapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'No sam file found'
fi
done
rm ./ref_for_remapping/$sampname*.fai

#And also map reads to reference sequences for vcf generation
printf "\n\nMapping reads to reference seqs hhv6A_ref_U1102 hhv6B_ref_z29 ... \n\n\n"
mkdir -p ./mapped_reads
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
eval $bowtie_command' -p '$SLURM_CPUS_PER_TASK' -S ./mapped_reads/'$sampname'_'$ref'.sam -x ./refs/'$ref
done

#Generate sorted bams for mapped reads
printf "\n\nMaking and sorting bam files ... \n\n\n"
for ref in hhv6A_ref_U1102 hhv6B_ref_z29; do
if [ -f './mapped_reads/'$sampname'_'$ref'.sam' ]
then
~/samtools-1.3.1/samtools view -bh -o './mapped_reads/'$sampname'_'$ref'.bam' './mapped_reads/'$sampname'_'$ref'.sam' -T './refs/'$ref'.fasta'  
rm './mapped_reads/'$sampname'_'$ref'.sam'
~/samtools-1.3.1/samtools sort -o './mapped_reads/'$sampname'_'$ref'.sorted.bam' './mapped_reads/'$sampname'_'$ref'.bam' 
rm './mapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'Mapping to '$ref 'failed. No sam file found'
fi
done
 
#Call R script to generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs_all
mkdir -p ./stats
Rscript --vanilla '/fh/fast/jerome_k/HHV6_WGS/hhv6_generate_consensus.R' s1=\"$sampname\"


#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka_6A
prokka --outdir './annotations_prokka_6A/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 6A' --species '' --proteins ./refs/HHV6_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_6A/'$sampname/*.fa
mkdir -p ./annotations_prokka_6B
prokka --outdir './annotations_prokka_6B/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 6B' --species '' --proteins ./refs/HHV6_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_6B/'$sampname/*.fa

