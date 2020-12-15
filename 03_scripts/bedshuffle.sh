#!/bin/bash
#Alan E Yocca
#bedshuffle
#10-14-18

#########*************JUST MAKING EXECUTABLE BASH SCRIPT SINCE TORQUE NO MORE RIP*********#############
#! /usr/bin/env bash
#PBS -N bedshuffle
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=5gb
#PBS -j oe
#PBS -o /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -e /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -m n

set -e

echo "STAT ATAC, DOES NOT SHUFFLE THE ORIGINAL ATAC ANNOTATIONS"
#but the code is in here if you wanted to: uncomment necessary lines
#change -b file in intersect to the shuffled atac file

#WKDIR where the cns annotations are
WKDIR="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/17_ATAC/02_bed/01_cns"
cd ${WKDIR}

date

#subsampling:
module load BEDTools/2.27.1
module load SAMtools/1.9
#echo "made it past bedtools"

#touch hello_there
#exit 1;

#make hash linking arabidopsis base to genome full path right? that would be handy
#declare -A line_to_genome=(); while read i; do tmp_array=($(echo ${i} | tr '	' '\n')); line_to_genome["${tmp_array[0]}"]="${tmp_array[1]}"; done < /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/acc_srr_conv_atac.txt
#key is Arabidopsis base, value is full path of SRR genome
#since we are looking throug files with SRR number, going to reverse that
#hmm but need to link srr to arab stuffs... make separate hash? yea, separate hash for srr id to full path to genome and srr to arab
declare -A SRR_to_genome=(); while read i; do tmp_array=($(echo ${i} | tr '	' '\n')); BASE=$(basename ${tmp_array[1]}); IFS='_' read -r -a tmp_srr_one <<< "${BASE}"; SRR_to_genome["${tmp_srr_one[0]}"]="${tmp_array[1]}"; done < /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/acc_srr_conv_atac.txt
declare -A SRR_to_arab=(); while read i; do tmp_array=($(echo ${i} | tr '	' '\n')); BASE=$(basename ${tmp_array[1]}); IFS='_' read -r -a tmp_srr_two <<< "${BASE}"; SRR_to_arab["${tmp_srr_two[0]}"]="${tmp_array[0]}"; done < /mnt/research/edgerpat_lab/AlanY/01_athal_cns/qsub/acc_srr_conv_atac.txt

#for i in "${!SRR_to_genome[@]}"
#do
#  echo "key  : $i"
#  echo "value: ${SRR_to_genome[$i]}"
#done

#for every bed file of cns annotations
#eg: SRR1945758_bwa_alt_3_masked_60_bs_chr_mi.bed
#this is also going to run on the bwa versions of the pbj genomes
#screw it just run it on them and ignore the results

#new thought do array submission so no for loop here, make i LINE like have done before:
#INFILE going to be
#INFILE="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/17_ATAC/02_bed/all_moved_overlap.txt
#Keep as variable defined in the sbatch array script in case want to change
#This file above is a list of the moved cns annotations AND all cns annotations for each accession, example file:
#SRR1945852_bwa_alt_3_pbj_2_rename_mx_all.bed
#SRR1945852_bwa_alt_3_pbj_2_rename_mx_moved.bed
#eg, I have 18 genomes, and 36 files in all_moved_overlap.txt
i=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`

#Extract SRR number, a little clunky but it seems to work
echo "starting on ${i}"
#remove final output if it exists
IFS='_' read -r -a SRR_tmp <<< "${i}"
#####used to have extra step where SRR[1] ends up being what I want
#####That is whats throughout the script so just left it
SRR[1]="${SRR_tmp[0]}"
echo "srr; ${SRR[1]}"
#now ${SRR[1]} is the SRR number
fai_check=$(ls ${SRR_to_genome["${SRR[1]}"]}.fai 2> /dev/null | wc -l)
if [ ${fai_check} != "0" ]; then
	echo "fai exists"
else
	echo "fai does not exist, creating"
	samtools faidx ${SRR_to_genome["${SRR[1]}"]}
fi
#######Add actual overlap, ie cns bed file, non shuffled, overlap with atac peaks, not shuffled
ibase=$(echo ${i} | sed "s/\.bed//")
rm -f ${WKDIR}/${ibase}_stat_overlap_count.txt
#Check if the atac_overlap file exists, if not, create it:
if [ ! -e ${WKDIR}/../../03_overlap/${ibase}_atac_overlap.txt ]; then
	echo "file: ${WKDIR}/../../03_overlap/${ibase}_atac_overlap.txt"; echo "does NOT exist, creating"
	bedtools intersect -nonamecheck -f 0.5 -u -a ${WKDIR}/${i} \
	-b ${WKDIR}/../02_atac_peak/${SRR_to_arab[${SRR[1]}]}_atac_peaks_sort.xls > \
	${WKDIR}/../../03_overlap/${ibase}_atac_overlap.txt
fi
OVERLAP=$(wc -l ${WKDIR}/../../03_overlap/${ibase}_atac_overlap.txt | sed "s/ .*//")
Annotation_count=$(wc -l ${WKDIR}/${i} | sed "s/ .*//")
echo "Annotation count: ${Annotation_count}"
########OVERLAP=$(wc -l ${WKDIR}/../../03_overlap/${i} | sed "s/ .*//")
for n in {1..1000}; do
	if ! (( ${n} % 100 )); then
		echo "Starting permutation ${n}"
		date
	fi
#	echo "srr: ${SRR[1]}"
	rm -f ${WKDIR}/${ibase}_shuffle_${n}.bed
#####	stat atac, not shuffling atac
#	rm -f ${WKDIR}/${ibase}_atac_shuffle_${n}.bed
	rm -f ${WKDIR}/${ibase}_intersect.txt
#	echo "shuffle 1"
	bedtools shuffle -i ${WKDIR}/${i} -g ${SRR_to_genome["${SRR[1]}"]}.fai > ${WKDIR}/${ibase}_shuffle_${n}.bed
#	echo "shuffle 2"
#####	stat atac, not shuffling atac
#	bedtools shuffle -i ${WKDIR}/../02_atac_peak/${SRR_to_arab[${SRR[1]}]}_atac_peaks_sort.xls -g ${SRR_to_genome[${SRR[1]}]}.fai > ${WKDIR}/${ibase}_atac_shuffle_${n}.bed
#	echo "intersect"
#####	stat atac, not shuffling atac
	bedtools intersect -nonamecheck -f 0.5 -u -a ${WKDIR}/${ibase}_shuffle_${n}.bed \
	-b  ${WKDIR}/../02_atac_peak/${SRR_to_arab[${SRR[1]}]}_atac_peaks_sort.xls > \
	${WKDIR}/${ibase}_intersect.txt

	inter_count=$(wc -l ${ibase}_intersect.txt | sed "s/ .*//")
#	echo "srr1: ${SRR[1]}"
#	echo "base: ${BASE}"
#	echo "${n}	${Annotation_count}	${inter_count}	${OVERLAP}"
	echo "${n}	${Annotation_count}	${inter_count}	${OVERLAP}" >> ${WKDIR}/${ibase}_stat_overlap_count.txt
#	exit 1
	rm -f ${WKDIR}/${ibase}_shuffle_${n}.bed
#####	stat atac, not shuffling atac
#	rm -f ${WKDIR}/${ibase}_atac_shuffle_${n}.bed
	rm -f ${WKDIR}/${ibase}_intersect.txt
done

date

#qstat -f ${PBS_JOBID}
