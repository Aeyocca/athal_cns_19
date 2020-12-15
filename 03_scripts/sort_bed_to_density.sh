#!/bin/bash
#sort_bed_to_density.sh
#Alan E. Yocca
#06-05-19
#read in sorted bed file on cmd line, generate density file

if [ -z "$1" ]; then
        echo "must provide input"
        exit
fi

#SRR1946047_bwa_alt_3_pbj_5_rename_tair_te_e4_sort.bed
TAG=${TAG:-"_tair_te_e4_sort.bed"}
BASE=$(basename $1 | sed "s/\.bed//")
SRR=$(basename $1 | sed "s/${TAG}//")
WKDIR="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/26_repeat/"
GENOMEDIR="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/16_final_genomes/04_all/"

cd ${WKDIR}

module load BEDTools/2.27.1

echo "using tag: ${TAG}"
echo "Removing existing outputs"
rm -f ${WKDIR}/01_bed/${BASE}_gc_bg.bed

echo "Starting ${SRR}"; date

#genome coverage
bedtools genomecov -bg \
-i ${1} \
-g ${GENOMEDIR}/${SRR}.fasta.fai \
> ${WKDIR}/01_bed/${BASE}_gc_bg.bed

/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/bed_density.pl \
--bed_a /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/17_ATAC/02_bed/01_cns/${SRR}_mx_moved.bed \
--genomecov ${WKDIR}/01_bed/${BASE}_gc_bg.bed \
--output 05_density/${BASE}_gc_bg_mx_moved_density.txt

/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/bed_density.pl \
--bed_a /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/17_ATAC/02_bed/01_cns/${SRR}_mx_coll.bed \
--genomecov ${WKDIR}/01_bed/${BASE}_gc_bg.bed \
--output 05_density/${BASE}_gc_bg_mx_coll_density.txt

echo "Finished ${SRR}"; date; echo ""
