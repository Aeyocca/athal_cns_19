#!/bin/bash
#blast_to_dist_vect.sh
#Alan E. Yocca
#06-03-19
#read in blast file on cmd line, generate distance vector file
#for proximate feature
#super specific,
#taking TAIR_te blast file to genome
#calculating proximate distance to collinear and PosV CNS

if [ -z "$1" ] ; then
        echo "must provide input"
        exit
fi

#blast to bed
#AT1TE52125|-|15827287|15838845|ATHILA2|LTR/Gypsy|11559  Chr1    97.07   6171    147     7       5399    11559   15237051        15230905        0.0     1.079e+04
#SRR1946047_bwa_alt_3_pbj_5_rename_tair_te_e4.blast
TAG=${TAG:-"_tair_te_e4.blast"}
BASE=$(basename $1 | sed "s/\.blast//")
SRR=$(basename $1 | sed "s/${TAG}//")
WKDIR="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/26_repeat/"

cd ${WKDIR}

module load BEDTools/2.27.1

echo "Using tag; ${TAG}"

echo "Removing existing outputs"
rm -f ${WKDIR}/01_bed/${BASE}.bed
rm -f ${WKDIR}/01_bed/${BASE}_sort.bed
rm -f ${WKDIR}/02_closest/${BASE}_sort_mx_moved_closest.txt
rm -f ${WKDIR}/02_closest/${BASE}_sort_mx_coll_closest.txt

echo "Starting ${SRR}"; date

awk -v OFS="	" -F "\t" \
'{if ($10 >= $9) \
{print $2,$9,$10,$1"_"$9"_"$10"_"$7"_"$8,"+","500"} \
else print $2,$10,$9,$1"_"$10"_"$9"_"$7"_"$8,"-","500"}' \
$1 > ${WKDIR}/01_bed/${BASE}.bed

#sort
bedtools sort -i ${WKDIR}/01_bed/${BASE}.bed \
> ${WKDIR}/01_bed/${BASE}_sort.bed

#echo "sort done"

#closest moved
bedtools closest -d -t all -nonamecheck \
-a /mnt/research/edgerpat_lab/AlanY/01_athal_cns/cns_calls/acc_moved/${SRR}_mx_moved.bed \
-b ${WKDIR}/01_bed/${BASE}_sort.bed \
> ${WKDIR}/02_closest/${BASE}_sort_mx_moved_closest.txt

#echo "closest moved done"

#closest coll
bedtools closest -d -t all -nonamecheck \
-a /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/17_ATAC/02_bed/01_cns/${SRR}_mx_coll.bed \
-b ${WKDIR}/01_bed/${BASE}_sort.bed \
> ${WKDIR}/02_closest/${BASE}_sort_mx_coll_closest.txt

#dist_vect moved
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/prox_feature_distance.pl \
--input ${WKDIR}/02_closest/${BASE}_sort_mx_moved_closest.txt \
--output ${WKDIR}/03_dist_vect/${BASE}_sort_mx_moved_closest_dist_vect.txt

#dist_vect coll
/mnt/research/edgerpat_lab/AlanY/01_athal_cns/perl/prox_feature_distance.pl \
--input	${WKDIR}/02_closest/${BASE}_sort_mx_coll_closest.txt \
--output ${WKDIR}/03_dist_vect/${BASE}_sort_mx_coll_closest_dist_vect.txt

echo "Finished ${SRR}"; date; echo ""
