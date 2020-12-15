#!/bin/sh -login
#Alan E Yocca
#05-07-18
#Full pipeline for intraspecific Arabidopsis thaliana CNS comparison, for assembly_golicz16
#SPECIFY THINGS ON SUBMISSION LINE (genome and accession)

#! /usr/bin/env bash
#PBS -N at_cns01_comp_pipe_cmd
#PBS -l walltime=5:00:00,nodes=40:ppn=1,mem=15gb
#PBS -j oe
#PBS -o /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -e /mnt/research/edgerpat_lab/AlanY/Error_files
#PBS -m n

#full pipe from genomes and GFF to CNS comparisons
cd /mnt/research/edgerpat_lab/AlanY/01_athal_cns/

echo Started:
date

#LINE=`/bin/sed -n ${PBS_ARRAYID}p ${INFILE}`
#leave this in so don't have to replace ${ACC} everywhere
#add bwa_all_tair10

#specify these upon submission
#ACC="Bur_0_SHORE_hq"
#GENOME="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/05_reference/Bur-0.SHORE.scaffolds.2010-09-30.500bp.fa"

if [ -z "$ACC" ] || [ -z "$GENOME" ]; then
	echo "Forgot to specify accession or genome when submitting, if want to hardcode or array job, try at_cns01_comp_pipe.sh / at_cns01_comp_pipe_assembly_golicz16.sh"
	exit 1;
fi

module load BLAST

formatdb -p F -o T -n "blast/database/${ACC}_db" -i ${GENOME}
echo Starting Blast:
date
blastall -p blastn -i /mnt/home/yoccaala/04_Edger/athal_cns/CNS_Brassicaceae_at_gene_nodups_mi.fa -d "blast/database/${ACC}_db" -e 10 -m 8 -b 5 -v 5 -o blast/${ACC}.blast -a 39
echo Blast finished:
date

perl/filter_cns_acc_blast_60_bs_chrom.pl -b blast/${ACC}.blast -c /mnt/home/yoccaala/04_Edger/athal_cns/CNS_Brassicaceae_at_gene_nodups_mi.fa -o blast/${ACC}_60_bs_chr.blast
perl/last_add_ref_metainfo.pl -i blast/${ACC}_60_bs_chr.blast -o blast/${ACC}_60_bs_chr_mi.blast
cp blast/${ACC}_60_bs_chr_mi.blast ~/01_VanBuren/mcscan/MCScanX/data/
perl/blast_to_gff_mcscanx.pl -b blast/${ACC}_60_bs_chr_mi.blast -o  MCScanX/data/${ACC}_60_bs_chr_mi.gff
cp MCScanX/data/${ACC}_60_bs_chr_mi.gff ~/01_VanBuren/mcscan/MCScanX/data/

~/01_VanBuren/mcscan/MCScanX/MCScanX ~/01_VanBuren/mcscan/MCScanX/data/${ACC}_60_bs_chr_mi

cp ~/01_VanBuren/mcscan/MCScanX/data/${ACC}_60_bs_chr_mi.collinearity MCScanX/data/

rm -f ~/01_VanBuren/mcscan/MCScanX/data/${ACC}_60_bs_chr_mi.blast
rm -f ~/01_VanBuren/mcscan/MCScanX/data/${ACC}_60_bs_chr_mi.gff
rm -f ~/01_VanBuren/mcscan/MCScanX/data/${ACC}_60_bs_chr_mi.collinearity

rm -rf cns_calls/acc_hits/${ACC}_hits.txt
rm -rf cns_calls/acc_missing/${ACC}_missing.txt
rm -rf cns_calls/acc_coll/${ACC}_mx.txt
rm -rf cns_calls/acc_moved/${ACC}_moved.txt

grep -v "^#" blast/${ACC}_60_bs_chr_mi.blast | cut -f7 -d"|" | sort | uniq | cat > cns_calls/acc_hits/${ACC}_hits.txt
comm -23 cns_calls/col_0_cns01_hits.txt cns_calls/acc_hits/${ACC}_hits.txt > cns_calls/acc_missing/${ACC}_missing.txt
grep -v "^#" MCScanX/data/${ACC}_60_bs_chr_mi.collinearity | cut -f21 -d"|" | sort | uniq | cat > cns_calls/acc_coll/${ACC}_mx.txt
comm -23 cns_calls/acc_hits/${ACC}_hits.txt cns_calls/acc_coll/${ACC}_mx.txt > cns_calls/acc_moved/${ACC}_moved.txt

echo Finished:
date

qstat -f ${PBS_JOBID}
