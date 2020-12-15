#!/bin/sh -login
#Alan E Yocca
#02-14-18
#atac_callpeak

#SBATCH --job-name=atac_callpeak
#SBATCH --mem=10Gb
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/research/edgerpat_lab/AlanY/Error_files/%x-%j.SLURMout

#cd /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/

if [ ! -z "$1" ]; then
        if [ $1 == "-h" ] || [ $1 == "--help" ]; then
                echo -e "\
                ############ Usage ##############\n\
                #                               #\n\
                #       BAM		        #\n\
                #                               #\n\
                #       OUTPUT                  #\n\
                #         output sam file       #\n\
                #                               #\n\
                #       FORCE                   #\n\
                #         overwrite output      #\n\
                #         if exists             #\n\
                #################################\n\
                "
                exit
        fi
fi

FORCE=${FORCE:-false}

if [ ${FORCE} != "false" ]; then
	rm -f ${OUTPUT}*
fi

if [ -z ${BAM+x} ] || [ -z ${OUTPUT} ]; then
	echo "required arguments not provided"
	exit 1
fi

source activate Python2

macs2 callpeak \
-f BAMPE \
-t ${BAM} \
-g 1.35e8 \
--keep-dup 1 \
-n ${OUTPUT}

echo -e "\nFINISHED"; date

#-keep-dup 1
##have macs2 keep 1 PCR duplicate
##since im not removing duplicates myself

#-f BAMPE
##my reads are paired-end, in bam format

#-g 1.35e8
##estimated genome size for athaliana

