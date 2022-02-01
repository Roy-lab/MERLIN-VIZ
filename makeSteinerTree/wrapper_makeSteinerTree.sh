INPUT=/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/output_net_0_8_sorted.txt

HEADDIR=/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100/steiner_tree/

PATHWAY=$1

INFILE=${HEADDIR}/input/${PATHWAY}.txt

EXE=/mnt/dv/wid/projects7/Roy-Aspergillus/Programs/makeSteinerTree/makeSteinerTree

mkdir ${HEADDIR}/${PATHWAY}

#while read -r line;
#do
#	echo "$line"
	
#	./makeSteinerTree ${INPUT} ${HEADDIR}/input/${PATHWAY}.txt $line ${HEADDIR}/${PATHWAY}/$line

#done < ${HEADDIR}/input/${PATHWAY}.txt

GENE=`head -1 ${INFILE}`

${EXE} ${INPUT} $INFILE ${HEADDIR}/${PATHWAY}/${PATHWAY}
