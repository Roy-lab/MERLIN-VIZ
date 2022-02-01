#vary this and see what effect it has on the output

#debbies python script
UTIL=/mnt/dv/wid/projects7/Roy-Aspergillus/Programs/aspergillusNetworkDiffusion/get_kernel_scores.py

gene_index=/mnt/dv/wid/projects7/Roy-Aspergillus/Programs/aspergillusNetworkDiffusion/network_files/aspergillus_network_I02_node_index.txt
edge_list=/mnt/dv/wid/projects7/Roy-Aspergillus/Programs/aspergillusNetworkDiffusion/network_files/aspergillus_network_I02_edge_list.txt

lambda=$1
original_scores=$2
outfile=$3

#check if scores exists
scores=nodeDiffusion_Input.txt

if [ -f $scores ];
then
rm $scores
fi


cat $gene_index | while IFS="$(echo t | tr t ' ')" read -r id gene
do
res1=$(grep -P "^${gene}\t" $original_scores | cut -f2)
if [ -z "$res1" ]
then
echo -e ${id}"\t"0 >> ${scores}
else
echo -e ${id}"\t"${res1} >> ${scores}
fi
done
cut -f1 -d' ' $scores > tmp.txt
mv tmp.txt $scores


#merge the scores and column one of the gene_index to get the format that is required by the node diffusion script

mkdir -p prep_node_diffusion
sif_format=prep_node_diffusion/influence_lr_top5_b0.3_edge_list.sif
awk -v OFS='\t' '{print $1,"binding",$2}' $edge_list > $sif_format


kernel=prep_node_diffusion/influence_lr_top5_b03_edge_list_kernel_${lambda}.npy
rm $kernel
if [ -f $kernel ];
then
#kernel exists
python2.7 $UTIL $scores $sif_format $lambda $kernel $outfile
else
python2.7 $UTIL $scores $sif_format $lambda $outfile
fi

rm $scores
rm -r prep_node_diffusion
