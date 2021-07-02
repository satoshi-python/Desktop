#! /lustre7/home/lustre3/satoshi/ALL_PCA/python/bash/pca.sh
#$ -e /lustre7/home/lustre3/satoshi/ALL_PCA/python/bash/error
#$ -o /lustre7/home/lustre3/satoshi/ALL_PCA/python/bash/out
#$ -S /usr/bin/sh
#$ -l s_vmem=512G,mem_req=512G

PATH=~/anaconda3/bin:$PATH
export OMP_NUM_THREADS=1
export TOOL=TOOL=/lustre7/home/lustre3/satoshi/sa/ALL/python
export PYTHONPATH=$TOOL:$PYTHONPATH

cd /lustre7/home/lustre3/satoshi/ALL_PCA/txt
prob_path="/lustre7/home/lustre3/satoshi/MED/aff4/prob.txt /lustre7/home/lustre3/satoshi/MED/eaf1/prob.txt /lustre7/home/lustre3/satoshi/MED/taf7/prob.txt /lustre7/home/lustre3/satoshi/MED/aff4_kai/prob.dat /lustre7/home/lustre3/satoshi/MED/eaf1_kai/prob.dat /lustre7/home/lustre3/satoshi/MED/taf7_kai/prob.dat"
trr_path="/lustre7/home/lustre3/satoshi/MED/aff4/test_all.trr /lustre7/home/lustre3/satoshi/MED/eaf1/test_all.trr /lustre7/home/lustre3/satoshi/MED/taf7/test_all.trr /lustre7/home/lustre3/satoshi/MED/aff4_kai/run_all.trr /lustre7/home/lustre3/satoshi/MED/eaf1_kai/run_all.trr /lustre7/home/lustre3/satoshi/MED/taf7_kai/run_all.trr"
pdb_path="/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb /lustre7/home/lustre3/satoshi/MED/eaf1/HEN.pdb /lustre7/home/lustre3/satoshi/MED/taf7/HEN.pdb /lustre7/home/lustre3/satoshi/MED/aff4_kai/aff4kai.pdb /lustre7/home/lustre3/satoshi/MED/eaf1_kai/eaf1kai.pdb /lustre7/home/lustre3/satoshi/MED/taf7_kai/taf7kai.pdb"
select="4 CA O"
sava_path="aff_20201211.txt eaf_20201211.txt taf_202001211.txt affkai_20201211.txt eafkai_20201211.txt tafkai_20201211.txt"
rm kai.txt
python /lustre7/home/lustre3/satoshi/ALL_PCA/python/pca1.py $trr_path $pdb_path $prob_path $select > kai.txt
python /lustre7/home/lustre3/satoshi/ALL_PCA/python/pca2.py kai.txt $prob_path  $save_path 
rm kai.txt
