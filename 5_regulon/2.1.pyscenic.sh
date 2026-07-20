#!/usr/bin/env bash
#SBATCH --job-name=scenic11
#SBATCH --account=fc_williamslab
#SBATCH --partition=savio3_bigmem
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ken_chen@berkeley.edu
export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles 
cd /global/scratch/users/kenchen/LAM/5_regulon/reg_results/pyscenic_runs
module load python
source activate regulon_3.9
pyscenic grn lam_ctrl10.loom pyscenic_files/allTFs_hg38.txt -o adj_1.tsv --num_workers 32

pyscenic ctx adj_1.tsv \
    pyscenic_files/motifs/rankings/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    pyscenic_files/motifs/rankings/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname pyscenic_files/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname lam_nonim10.loom \
    --output reg_1.csv \
    --mask_dropouts \
    --num_workers 32
    
pyscenic aucell \
    lam_ctrl10.loom \
    reg_1.csv \
    --output output_1.loom\
    --num_workers 32