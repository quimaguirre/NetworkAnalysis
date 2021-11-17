#!/bin/bash
#SBATCH --mem=50000
#SBATCH -p long
#SBATCH -o /users/sbi/quim/BIANA_scripts/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020.out
#SBATCH -e /users/sbi/quim/BIANA_scripts/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020.err
module load Python/2.7.11
python /users/sbi/quim/BIANA_scripts/scripts/generate_network_interactomix.py -radius 0 -taxid 9606 -edge /users/sbi/quim/BIANA_scripts/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020.txt -node /users/sbi/quim/BIANA_scripts/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020_nodes.txt -trans /users/sbi/quim/BIANA_scripts/outputs/BIANA_2020_geneID_seqtax_drugtarget/human_network_biana_2020_translation.txt -ttype geneid -format multi-fields