#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --gpus=1
#SBATCH --partition=gpu_mig
#SBATCH --time=01:00:00
#SBATCH --mem=20G

module load 2023
module load Anaconda3/2023.07-2
conda init bash
source activate ProteinMPNN # initialize your conda environment on the cluster or your local machine

folder_pdbs="<>" # <> means put the path to the folder containing the pdbs
out_dir="<>"
out_seqs=$out_dir'/<>'

if [ ! -d $out_dir ]
then
        mkdir -p $out_dir
fi

if [ ! -d $out_seqs ]
then
        mkdir -p $out_seqs
fi

path_parsed_chains=$out_seqs"/parsed_pdbs.jsonl"
path_for_fixed_positions=$out_seqs"/fixed_pdbs.jsonl"
path_for_assigned_chains=$out_seqs"/assigned_pdbs.jsonl"
chains_to_design="A"
fixed_positions="8 9 15 40 41 42 43 44 45 46 47 48 49 51 52 61 79 80 100 101 102 103 106 107 109 110 111 112 113 114 115 116 117 118 158 159 160 161 163 164 165 166 280 281 285 287 289 291 319 321 371 372 373 374 208 248 249 273 274 317 342 344"

python <>/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_pdbs --output_path=$path_parsed_chains
python <>/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python <>/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python <>/ProteinMPNN/protein_mpnn_run.py \
			--use_soluble_model \ 
			--jsonl_path $path_parsed_chains \
			--chain_id_jsonl $path_for_assigned_chains \
			--fixed_positions_jsonl $path_for_fixed_positions \
			--out_folder $out_seqs \
			--num_seq_per_target 150 \
			--sampling_temp "0.1" \
			--seed 159 \
			--batch_size 1