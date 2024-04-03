import argparse
import torch
import esm
import os
import time

def main(args):
	model = esm.pretrained.esmfold_v1()
	model = model.eval().cuda()

	prot_name = args.protein_name
	input_fasta = args.input_fasta
	out_path_predictions = args.output_folder

	list_out_all = os.listdir(out_path_predictions)
	if out_path_predictions not in list_out_all:
			os.mkdir(out_path_predictions)

	print('input_fasta:', input_fasta)
	print('folder_out:', out_path_predictions)

	with open(input_fasta) as fasta:
		lines = fasta.readlines()
	seqs = [lines[line].strip('\n') for line in range(len(lines)) if line % 2 != 0]

	numb_models = len(seqs)
	count = 1
	start = time.time()

	for index in range(len(seqs)):
		with torch.no_grad():
			output = model.infer_pdb(seqs[index])
		with open(out_path_predictions + '/' + '{}_{}.pdb'.format(index, prot_name), 'w') as f:
			f.write(output)
			print('predicted', count, '/', numb_models, 'structures')
			count += 1

	end = time.time()
	print('total execution time:', round(end - start, 2) / 60)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="ESMfold structure prediction")
	parser.add_argument("--protein_name", type=str, help="Name of the protein")
	parser.add_argument("--input_fasta", type=str, help="Path to the input fasta file")
	parser.add_argument("--output_folder", type=str, help="Path to the output folder")
	args = parser.parse_args()
	main(args)