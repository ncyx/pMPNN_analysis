import random
from mpnn_analysis import *

class additional_tools:

	def __init__(self, fasta, folder_out, msa_wt=None, pdb=None, filename=None, str_predict=None):
		self.fasta = fasta
		self.folder_out = folder_out
		self.pdb = pdb
		self.msa_wt = msa_wt
		self.filename = filename
		self.str_predict = str_predict

	def gen_fa_singleseq(self):
		'''
		Generates single-sequence fasta files of the designs from pMPNN from the output FASTA file.
		These files can be used for single-sequence AF2 prediction (batch).

		Parameters:
				None

		Returns:
				None
		'''
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		all_seqs = [lines[i].strip('\n') for i in range(len(lines)) if i %2 != 0]
		count = 0
		for sequence in all_seqs:
			with open(self.folder_out+'/'+'{}_{}.a3m'.format(count, self.filename), 'w') as single_seq:
				single_seq.write('>{}_{}'.format(self.filename, count)+'\n'+str(sequence))
			count+=1

	def swap_msa_seq(self):
		'''
		Replaces the parent sequence in the Multiple Sequence Alignment (MSA) with the designed sequence for subsequent AF2 "msa|swap" prediction.
		
		Parameters:
				self.msa_wt - MSA in .a3m format which has been computed for the native (parent) sequence

		Returns:
				None
		'''
		with open(self.fasta, 'r') as fasta_mpnn:
			lines_fasta = fasta_mpnn.readlines()
		seqs_des = [lines_fasta[index].strip('\n') for index in range(len(lines_fasta)) if index%2 !=0]
		with open(self.msa_wt, 'r') as msa_wt_input:
			lines_msa = msa_wt_input.readlines()
		for index in range(len(seqs_des)): #+1: the wt-msa as well to check and compare
			if index == 0:
				with open(self.folder_out +'/'+'{}_{}.a3m'.format(index, self.filename), 'w') as out_msa_changed:
					for line in range(len(lines_msa)):
						out_msa_changed.write(lines_msa[line])
			else:
				with open(self.folder_out +'/'+'{}_{}.a3m'.format(index, self.filename), 'w') as out_msa_changed:
					for line in range(len(lines_msa)):
						if line == 2:
							out_msa_changed.write(seqs_des[index]+'\n')
						else:
							out_msa_changed.write(lines_msa[line])

	def scramble_design(self, percentage = None):
		'''
		Scrambles the designed sequences by replacing a specified percentage of positions with random residues.

		Parameters:
				percentage (float): The percentage of positions to mutate [0, 1]. If None, all positions differing from the native sequence are mutated.
				
		Returns:
				list: A list containing the scrambled sequences.

		If percentage is None, each position that differs from the original sequence in the designs is mutated by replacing it with a random residue. 
		The mutated sequences are then written to separate files based on the prediction type (AF2 or ESMfold) and saved in the specified output folder.

		If percentage is provided, a specified percentage of positions in each sequence is randomly selected and mutated, while the remaining positions are kept unchanged. 
		The native sequence is kept unchanged in the output files for further comparison during structure predictions.
		Mutated sequences are written to separate files based on the prediction type and the specified mutation percentage.
		'''
		designed_res_indices = mpnn_analysis.find_redes_residues(self) #list of redesigned positions
		picked_letters = []
		list_scrambles = []
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		all_seqs = [lines[i].strip('\n') for i in range(len(lines)) if i %2 != 0]
		orig_seq = all_seqs[0] #the first sequence will be kept unchanged since there are no differing residues
		if percentage == None:
			for i in range(len(all_seqs)):
				scramble_design = '' 
				for index, residue in enumerate(orig_seq):
					if index in designed_res_indices[i]:
						guess_letter = random.choice(letters)
						picked_letters.append(guess_letter)
						scramble_design += guess_letter
					else:
						scramble_design += residue
				list_scrambles.append(scramble_design)
			if self.str_predict == 'AF2':
				for seq_id in range(len(list_scrambles)):
					with open(self.folder_out+'AF2_'+'/'+'{}_{}.a3m'.format(seq_id, self.filename), 'w') as single_seq:
						single_seq.write('>{}={}\n'.format(self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
			if self.str_predict == 'ESMfold':
				with open(self.folder_out+'ESMfold_'+'/'+'{}.fa'.format(self.filename), 'w') as seqs_fa:
					for seq_id in range(len(list_scrambles)):
						seqs_fa.write('>{}={}\n'.format(self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
			return list_scrambles
		else: 
			number_pos_mutate = int(round(percentage*len(orig_seq),0))
			for i in range(len(all_seqs)):
				scramble_design = '' 
				indices_mutate = sorted(random.sample(range(len(orig_seq)),number_pos_mutate)) #selecting random positions to mutate
				if i == 0: #native sequence
					for j in range(len(orig_seq)):
						scramble_design += orig_seq[j] #keeping native sequence untouched for further comparison in structure predictions
					list_scrambles.append(scramble_design)
				else:
					for j in range(len(orig_seq)):
						if j in indices_mutate:
							guess_letter = random.choice(letters)
							scramble_design += guess_letter
						else:
							scramble_design += orig_seq[j]
					list_scrambles.append(scramble_design)
			if self.str_predict == 'AF2':
				with open(self.folder_out + '/' + '{}_scr_AF2_{}.fa'.format(percentage, self.filename), 'w') as seqs_fa:					
					for seq_id in range(len(list_scrambles)):
						if seq_id == 0:
							seqs_fa.write('>native_{}={}\n'.format(self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
						else:
							seqs_fa.write('>percent_scr_{}_{}={}\n'.format(percentage, self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
			if self.str_predict == 'ESMfold':
				with open(self.folder_out + '/' + '{}_scr_ESMfold_{}.fa'.format(percentage, self.filename), 'w') as seqs_fa:
					for seq_id in range(len(list_scrambles)):
						if seq_id == 0:
							seqs_fa.write('>native_{}={}\n'.format(self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
						else:
							seqs_fa.write('>percent_scr_{}_{}={}\n'.format(percentage, self.filename, seq_id)+str(list_scrambles[seq_id])+'\n')
			return list_scrambles
		
	def masking_pdb(self, scramble='G'):
		'''
		Takes as input clean pdb containing only atoms and renames residues to either GLY or to any random residue; used for checking aa identity influence on inference of pMPNN
		is perhaps not bug-free, has to be benchmarked
		'''
		_3letters = [i for i in aas_alphabet.keys()]
		masked_pdb = []
		resids_lst = []    
		if scramble != 'G':
			with open(self.pdb,'r') as pdb_input:
				lines = pdb_input.readlines()
				random_aa = random.choices(_3letters)[0] # initializing the first random aa choice 
				aa_keep = random_aa
			for line in range(len(lines)):
				if line > 0: # index 0 == REMARK, should be ignored 
					res_numb = lines[line][22:26] # residue number 
					resids_lst.append(res_numb.split()[0]) # append of the residue number to keep counting whether it changes or not in the next iter step
					if len(resids_lst) > 0:
						if resids_lst[line-1] == resids_lst[line-2]: # if residue doesnt change 
							masked_pdb.append(lines[line][:17]+str(aa_keep)+lines[line][20:]) # keep the same aa as it was before
						else:
							random_aa = random.choices(_3letters)[0]
							aa_keep = random_aa
							masked_pdb.append(lines[line][:17]+str(aa_keep)+lines[line][20:])
			with open(self.folder_out+'/'+'masked_pdb.pdb', 'w') as masked_pdb_file:
				masked_pdb_file.write(''.join(masked_pdb))
			return masked_pdb
		else: 
			with open(self.pdb,'r') as pdb_input:
				lines = pdb_input.readlines()
				aa_scr = 'GLY' 
				print(lines[0])# initializing the first random aa choice 
			for line in range(len(lines)):
				if line != 'REMARK':
					masked_pdb.append(lines[line][:17]+str(aa_scr)+lines[line][20:])
				if line > 0: # index 0 == REMARK, should be ignored 
					res_numb = lines[line][22:26] # residue number 
					resids_lst.append(res_numb.split()[0]) # append of the residue number to keep counting whether it changes or not in the next iter step
					if len(resids_lst) > 0:
						if resids_lst[line-1] == resids_lst[line-2]: # if residue doesnt change 
							masked_pdb.append(lines[line][:17]+str(aa_scr)+lines[line][20:]) # keep the same aa as it was before
						else:
							masked_pdb.append(lines[line][:17]+str(aa_scr)+lines[line][20:])
			with open(self.folder_out+'/'+'maskedgly_pdb.pdb', 'w') as masked_pdb_file:
				masked_pdb_file.write(''.join(masked_pdb))