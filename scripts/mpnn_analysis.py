import sys
import pandas as pd
import numpy as np
import re
import json
import os
import freesasa as fa
import MDAnalysis
import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import time
import MDAnalysis.analysis.rms
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import seq1
from Bio.Align import AlignInfo
import biotite.structure.io as bsio

# Definition of the global variables 

letters = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
					"M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

_1to3_ = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
					'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU',
					'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 
					'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

polarity_aas = {'hydrophilic':['ARG', 'LYS', 'SER', 'THR', 'GLU', 'ASP', 'GLN', 'ASN', 'HIS', 'TYR'],
					 'hydrophobic':['GLY', 'PRO', 'PHE', 'ALA', 'VAL', 'LEU', 'ILE', 'TRP', 'MET', 'CYS']}

aas_alphabet = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
						 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
						 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
						 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

charged_aas = {'positive':['ARG', 'LYS'],
	 'negative':['ASP', 'GLU']}

class mpnn_analysis:

	def __init__(self, fasta, json_folder, pdb_folder, str_predict=None, figs_folder=None):
		self.fasta = fasta
		self.json_folder = json_folder
		self.pdb_folder = pdb_folder
		self.str_predict = str_predict
		self.figs_folder = figs_folder
	
	def seq_identity(self):
		"""
		Calculates the sequence identity scores between the original sequence and redesigned sequences from the fasta output of pMPNN.

		Returns a list of sequence identity scores, representing the proportion of amino acids in each 
		redesigned sequence that match the corresponding amino acids in the original sequence.

		Parameters:
				self.fasta (str): The path to the input FASTA file containing multiple sequences.

		Returns:
				list: A list of sequence identity scores, rounded to four decimal places, for each redesigned 
							sequence in comparison to the original sequence.
		"""
		with open(self.fasta) as fasta_input:
				lines = fasta_input.readlines()
		matched_aas = ''
		all_matched_aas = []
		all_seqs = [lines[i].strip('\n') for i in range(len(lines)) if i %2 != 0]
		orig_seq = all_seqs[0] #assuming that the first sequence is the input sequence; index[1] because [0] in .fa is the header of the sequence 
		seqs_to_match = all_seqs[1:] #all redesigned sequences
		for seq_match in seqs_to_match:
			for index in range(len(orig_seq)):
					if orig_seq[index] == seq_match[index]:
							matched_aas += seq_match[index]
			all_matched_aas.append(matched_aas)
			matched_aas = '' 
		scores_identity = [round(len(i)/len(orig_seq),4) for i in all_matched_aas]
		return scores_identity
	
	def seq_recovery_extract(self):
		"""
		Extracts information about sequence recovery and sample names from the fasta output from pMPNN.

		Parameters:
		self.fasta (str): The path to the input FASTA file containing multiple sequences.

		Returns two lists:
		'sample_name': A list containing sample names extracted from the 'sample=' pattern in the header 
					 lines.
		'recovery': A list containing sequence recovery extracted from the 
					'seq_recovery=' pattern in the header lines.
		"""
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		recovery = []
		sample_name = []
		recovery_pattern = re.compile(r'seq_recovery=(.*)')
		sample_pattern = re.compile(r'sample=\d*')
		for index in range(len(lines[2:])): #0, 1 will be the header and the orig seq
			if index % 2 == 0: #every second line is the header of the sequence
				recovery.append(recovery_pattern.findall(lines[2:][index])[0]) 
				sample_name.append(sample_pattern.findall(lines[2:][index])[0])
		return sample_name, recovery
	
	def find_redes_residues(self):
		"""
		Finds indices at which each the original sequence identity was mutated

		Parameters:
			self.fasta (str): The path to the input FASTA file containing multiple sequences.

		Returns:
			list: A list of lists, where each inner list contains indices of differing residues for 
					corresponding redesigned sequences.
		"""
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		all_seqs = [lines[line].strip('\n') for line in range(len(lines)) if line %2 != 0]
		all_seq_redes = []
		orig_input = all_seqs[0] 
		for seq in all_seqs:
			index_redes = [] #each iteration will reset the list of differing indices
			for index in range(len(seq)):
				if seq[index] != orig_input[index]:
					index_redes.append(index)
			all_seq_redes.append(index_redes)
		return all_seq_redes #0-index is empty because seqs[0] is input itself -> identical sequence, no differences
	
	def get_pdbs(self):
		"""
		Retrieves a list of numerically sorted PDB files from the specified folder based on the protein structure prediction network.

		Parameters: 
			self.fasta mentioned above
			self.str_predict (str) = 'AF2', 'ESMfold' - structure prediction method

		Returns:
				list: A sorted list of PDB files based on the prediction type.
							If the prediction type is 'RF2', an empty list is returned.

		Example:
				If 'self.str_predict' is 'AF2' and the folder contains the following files:
				- .model_001_rank_001.pdb
				- model_001_rank_001.pdb
				- model_002_rank_002.pdb
				- .DS_Store
				The function will return ['model_001_rank_001.pdb', 'model_002_rank_002.pdb'].
		"""
		def sort_models(list_models):
			only_int = lambda text: int(text) if text.isdigit() else text
			alphanum = lambda key: [only_int(i) for i in re.split('([-+]?[0-9]*\.?[0-9]*)', key)]
			list_models.sort(key=alphanum)
			return list_models
		if self.str_predict == 'AF2': 
			files = os.listdir(self.pdb_folder)
			new_files = [name for name in files if name[0]!='.']
			toprank_pattern = re.compile(r'_rank_00(.)')
			models_full = [model for model in new_files if '_000' in model]
			toprank_models = [toprank for toprank in models_full if toprank_pattern.findall(str(toprank))[0] == '1']
			if len(toprank_models) != 0:
				sorted_toprank_pdb = sort_models(toprank_models)
				return sorted_toprank_pdb
			else: 
				sorted_toprank_pdb = sort_models(models_full)
				return sorted_toprank_pdb
		if self.str_predict == 'ESMfold':
			files = os.listdir(self.pdb_folder)
			new_files = [name for name in files if name[0]!='.']
			all_models_ESMfold = [model for model in new_files]
			sorted_models_ESMfold = sort_models(all_models_ESMfold)
			return sorted_models_ESMfold
		
	def get_plddts(self):
		"""
		Retrieves plddt scores from JSON files or calculates them from PDB structures (bfactors).

		Depending on the value of 'self.str_predict', this function retrieves plddt scores from JSON files
		or calculates them from PDB structures. The function performs different actions based on different
		prediction types:
		- For 'AF2' prediction type:
			- Retrieves plddt scores from JSON files in the folder specified by 'self.json_folder'.
			- Identifies top-ranked models based on filenames containing '_rank_001'.
			- Returns the plddt score of the input model and plddt scores of top-ranked models.
		- For 'ESMfold' prediction type:
			- Retrieves PDB files using the 'get_pdbs' function.
			- Calculates plddt scores from PDB structures (bfactors) and returns them.

		Returns:
				tuple: A tuple containing:
							 - The plddt score of the input model.
							 - A list of plddt scores of other models.
		"""
		def sort_models(list_models):
			only_int = lambda text: int(text) if text.isdigit() else text
			alphanum = lambda key: [only_int(i) for i in re.split('([-+]?[0-9]*\.?[0-9]*)', key)]
			list_models.sort(key=alphanum)
			return list_models
		plddts = []
		if self.str_predict == 'AF2':
			files = os.listdir(self.json_folder)
			new_files = [name for name in files if name[0]!='.']
			toprank_pattern = re.compile(r'_rank_00(.)') #will always have the following naming
			models_full = [model for model in new_files if '_000' in model]
			toprank_models = [toprank for toprank in models_full if toprank_pattern.findall(str(toprank))[0] == '1']
			if len(toprank_models) != 0:
				sorted_toprank_json = sort_models(toprank_models)
			else: 
				sorted_toprank_json = sort_models(models_full)
			for index in sorted_toprank_json:
				with open(self.json_folder+'/'+str(index), 'r') as model:
					ranking_dict = json.load(model)
					plddts.append(round(np.mean(ranking_dict['plddt']),1))
			print('input_plddt:', plddts[0]) #printing the plddt of the input
			return plddts[0], plddts[1:] #return only the models, not the input  
		if self.str_predict == 'ESMfold':
			files = mpnn_analysis.get_pdbs(self)
			new_files = [name for name in files if name[0]!='.']
			for model in new_files: 
				struct = bsio.load_structure(self.pdb_folder+'/'+model, extra_fields=["b_factor"])
				plddt = round(struct.b_factor.mean(),2)
				plddts.append(plddt)
			print('input_plddt:', plddts[0])
			return plddts[0], plddts[1:]
		
	def get_sasa(atomgroup, selection='R'):
		"""
		Calculates solvent accessible surface area (SASA) for atoms or residues given a universe object from MDanalysis.

		This function calculates the solvent accessible surface area (SASA) for atoms or residues in an 
		MDAnalysis atomgroup. It uses the Shrake-Rupley algorithm with a probe radius of 1.4 Å. The SASA 
		is calculated based on the atomic coordinates and radii, and it can be computed for individual atoms 
		('A'), residues ('R'), or the entire structure ('S').

		Parameters:
				atomgroup (MDAnalysis AtomGroup): The MDAnalysis AtomGroup for which SASA is to be calculated.
				selection (str, optional): The selection for which SASA is calculated. 
																		Options: 'A' (atoms), 'R' (residues), 'S' (entire structure).
																		Default is 'R'.
		Returns:
				numpy.ndarray: An array containing the SASA values calculated based on the specified selection.
		"""
		#set algorithm to shrake-rupley, set probe radius if neccesary
		param = fa.Parameters()
		param.setAlgorithm('ShrakeRupley')
		param.setProbeRadius(1.4)
		# standard radii taken from the biopython implementation of shake-rupley
		radii_atomtypes= {"H": 1.200, "HE": 1.400, "C": 1.700, "N": 1.550, "O": 1.520, "F": 1.470,
									"NA": 2.270, "MG": 1.730, "P": 1.800, "S": 1.800, "CL": 1.750, "K": 2.750,
									"CA": 2.310, "NI": 1.630, "CU": 1.400, "ZN": 1.390, "SE": 1.900, "BR": 1.850,
									"CD": 1.580, "I": 1.980, "HG": 1.550}
		# calcuating sasa requires coordinates and radii per atom
		atomtypes = atomgroup.atoms.types
		radii = np.vectorize(radii_atomtypes.get)(atomtypes)
		coords = atomgroup.atoms.positions.flatten()
		sasa_result = fa.calcCoord(coords, radii, param)
		sasa_atom = np.array([sasa_result.atomArea(i) for i, a in enumerate(atomgroup.atoms)])
		# sasa per atom
		if selection == 'A':
			return sasa_atom
		# sasa per residue
		if selection == 'R':
			# if array is not equal to array shifted by one, then a transition took place
			resids = atomgroup.atoms.resids
			resids_shifted = np.concatenate((resids[-1:], resids[0:-1]))
			transitions = np.concatenate((np.where(resids != resids_shifted)[0], [len(resids)]))
			# use transitions to get start/end of atoms with same resid to split up by residues
			return np.array([np.sum(sasa_atom[start:end]) for start, end in zip(transitions[0:-1], transitions[1:])])
		# sasa for entire structure
		if selection == 'S':
			return np.sum(sasa_atom)
	
	def data_universe(universe):
		"""
		Processes MDAnalysis Universe object to extract per-residue data.

		This function processes an MDAnalysis Universe object to extract residue data, including residue 
		name, residue number, chain name, solvent accessible surface area (SASA), and hydrophobicity.
		The SASA is calculated using the 'get_sasa' function from the 'mpnn_analysis' module. The 
		hydrophobicity is determined based on the polarity of amino acids.

		Args:
				universe (MDAnalysis Universe): MDAnalysis Universe object containing molecular data.

		Returns:
				pandas.DataFrame: A DataFrame containing the extracted residue data, including residue name, 
													residue number, chain name, SASA, hydrophobicity, and fraction of SASA 
													compared to total SASA for each residue for a given universe obj

		polarity of aas classified based on https://doi.org/10.1021/bi00405a042
		(kcal/mol) of transfer from dilute solution in cyclohexane to dilute aqueous solution at pH 7.0
		used in: https://doi.org/10.1186/1471-2105-9-272

		free-in-solution SASA of amino acids taken from: 10.1007/s00894-009-0454-9
		"""
	
		total_sasa_res = {'ALA':209.02, 'ARG':335.73, 'ASN':259.85, 'ASP':257.99, 'CYS':240.5, 'GLN':286.76, 'GLU':285.03, 
									'GLY':185.15, 'HIS':290.04, 'ILE':273.46, 'LEU':278.44, 'LYS':303.43, 'MET':291.52, 'PHE':311.30, 
									'PRO':235.41, 'SER':223.04, 'THR':243.55, 'TRP':350.68, 'TYR':328.82, 'VAL':250.09}

		res_seq_pattern = re.compile(r'<Residue\s(\w{3}),') #residue name pattern as 3-letter code
		res_numb_pattern = re.compile(r',\s{1}\d+') #residue number pattern
		res_all = [i for i in universe.residues if str(res_seq_pattern.findall(str(i))[0]) in total_sasa_res.keys()]
		res_numb = [j.split()[1] for j in res_numb_pattern.findall(str(res_all))]
		res_name = res_seq_pattern.findall(str(res_all))
		hydrophobicity = []
		
		chain_id = [ch_id for ch_id in universe.residues.segids]
		sasa_values = mpnn_analysis.get_sasa(universe, selection='R') #calls a function to calculate per-residue SASA 

		for index, res in enumerate(res_name):
			if res in list(polarity_aas.values())[0]: #if in hydrophilic values -> append hydrophilic 
				hydrophobicity.append('hydrophilic')
			else:
				hydrophobicity.append('hydrophobic')

		dict_universe = {'resname':res_name, 
										'res_numb':res_numb,
										'chain_name':chain_id,
										'sasa_value':sasa_values, 
										'hydrophobicity':hydrophobicity}
		
		universe_df = pd.DataFrame(data=dict_universe)
		universe_df['res_numb'] = universe_df['res_numb'].astype(int)
		fraction_sasa = {}
		for index_df in range(len(universe_df['sasa_value'])):
			if index_df in universe_df.index:
				dict = {universe_df['resname'][index_df]:universe_df['sasa_value'][index_df]} #creating amino_acid:SASA_value pairs in dict from the pd df
				for amino_acid in list(total_sasa_res.keys()):
					if amino_acid == list(dict.keys())[0]:
						frac = dict[amino_acid] / total_sasa_res[amino_acid]
						fraction_sasa[index_df] = frac
		universe_df['fraction_sasa'] = fraction_sasa.values()
		return universe_df

	def get_data(self): 
		"""
		Retrieves data from MDAnalysis Universe objects for multiple models.

		This function retrieves data from MDAnalysis Universe objects for multiple models obtained from the 
		'get_pdbs' function. For each model, it creates an MDAnalysis Universe object from the corresponding 
		PDB file and extracts residue data using the 'data_universe' function. The extracted data includes 
		residue name, residue number, chain name, SASA, and hydrophobicity, as calculated by the 'data_universe' 
		function. 

		Returns:
				list: A list containing extracted residue data for each model.
		"""
		list_pdb_files = mpnn_analysis.get_pdbs(self)
		data_models = []
		for model in list_pdb_files:
			model_universe = MDAnalysis.Universe(self.pdb_folder+'/'+str(model), dt=1)
			data_models.append(mpnn_analysis.data_universe(model_universe))
		return data_models
	
	def designed_res_exposure_hydrophilicity(self):
		"""
		Calculates the fraction of solvent-exposed redesigned residues and the fraction of hydrophilic residues among 
		the exposed redesigned residues.

		Returns:
				tuple: A tuple containing:
							 - List of fractions of exposed residues in the reference model.
							 - List of fractions of exposed residues in the redesigned models at the corresponding redesigned indices from reference.
							 - List of fractions of hydrophilic exposed residues among the designed ones in the reference model.
							 - List of fractions of hydrophilic exposed residues among the designed ones in the redesigned models.

		This function calculates the fraction of solvent-exposed redesigned residues and the fraction of hydrophilic 
		residues among the exposed redesigned residues. It first retrieves data for all models using the get_data() 
		function. Then, for each model, it calculates the fraction of exposed redesigned residues and the fraction 
		of hydrophilic exposed residues among the designed ones. Finally, it returns these fractions.

		"""
		start = time.time()
		print('--- reading pdbs ---')
		list_data_models = mpnn_analysis.get_data(self)
		designed_res_indices = mpnn_analysis.find_redes_residues(self)
		frac_exposed_ref = []
		frac_exposed_redes = []
		frac_hydrophilic_subset_ref = []
		frac_hydrophilic_subset_redes = []
		reference_universe = list_data_models[0]
		to_compare_universes = list_data_models[1:]
		print('--- calculating fraction of exposed redesigned residues ---')
		for index in range(len(to_compare_universes)):
			#iteratively based on the index of designed residues subsets the reference AF2 prediction and calculate fraction of hydrophilic residues 
			exposed_condition_ref = reference_universe.iloc[designed_res_indices[1:][index]]['fraction_sasa'] > 0.2 #[1:] bc [0] is empty -> wt sequence
			exposed_condition_des = to_compare_universes[index].iloc[designed_res_indices[1:][index]]['fraction_sasa'] > 0.2
			exposed_fraction_ref = round(len(reference_universe.iloc[designed_res_indices[1:][index]][exposed_condition_ref])/len(reference_universe.iloc[designed_res_indices[1:][index]]),3)
			exposed_fraction_design = round(len(to_compare_universes[index].iloc[designed_res_indices[1:][index]][exposed_condition_des])/len(to_compare_universes[index].iloc[designed_res_indices[1:][index]]),3)    
			frac_exposed_ref.append(exposed_fraction_ref)
			frac_exposed_redes.append(exposed_fraction_design)
			#same here but then calculating frac of hydrophilic exposed residues
			hydrophil_subset_ref = reference_universe.iloc[designed_res_indices[1:][index]][exposed_condition_ref].groupby('hydrophobicity')['chain_name'].count() #grouping by chain -> benchmarked for monomers
			hydrophil_subset_redes = to_compare_universes[index].iloc[designed_res_indices[1:][index]][exposed_condition_des].groupby('hydrophobicity')['chain_name'].count()
			frac_hydrophilic_subset_ref.append(round(hydrophil_subset_ref[0]/(hydrophil_subset_ref[0]+hydrophil_subset_ref[1]),3))
			frac_hydrophilic_subset_redes.append(round(hydrophil_subset_redes[0]/(hydrophil_subset_redes[0]+hydrophil_subset_redes[1]),3))

		end = time.time()
		print('--- finished after', round(end-start, 1), 'seconds ---')
		return frac_exposed_ref, frac_exposed_redes, frac_hydrophilic_subset_ref, frac_hydrophilic_subset_redes

	def find_exposed_reference(self):
		"""
		Identifies solvent-exposed residues in the reference model.

		Returns:
				list: A list of indices corresponding to solvent-exposed residues in the reference model.

		"""
		list_pdb_files= mpnn_analysis.get_pdbs(self)
		reference_model = MDAnalysis.Universe(self.pdb_folder+'/'+str(list_pdb_files[0]), dt=1)
		reference_universe = mpnn_analysis.universe_sasa(reference_model)
		exposed_condition = reference_universe[reference_universe['fraction_sasa'] > 0.2]
		return list(exposed_condition.index)
	
	def find_core_size(self):
		"""
		This function calculates the size of the protein core by counting the number of hydrophobic residues with a 
		solvent-accessible surface area (SASA) value of 0.
		
		Returns:
				tuple: A tuple containing:
							- pandas.DataFrame:: A DataFrame containing the count of hydrophobic residues in the core and annotations 
													('wt' for the wild-type structure and 'des' for the designed structures).
							- list: A list containing the count of hydrophobic residues in the core for designed structures.

		"""
		start=time.time()
		print('--- calculating core size ---')
		list_data_models = mpnn_analysis.get_data(self)
		count_hb = []
		for structure in list_data_models:
			count_hb.append(len(structure[structure['sasa_value']==0])) #relative definition of the "core" is SASA == 0; benchmarked only for medium-sized proteins
		annot_wt = ['wt' for i in range(len(count_hb)) if i == 0]
		annot_des = ['des' for i in range(len(count_hb)) if i != 0]
		annot = annot_wt + annot_des
		dict_hb_core = {'count_hb':count_hb,
										'annot':annot}
		df_hb_core = pd.DataFrame(dict_hb_core, columns=dict_hb_core.keys())
		end=time.time()
		print('--- finished after', round(end-start,1),'---')
		return df_hb_core, count_hb[1:] #count_hb[1:] = only for designed variants; [0] is wt
	
	def rmsd_calculate(self, struct_compare=None):
		"""
		Calculates the root-mean-square deviation (RMSD) of the designed structures compared to the reference structure.

		Args:
				struct_compare (str, optional): Path to the reference structure. If None, the first structure in the list 
																				of designed structures will be used as the reference. Defaults to None.

		Returns:
				list: A list containing the RMSD values of the designed structures compared to the reference structure.

		This function calculates the RMSD of the backbone atoms of the designed structures compared to the reference 
		structure. It first retrieves the list of PDB files using the get_pdbs() function. Then, it initializes MDAnalysis 
		universes for each PDB file. If the struct_compare argument is provided, it sets the reference structure to the 
		structure specified. Otherwise, it uses the first structure in the list of designed structures as the reference. 
		For each designed structure, it calculates the RMSD of the backbone atoms relative to the reference structure. 
		The RMSD values are stored in a list and returned.
		"""
		rmsd_list = [] 
		list_pdb_files = mpnn_analysis.get_pdbs(self)
		universes = []
		for model in list_pdb_files: 
			universe = MDAnalysis.Universe(self.pdb_folder+'/'+str(model), dt=1)
			universes.append(universe)
		print('--- calculating backbone rmsd of the designed structures to the reference ---')
		start = time.time()
		if struct_compare == None:
			reference_str = universes[0]
			to_compare_str = universes[1:]
			for index in range(len(to_compare_str)):
				R = MDAnalysis.analysis.rms.RMSD(to_compare_str[index], reference_str,
					select="backbone")
				R.run()
				rmsd = R.results.rmsd.T
				rmsd_list.append(rmsd[2][0]) #retrieving rmsd value from the transposed rmsd vector 
			end = time.time()
			print('--- finished after', round(end-start, 2), 'seconds ---')
			return rmsd_list
		else:
			reference_str = MDAnalysis.Universe(struct_compare)
			ref_str_noH = reference_str.select_atoms("protein and not (name H*)")
			to_compare_str = universes[1:]
			for index in range(len(to_compare_str)):
				R = MDAnalysis.analysis.rms.RMSD(to_compare_str[index], ref_str_noH,
					select="backbone")
				R.run()
				rmsd = R.rmsd.T
				rmsd_list.append(rmsd[2][0])
			end = time.time()
			print('--- finished after', round(end-start, 2), 'seconds ---')
			return rmsd_list

	def charge_analysis(self, selection='total'):
		"""
		Analyzes the charge of the pMPNN-designed variants.

		Args:
			selection (str, optional): Specifies the type of charge analysis to perform. 
			'total' to calculate the total charge of the structure
			'exposed' to calculate the charge of exposed residues only.

		Returns:
				tuple: A tuple containing charge information depending on the selection:
							 - For 'total' selection: (charge_fullstr_ref, charge_fullstr_des)
							 - For 'exposed' selection: (exposed_charge_ref, exposed_charge_des)
		"""
		designed_res_indices = mpnn_analysis.find_redes_residues(self)
		start = time.time()
		print('--- calculating total charge ---')
		list_data_models = mpnn_analysis.get_data(self)
		reference_universe = list_data_models[0]
		to_compare_universes = list_data_models[1:]
		if selection == 'exposed':
			exposed_charge_ref = [] 
			exposed_charge_des = []
			for index in range(len(to_compare_universes)):
				exposed_condition_ref = reference_universe.iloc[designed_res_indices[1:][index]]['fraction_sasa'] > 0.2
				exposed_condition_redes = to_compare_universes[i].iloc[designed_res_indices[1:][index]]['fraction_sasa'] > 0.2
				charge_ref_subset = reference_universe.iloc[designed_res_indices[1:][index]][exposed_condition_ref]
				charge_redes_subset = to_compare_universes[i].iloc[designed_res_indices[1:][index]][exposed_condition_redes]
				charge_ref = []
				charge_des = []
				for res_name in charge_ref_subset['resname']:
					if res_name in list(charged_aas.values())[0]:
						charge_ref.append(1)
					if res_name in list(charged_aas.values())[1]:
						charge_ref.append(-1)
				exposed_charge_ref.append(sum(charge_ref))
				for res_name_des in charge_redes_subset['resname']:
					if res_name_des in list(charged_aas.values())[0]:
						charge_des.append(1)
					if res_name_des in list(charged_aas.values())[1]:
						charge_des.append(-1)
				exposed_charge_des.append(sum(charge_des))
			return exposed_charge_ref, exposed_charge_des
		if selection=='total': 
			charge_fullstr_ref = 0
			charge_indiv_ref = []
			for res_name in reference_universe['resname']:
				if res_name in list(charged_aas.values())[0]:
					charge_indiv_ref.append(1)
				if res_name in list(charged_aas.values())[1]:
					charge_indiv_ref.append(-1)
			charge_fullstr_ref=sum(charge_indiv_ref)        
			charge_fullstr_des = []

			for i in range(len(to_compare_universes)):
				des = to_compare_universes[i]
				charge_indiv_des = []
				for res_name in des['resname']:
					if res_name in list(charged_aas.values())[0]:
						charge_indiv_des.append(1)
					if res_name in list(charged_aas.values())[1]:
						charge_indiv_des.append(-1)
				charge_fullstr_des.append(sum(charge_indiv_des))
			end = time.time()
			print('--- finished after',round(end-start,1),'---')
			return charge_fullstr_ref, charge_fullstr_des
			
	def gen_pssm(self, heatmap='no'):
		"""
		Generates a Position-Specific Scoring Matrix (PSSM) for the input sequence and all designed sequences (R = Lx20; L=len(seq))

		Args:
			self: An instance of the class containing the gen_pssm function.
			heatmap (str, optional): Specifies whether to return PSSM matrices for subsequent heatmap generation. 

		Returns:
			tuple: A tuple containing PSSM matrices based on the specified heatmap option:
							 - If heatmap is 'yes': (orig_pssm, norm_pssm_msa_df)
							 - If heatmap is 'no': (orig_pssm_copy, norm_pssm_msa_df)

		Note:
			This function assumes that the 'letters' global variable is defined, representing the 20 standard amino acids.
		"""
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		all_seqs = [lines[i].strip('\n') for i in range(len(lines)) if i %2 != 0]
		orig_seq = all_seqs[0]
		records = [SeqRecord(Seq(seq), id=f"seq_{i+1}", description=f"Sequence {i+1}") for i, seq in enumerate(all_seqs[1:])]
		alignment = MultipleSeqAlignment(records)
		summary_align = AlignInfo.SummaryInfo(alignment)
		pssm = summary_align.pos_specific_score_matrix()
		total_values = []
		for pos in range(len(orig_seq)):
			values_list = []
			for amino_acid in letters:
				values_list.append(pssm[pos][amino_acid])
			total_values.append(values_list)
		pssm_msa_df = pd.DataFrame(total_values, columns=[i for i in letters])
		norm_pssm_msa_df = (pssm_msa_df-pssm_msa_df.min())/(pssm_msa_df.max()-pssm_msa_df.min())
		
		#part for generation of pssm of only original sequence
		mat_orig =np.zeros((len(orig_seq), 20))
		for pos in range(len(orig_seq)):
			letter = orig_seq[pos]
			substitute_pos = letters.index(letter) #global var letters
			mat_orig[pos][substitute_pos] += 1
		orig_pssm = pd.DataFrame(mat_orig, columns = [i for i in letters])
		orig_pssm_copy = orig_pssm.copy()
		orig_pssm.replace(0, np.nan, inplace=True)
		if heatmap == 'yes':
			return orig_pssm, norm_pssm_msa_df
		else:
			return orig_pssm_copy, norm_pssm_msa_df

	def aggregate_data(self):
		"""
		This function aggregates different data metrics obtained from various analysis functions into a single DataFrame 
		for comprehensive analysis. The aggregated metrics include sample name, recovery, sequence identity, pLDDT values, 
		fraction of exposed residues (both reference and designed), hydrophilicity of exposed residues (both reference and designed), 
		total charges (both reference and designed), and core size. (R = L x 10 where 10 = columns corresponding to the mentioned metrics)

		Returns:
				pandas.DataFrame: DataFrame containing aggregated data metrics including sample name, recovery, sequence identity,
									predicted pLDDT values, fraction of exposed residues (reference and designed), hydrophilicity of exposed residues,
									total charges, and core size.

		"""
		start = time.time()
		sample_name, recovery = mpnn_analysis.seq_recovery_extract(self)
		seq_identity = mpnn_analysis.seq_identity(self)
		plddt_input, plddt_des = mpnn_analysis.get_plddts(self)
		frac_exposed_ref, frac_exposed_des, hydrophil_exposed_ref, hydrophil_exposed_des = mpnn_analysis.designed_res_exposure_hydrophilicity(self)
		charges_total_ref, charges_total_des = mpnn_analysis.charge_analysis(self, selection='total')
		df_core, core_sizes = mpnn_analysis.find_core_size(self)
		dict_mpnn = {'sample_name':sample_name,
						'seq_recovery':recovery,
						'seq_identity':seq_identity,
						'plddt':plddt_des,
						'frac_exposed_ref':frac_exposed_ref,
						'frac_exposed_redes':frac_exposed_des, 
						'hydrophilic_exposed_ref':hydrophil_exposed_ref, 
						'hydrophilic_exposed_redes':hydrophil_exposed_des,
						'charges_total':charges_total_des,
						'core_size':core_sizes}
		out_df = pd.DataFrame(dict_mpnn)
		out_df['seq_recovery'] = out_df['seq_recovery'].astype(float)
		end = time.time()
		print('--- total running time of the script:', round(end-start,1),'seconds ---')
		return out_df
