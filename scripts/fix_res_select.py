import numpy as np 
import MDAnalysis
from mpnn_analysis import *

class res_selection:

	def __init__(self, start_structure, cofactor):
		'''
		start_structure (list): list of strings encoding hetatm entries
		'''
		self.start_structure = start_structure
		self.cofactor = cofactor

	def find_hetatoms(self):
		'''
		Finds and returns a list of unique (non-water) heteroatom names present in the PDB file.

		Parameters:
				- self: The object itself.

		Returns:
				- list: A list of unique heteroatom names found in the PDB file.
		'''
		with open (self.start_structure) as pdb:
			pdb_lines = pdb.readlines()
			hetatm_name = []
			for line in pdb_lines:
				if line.startswith('HETATM'):
					if line[17:20] != 'HOH' and line[17:20] not in hetatm_name:
						hetatm_name.append(line[17:20].strip(' '))
		if len(hetatm_name) != 0:
			print('--- heteroatoms found:',hetatm_name, ' ---')
		return hetatm_name
				
	def fix_residues(self, angstroem):
		'''
		Identifies residues in close proximity to the specified HETATM (cofactor or a substrate) 
		and returns their indices (0-indexed).

		Parameters:
				- angstroem (int): The distance in angstroms to consider for identifying nearby residues.

		Returns:
				- list: A list of indices of residues in close proximity to the specified cofactor.
		'''
		hetatoms_list = res_selection.find_hetatoms(self)
		if len(self.cofactor) == 1:
			hetatm = self.cofactor[0]
			if hetatm in hetatoms_list:
				index_cofactor = hetatoms_list.index(hetatm)
			universe_start_struct = MDAnalysis.Universe(self.start_structure)
			selection_fix = f'protein and byres around {angstroem} resname {str(hetatoms_list[index_cofactor])}'
			fix_res = universe_start_struct.select_atoms(selection_fix)
			fix_res_data = mpnn_analysis.data_universe(fix_res)
			print('--- selected fixed residues for pMPNN design around {} ---'.format(hetatm)) 
			indices_fix = [i for i in fix_res_data['res_numb']]
			return list(dict.fromkeys(indices_fix))
		else:
			indices_fix = []
			for i in range(len(hetatoms_list)):
				if hetatoms_list[i] in self.cofactor:
					hetatm = hetatoms_list[i]
					index_cofactor = hetatoms_list.index(hetatm)
					universe_start_struct = MDAnalysis.Universe(self.start_structure)
					selection_fix = f'protein and byres around {angstroem} resname {str(hetatoms_list[index_cofactor])}'
					fix_res = universe_start_struct.select_atoms(selection_fix)
					fix_res_data = mpnn_analysis.data_universe(fix_res)
					indices_iter = [i for i in fix_res_data['res_numb']]
				print('--- selected fixed residues for pMPNN design around {} ---'.format(hetatm)) 
				for i in indices_iter:
					indices_fix.append(i)
			return list(dict.fromkeys(indices_fix))