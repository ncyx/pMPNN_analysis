import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr
from statannotations.Annotator import Annotator
from mpnn_analysis import *

class plot_mpnn_analysis:

	def __init__(self, dataframe_analysis, fasta, json_folder, pdb_folder, figs_out, str_predict=None):
		self.dataframe = dataframe_analysis
		self.fasta = fasta
		self.json_folder = json_folder
		self.pdb_folder = pdb_folder
		self.figs_out = figs_out
		self.str_predict = str_predict

	def plot_seq_recovery(self, save='yes'):
		'''
		Function to plot sequence recovery vs sequence identity.

		Args:
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				None
		'''
		seq_identity = mpnn_analysis.seq_identity(self)
		name_fig = 'seq_recovery_designs'
		fig = plt.figure(figsize=(6.4 ,4.8))
		ax1 = fig.add_subplot(111)
		ax1.hist(self.dataframe['seq_recovery'], alpha=0, color='#A7D5B9',label='native_recovery', edgecolor='white', 
						linewidth=1.0)
		ax2=ax1.twiny()
		ax2.xaxis.set_ticks_position("top")
		ax2.hist(seq_identity,  alpha=0.8, color='#A7D5B9', label='seq_identity', edgecolor='white', 
						linewidth=1.0)
		ax1.set_ylabel('Count', labelpad=10, fontsize=14)
		ax1.set_xlabel('Native sequence recovery', labelpad=10, fontsize=14)
		ax2.set_xlabel('Sequence identity', labelpad=10, fontsize=14)
		ax1.tick_params(axis='x', labelsize=12)
		ax1.tick_params(axis='y', labelsize=12)
		ax2.tick_params(axis='x', labelsize=12)
		if save == 'yes':
			plt.savefig(self.figs_out+'/'+name_fig+'.png', 
									dpi=300, bbox_inches='tight')
			plt.clf()
			plt.show()
		else: 
			plt.show()

	def redes_exposed_for_model(self, sample_numb, plot='yes', save='yes'):
		'''
		Plots a histogram of the exposed input and redesigned residues for a specific model sorted according to hydrophilicity.

		Parameters:
				- sample_numb (int): The index of the model/sample to inspect.
				- plot (str, optional): Whether to plot the data. Default is 'yes'.
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				- residues_list (list): list of redesigned residue numbers in the given model
		'''
		designed_res_indices = mpnn_analysis.find_redes_residues(self)
		list_pdb_files= mpnn_analysis.get_pdbs(self)
		print('plotting data for model', str(list_pdb_files[sample_numb]))
		ref_universe = MDAnalysis.Universe(self.pdb_folder+'/'+str(list_pdb_files[0]), dt=1)
		des_universe = MDAnalysis.Universe(self.pdb_folder+'/'+str(list_pdb_files[sample_numb]), dt=1)
		ref_structure_data = mpnn_analysis.data_universe(ref_universe)
		des_structure_data = mpnn_analysis.data_universe(des_universe)
		
		subset_ref = ref_structure_data.iloc[designed_res_indices[1:][sample_numb-1]]
		subset_des = des_structure_data.iloc[designed_res_indices[1:][sample_numb-1]]

		residues_list = list(np.array(designed_res_indices[1:][sample_numb-1])+1) #residue in pdb=index+1
		print('list of redesigned residues:\n',residues_list,end=',') #prints RESIDUE number at designed positions
		if plot == 'yes': 
			exposed_ref = subset_ref[subset_ref['fraction_sasa'] > 0.2]
			exposed_des = subset_des[subset_des['fraction_sasa'] > 0.2]
			count_byres_ref = exposed_ref.groupby('resname')['res_numb'].count()
			count_byres_des = exposed_des.groupby('resname')['res_numb'].count()
			dict_count_ref = {i:count_byres_ref[i] for i in list(count_byres_ref.index)}
			dict_count_des = {j:count_byres_des[j] for j in list(count_byres_des.index)}
			
			for amino_acid in aas_alphabet.keys():
				if amino_acid not in dict_count_ref.keys():
					dict_count_ref[amino_acid] = 0        
			for amino_acid in aas_alphabet.keys():
				if amino_acid not in dict_count_des.keys():
					dict_count_des[amino_acid] = 0
			count_ref_residues_sorted = {i:dict_count_ref[i] for i in sorted(dict_count_ref)}
			count_des_residues_sorted = {j:dict_count_des[j] for j in sorted(dict_count_des)}
			count_exposed_df = pd.DataFrame([count_ref_residues_sorted, count_des_residues_sorted], index=['ref', 'redes'])
			plot_df = count_exposed_df.T
				
			xticks_gen = [i for i in range(0, 110, 5)]
			max_ref = max(plot_df['ref'])
			max_redes = max(plot_df['redes'])
			if max_ref > max_redes:
				for value in range(len(xticks_gen)):
					if max_ref < xticks_gen[value]:
						upper_limit_ticks = xticks_gen[value+1]
						break
			else: 
				for value in range(len(xticks_gen)):	
					if max_redes < xticks_gen[value]:
						upper_limit_ticks = xticks_gen[value+1]
						break

			order_hydrophilic = []
			order_hydrophobic = []
			for i in polarity_aas['hydrophilic']:
				order_hydrophilic.append(i)
			for j in polarity_aas['hydrophobic']:
				order_hydrophobic.append(j)
			combined_order = sorted(order_hydrophobic)+sorted(order_hydrophilic)
			plot_df = plot_df.reindex(index = combined_order)
			
			fig, axes = plt.subplots(figsize=(10,5), ncols=2, sharey=True) 
			colors = ["lightgrey" for i in list(plot_df.index)]
			fig.tight_layout()
			axes[0].barh(plot_df.index, plot_df['ref'], align='center', color=colors, zorder=10)
			axes[0].set(yticks=plot_df.index, yticklabels=plot_df.index)
			axes[0].invert_xaxis() 
			axes[0].xaxis.set_tick_params(labelsize=14)
			axes[0].set_xticks([i for i in range(0, upper_limit_ticks, 5)])
			axes[0].yaxis.set_tick_params(labelsize=14)
			axes[0].yaxis.tick_left()
			axes[0].set_title('Exposed input residues', color='grey', fontsize=18)
			# axes[0].grid(axis = 'y', alpha=0.3)
			# axes[0].tick_params(tick1On=False)
			axes[1].barh(plot_df.index, plot_df['redes'], align='center', color ='#AA2023', zorder=10) #original #A5BCE5
			axes[1].set_xticks([i for i in range(5, upper_limit_ticks, 5)])
			axes[1].xaxis.set_tick_params(labelsize=14)
			axes[1].yaxis.set_ticks_position('none') 
			axes[1].set_title('Exposed redesigned residues', color='#AA2023', fontsize=18) #original color #8FADE3
			# axes[1].grid(axis = 'y', alpha=0.3)
			plt.gca().invert_yaxis()
			axes[0].spines['top'].set_visible(False)
			axes[0].spines['right'].set_visible(False)
			axes[0].spines['left'].set_visible(False)
			axes[1].spines['top'].set_visible(False)
			axes[1].spines['right'].set_visible(False)
			axes[1].spines['left'].set_visible(False)
			plt.subplots_adjust(wspace=0, top=1.0, bottom=0.1, left=0.18, right=0.95)
			if save == 'yes':
				plt.savefig(self.figs_out+'/'+'count_residues_des_new{}.png'.format(sample_numb), 
									dpi=300, bbox_inches='tight')
				plt.show()
		else:
			return residues_list
		
	def plot_core_size(self, save='yes'):
		'''
		Plots a catplot with the number of buried residues in reference and each of the designed structures.

		Args:
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				pandas.DataFrame: DataFrame containing the count of buried (solvent excluded) residues in the reference and each designed structure.
		'''
		name_plot = 'hydrophobic_core_size'
		list_data_models = mpnn_analysis.get_data(self)
		count_hb = []
		for structure in list_data_models:
			count_hb.append(len(structure[structure['sasa_value'] == 0]))
		annot_wt = ['wt' for i in range(len(count_hb)) if i == 0]
		annot_des = ['des' for i in range(len(count_hb)) if i != 0]
		annot = annot_wt + annot_des
		dict_hb_core = {'count_hb':count_hb,
										'annot':annot}
		df_hb_core = pd.DataFrame(dict_hb_core, columns=dict_hb_core.keys())
		custom_palette = ["#92D1B9", "#A3AFCE", '#BDA3CE', '#CEC4A3', '#EEC3C3']
		sns.catplot(data=df_hb_core, x='annot', y='count_hb', legend=False, linewidth=0.6, palette=custom_palette, s = 7) #change size to any desired value; this will appear pretty big
		plt.ylabel('Count', labelpad=10, fontsize=14)
		plt.xlabel('')
		plt.title('Calculated core size', fontsize=14)   
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)
		if save=='yes':
			plt.savefig(self.figs_out+'/'+name_plot+'.png', 
												dpi=300, bbox_inches='tight')
			plt.clf()
		return df_hb_core
	
	def violin_fraction_exposed(self, save='yes'):
		'''
		Plots a combined stip (scatter points) and violin plot to visualize the fraction of exposed residues in the input and designed structures.

		Parameters:
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				pandas.DataFrame: dataframe containing fractions of input and designed exposed values
		'''
		name_plot = 'frac_exposed_redes'
		ref_annot = ['input' for i in range(len(self.dataframe))]
		redes_annot = ['design' for i in range(len(self.dataframe))]
		fractions = pd.concat([self.dataframe['frac_exposed_ref'], self.dataframe['frac_exposed_redes']], ignore_index=True)
		annots = ref_annot + redes_annot
		group1=self.dataframe['frac_exposed_ref']
		group2=self.dataframe['frac_exposed_redes']
		mw_statistic, mw_pvalue = stats.mannwhitneyu(group1, group2)
		pvalue_formatted=[f'p={mw_pvalue:.2e}']
		pvalue_num = [mw_pvalue]
		
		fraction_exposed_violin = {'fraction_exposed':fractions, 
									'annot':annots}
		df_exposed_violin = pd.DataFrame(fraction_exposed_violin)
		custom_palette=['#A7D5B9','#A3AFCE']
		plotting_parameters = {
		'data':    df_exposed_violin,
		'x':       'annot',
		'y':       'fraction_exposed',
		'palette': custom_palette
		}
		pairs = [('input', 'design')]
		sns.catplot(**plotting_parameters,legend=False, linewidth=0.4,
							kind='strip')
		ax = sns.violinplot(**plotting_parameters, linewidth = 0)
		for violin, alpha in zip(ax.collections, [1.0, 1.0, 0.6, 0.6, 0.6]):
			violin.set_alpha(alpha)
		annotator = Annotator(ax, pairs, **plotting_parameters)
		annotator.set_custom_annotations(pvalue_formatted)
		annotator.set_pvalues(pvalue_num)
		annotator.annotate()
		# ax.spines['left'].set_linewidth(1.5)
		# ax.spines['bottom'].set_linewidth(1.5)
		plt.yticks(fontsize=12)
		plt.xticks(fontsize=12)
		plt.ylabel('Fraction of exposed resdiues', labelpad = 10, fontsize=14)
		plt.xlabel('')
		plt.title('Solvent exposure of designed residues', fontsize=14)
		if save=='yes':
			plt.savefig(self.figs_out+'/'+name_plot+'.png', 
												dpi=300, bbox_inches='tight')
			plt.clf()
		return df_exposed_violin

	def violin_hydrophilic_exposed(self, save='yes'):
		'''
 		Plots a combined stip (scatter points) and violin plot to visualize the fraction of exposed hydrophilic residues in the input and designed structures.

		Args:
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				pandas.DataFrame: DataFrame containing the count of buried (solvent excluded) residues in the reference and each designed structure.
		'''
		name_plot = 'frac_hydrophilic_exposed_redes'
		ref_annot = ['input' for i in range(len(self.dataframe))]
		redes_annot = ['design' for i in range(len(self.dataframe))]
		fractions = pd.concat([self.dataframe['hydrophilic_exposed_ref'], self.dataframe['hydrophilic_exposed_redes']], ignore_index=True)
		annots = []
		annots.extend(ref_annot)
		annots.extend(redes_annot)
		group1=self.dataframe['hydrophilic_exposed_ref']
		group2=self.dataframe['hydrophilic_exposed_redes']
		mw_statistic, mw_pvalue = stats.mannwhitneyu(group1, group2)
		pvalue_formatted=[f'p={mw_pvalue:.2e}']
		pvalue_num = [mw_pvalue]
		
		dict_hydrophilic_violin = {'hydrophilic_frac_exposed':fractions, 
									'annot':annots}
		df_hydrophilic_violin = pd.DataFrame(dict_hydrophilic_violin)
		custom_palette=['#A7D5B9','#A3AFCE']
		plotting_parameters = {
		'data':    df_hydrophilic_violin,
		'x':       'annot',
		'y':       'hydrophilic_frac_exposed',
		'palette': custom_palette
		}
		pairs = [('input', 'design')]
		sns.catplot(**plotting_parameters,legend=False, linewidth=0.4,
											kind='strip')
		ax = sns.violinplot(**plotting_parameters, linewidth = 0, alpha=0.5)
		for violin, alpha in zip(ax.collections, [1.0, 1.0, 0.6, 0.6, 0.6]):
			violin.set_alpha(alpha)
		annotator = Annotator(ax, pairs, **plotting_parameters)
		annotator.set_pvalues(pvalue_num)
		annotator.annotate()
		# ax.spines['left'].set_linewidth(1.5)
		# ax.spines['bottom'].set_linewidth(1.5)
		plt.yticks(fontsize=12)
		plt.xticks(fontsize=12)
		plt.xlabel('')
		plt.title('Hydrophilicity of solvent-exposed designed residues', fontsize=14)
		plt.ylabel('Fraction of exposed hydrophilic resdiues', labelpad = 10, fontsize=14)
		if save=='yes':
			plt.savefig(self.figs_out+'/'+name_plot+'.png', 
												dpi=300, bbox_inches='tight')
			plt.clf()
		return df_hydrophilic_violin
	
	def histogram_redes_res(self, save='yes', selection='exposed'):
		'''
		Plots a bar plot to visualize the frequency of exposed or all designed residues in the designed structures. Y-axis: normalized frequency, x-axis: aa identity

		Parameters:
				- save (str, optional): Whether to save the plot. Default is 'yes'.
				- selection (str, optional): Selection of residues to plot. Options are 'exposed' (default) for exposed residues and 'all' for all designed residues.

		Returns:
				- count_redes_exposed_residues_sorted (dict): A dictionary containing the count of exposed residues sorted by residue name if selection='exposed'. 
				- count_redes_all_residues_sorted (dict): A dictionary containing the count of all designed residues sorted by residue name if selection='all'.
		'''
		list_pdb_files= mpnn_analysis.get_pdbs(self)
		designed_res_indices = mpnn_analysis.find_redes_residues(self)
		list_res_sasas = []
		for model in list_pdb_files:
			universe_per_res = MDAnalysis.Universe(self.pdb_folder+'/'+str(model), dt=1)
			list_res_sasas.append(mpnn_analysis.data_universe(universe_per_res))
		reference_universe = list_res_sasas[0]
		to_compare_universes = list_res_sasas[1:]
		if selection == 'exposed':
			name_plot = 'frequency_exposed_designed_res'
			count_redes_exposed_residues = {}
			for i in range(len(to_compare_universes)):
				exp_condition_redes = to_compare_universes[i].iloc[designed_res_indices[1:][i]]['fraction_sasa'] > 0.2
				condition_redes_subset = to_compare_universes[i].iloc[designed_res_indices[1:][i]][exp_condition_redes]
				exp_subset_redes_count = condition_redes_subset.groupby('resname')['res_numb'].count()
				dict_count_redes = {j:exp_subset_redes_count[j] for j in list(exp_subset_redes_count.index)}
				for amino_acid in list(dict_count_redes.keys()):
					if amino_acid in list(count_redes_exposed_residues.keys()):
						count_redes_exposed_residues[amino_acid] += dict_count_redes[amino_acid]
					else:
						count_redes_exposed_residues[amino_acid] = dict_count_redes[amino_acid]
			for amino_acid in aas_alphabet.keys():
				if amino_acid not in count_redes_exposed_residues.keys():
					count_redes_exposed_residues[amino_acid] = 0
			count_redes_exposed_residues_sorted = {k:count_redes_exposed_residues[k] for k in sorted(count_redes_exposed_residues)}
				
			hydrophobicity_annot = []
			for key in list(count_redes_exposed_residues_sorted.keys()):
				if key in polarity_aas[list(polarity_aas.keys())[0]]:
					hydrophobicity_annot.append(list(polarity_aas.keys())[0])
				else:
					hydrophobicity_annot.append(list(polarity_aas.keys())[1])
			normalized_values = [i/sum(list(count_redes_exposed_residues_sorted.values())) for i in list(count_redes_exposed_residues_sorted.values())]
			dict_freq_redes_exp = {'resname':list(count_redes_exposed_residues_sorted.keys()),
												'freq':normalized_values,
												'hydrophobicity':hydrophobicity_annot}
			df_freq_redes_exp = pd.DataFrame(dict_freq_redes_exp, columns=dict_freq_redes_exp.keys())
			color_hydrophobic = ['#8FADE3', '#CBCD91']
			fig = plt.figure(figsize=(6.4 ,4.8))
			ax = sns.barplot(data=df_freq_redes_exp.sort_values(['hydrophobicity', 'resname']), 
									x='resname', y='freq', hue='hydrophobicity', palette=color_hydrophobic, dodge=False)
			# plt.bar(list(count_redes_exposed_residues_sorted.keys()), normalized_values, color='lightsteelblue', edgecolor='black')
			plt.ylabel('Frequency')
			plt.title('exposed residues designed by pMPNN')
			plt.xticks(rotation=90)
			if save=='yes':
				plt.savefig(self.figs_out+'/'+name_plot+'.png', 
													dpi=300, bbox_inches='tight')
				plt.clf()
			return count_redes_exposed_residues_sorted
		else: 
			name_plot = 'frequency_all_designed_res_excl'
			count_redes_residues = {}
			for i in range(len(to_compare_universes)):
				redes_res = to_compare_universes[i].iloc[designed_res_indices[1:][i]]
				redes_res_count = redes_res.groupby('resname')['res_numb'].count()
				df_freq_redes_exp = {j:redes_res_count[j] for j in list(redes_res_count.index)}
				for amino_acid in list(dict.keys()):
					if amino_acid in list(count_redes_residues.keys()):
						count_redes_residues[amino_acid] += df_freq_redes_exp[amino_acid]
					else:
						count_redes_residues[amino_acid] = df_freq_redes_exp[amino_acid]
			for amino_acid in aas_alphabet.keys():
				if amino_acid not in count_redes_residues.keys():
					count_redes_residues[amino_acid] = 0
			count_redes_all_residues_sorted = {k:count_redes_residues[k] for k in sorted(count_redes_residues)}
			normalized_values = [i/sum(list(count_redes_all_residues_sorted.values())) for i in list(count_redes_all_residues_sorted.values())]
			hydrophobicity_annot = []
			for key in list(count_redes_all_residues_sorted.keys()):
				if key in polarity_aas[list(polarity_aas.keys())[0]]:
					hydrophobicity_annot.append(list(polarity_aas.keys())[0])
				else:
					hydrophobicity_annot.append(list(polarity_aas.keys())[1])

			dict_freq_redes = {'resname':list(count_redes_all_residues_sorted.keys()),
												'freq':normalized_values,
												'hydrophobicity':hydrophobicity_annot}
			df_freq_redes = pd.DataFrame(dict_freq_redes, columns=dict_freq_redes.keys())
			color_hydrophobic = ['#8FADE3', '#CBCD91']
			fig = plt.figure(figsize=(6.4 ,4.8))
			ax = sns.barplot(data=df_freq_redes.sort_values(['hydrophobicity', 'resname']), 
									x='resname', y='freq', hue='hydrophobicity', palette=color_hydrophobic, dodge=False)
			# plt.bar(list(count_redes_residues_sorted.keys()), normalized_values, color='lightsteelblue', edgecolor='black')
			plt.ylabel('Frequency')
			plt.title('designed residues by pMPNN')
			plt.xticks(rotation=90)
			if save=='yes':
				plt.savefig(self.figs_out+'/'+name_plot+'.png', 
													dpi=300, bbox_inches='tight')
				plt.clf()
			return count_redes_all_residues_sorted

	def plot_charges(self, selection='total', save='yes'):
		'''
		Plots a violin plot or a histogram to visualize the total charge or charge of exposed residues in the designed structures.

		Parameters:
				- selection (str, optional): Selection of charges to plot. Options are 'total' (default) for total charge and 'exposed' for charge of exposed residues.
				- save (str, optional): Whether to save the plot. Default is 'yes'.

		Returns:
				- df_charges (DataFrame): DataFrame containing the charges and annotations if selection='exposed'.
				- charge_fullstr_des (list): List containing the total charges of the designed structures if selection='total'.
		'''
		list_pdb_files = mpnn_analysis.get_pdbs(self)
		designed_res_indices =  mpnn_analysis.find_redes_residues(self)
		count_redes_residues_sasa = []
		list_res_sasas = []
		for model in list_pdb_files:
			universe_per_res = MDAnalysis.Universe(self.pdb_folder+'/'+str(model), dt=1)
			list_res_sasas.append(mpnn_analysis.data_universe(universe_per_res))
		reference_universe = list_res_sasas[0]
		to_compare_universes = list_res_sasas[1:]
		if selection == 'exposed':
			exposed_charge_ref = [] 
			exposed_charge_des = []

			for i in range(len(to_compare_universes)):
				charged_condition_ref = reference_universe.iloc[designed_res_indices[1:][i]]['fraction_sasa'] > 0.2
				charged_condition_redes = to_compare_universes[i].iloc[designed_res_indices[1:][i]]['fraction_sasa'] > 0.2
				charge_ref_subset = reference_universe.iloc[designed_res_indices[1:][i]][charged_condition_ref]
				charge_redes_subset = to_compare_universes[i].iloc[designed_res_indices[1:][i]][charged_condition_redes]
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

			ref_annot = ['input' for i in range(len(exposed_charge_ref))]
			redes_annot = ['design' for i in range(len(exposed_charge_des))]
			annots = ref_annot + redes_annot
			charges = exposed_charge_ref + exposed_charge_des
			group1=exposed_charge_ref
			group2=exposed_charge_des
			mw_statistic, mw_pvalue = stats.mannwhitneyu(group1, group2)
			pvalue_formatted=[f'p={mw_pvalue:.2e}']
			pvalue_num = [mw_pvalue]
			dict_charges = {'charges':charges, 
											'annot':annots}
			df_charges = pd.DataFrame(dict_charges)
			name_plot = 'charges_{}_redesigned'.format(selection)
			custom_palette=['#A7D5B9','#A3AFCE']
			plotting_parameters = {
			'data':    df_charges,
			'x':       'annot',
			'y':       'charges',
			'palette': custom_palette
			}
			pairs = [('input', 'design')]
			sns.catplot(**plotting_parameters,legend=False, linewidth=0.4,
												kind='strip')
			ax = sns.violinplot(**plotting_parameters, linewidth = 0, alpha=0.5)
			for violin, alpha in zip(ax.collections, [1.0, 1.0, 0.6, 0.6, 0.6]):
					violin.set_alpha(alpha)
			annotator = Annotator(ax, pairs, **plotting_parameters)
			annotator.set_pvalues(pvalue_num)
			annotator.annotate()
			# ax.spines['left'].set_linewidth(1.5)
			# ax.spines['bottom'].set_linewidth(1.5)
			plt.yticks(fontsize=12)
			plt.xticks(fontsize=12)
			plt.xlabel('')
			plt.title('Total charge of solvent-exposed redesigned residues', fontsize=14)
			plt.ylabel('Charge', labelpad = 10, fontsize=14)
			if save=='yes':
				plt.savefig(self.figs_out+'/'+name_plot+'.png', 
													dpi=300, bbox_inches='tight')
				plt.clf()
			return df_charges
		if selection=='total': 
			start = time.time()
			total_charge_ref = []
			total_charge_des = []
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
				charge_indiv = []
				for res_name in des['resname']:
					if res_name in list(charged_aas.values())[0]:
						charge_indiv.append(1)
					if res_name in list(charged_aas.values())[1]:
						charge_indiv.append(-1)
				charge_fullstr_des.append(sum(charge_indiv))
			name_plot = '{}_charges'.format(selection)
			# didn't benchmark this part too thoroughly
			if charge_fullstr_ref < max(charge_fullstr_des):
				min_plot = int(max(charge_fullstr_des))+2
			else:
				min_plot = int(charge_fullstr_ref)+2
			if min(charge_fullstr_des) < charge_fullstr_ref:
				max_plot = int(min(charge_fullstr_des))-2
			else:
				max_plot = int(charge_fullstr_ref)-2
			fig = plt.figure(figsize=(6.4 ,4.8))
			ax = plt.hist(charge_fullstr_des, bins=9, color='lightsteelblue', edgecolor='white', label='net_charge_designs', alpha=0.8)
			plt.axvline(charge_fullstr_ref, alpha=0.8, linestyle ="--", color='salmon', label='net_charge_input')
			plt.xlim(max_plot, min_plot)
			plt.title('Net charge of designs', fontsize=14)
			plt.ylabel('Count', labelpad=10, fontsize=14)
			plt.xlabel('Charge', fontsize=14)
			plt.xticks(fontsize=12)
			plt.yticks(fontsize=12)
			plt.legend(loc='best',prop = {"size": 8})
			if save=='yes':
				plt.savefig(self.figs_out+'/'+name_plot+'.png', 
													dpi=300, bbox_inches='tight')  
				plt.clf()
			end = time.time()
			print('--- finished after:', end-start, '---')
			return charge_fullstr_des
		
	def heatmap_exchanged_residues(self, save='yes', sort_hydrophobic=False, non_wt_only = True):
		'''
		Plots a 20x20 (des identity = f(orig identity)) heatmap to visualize the exchanged residues between original and designed sequences

		Parameters:
				save (str, optional): Whether to save the plot. Default is 'yes'.
				sort_hydrophobic (bool, optional): Whether to sort residues by hydrophobicity. Default is False.
				non_wt_only (bool, optional): Whether to consider only non-wild-type (non-wt) exchanges. Default is True.

		Returns:
				pd.DataFrame: A DataFrame containing the normalized heatmap data.

		The function generates a heatmap to visualize the exchanged residues between the original and designed sequences.
		It counts the exchanges of each residue type between the original and designed sequences and plots a heatmap based on the frequency of exchanges.

		If non_wt_only is True, the heatmap only includes exchanges between non-wild-type residues, leaving the diagonal empty.
		If sort_hydrophobic is True, residues are sorted by hydrophobicity in the heatmap.
		The plot can be saved as a PNG image if save is set to 'yes'.
		'''
		dict_designed = {}
		index_redes = []
		for amino_acid in letters:
			dict_designed[amino_acid] = [[i,0] for i in letters]
		with open(self.fasta) as fasta_input:
			lines = fasta_input.readlines()
		all_seqs = [lines[i].strip('\n') for i in range(len(lines)) if i %2 != 0]
		orig_input = all_seqs[0]
		for seq in all_seqs[1:]:
			for index in range(len(seq)):
				if seq[index] != orig_input[index]:
					index_redes.append(index)
					for entry in range(len(dict_designed[orig_input[index]])):
						if seq[index] == dict_designed[orig_input[index]][entry][0]: #iterate over the residues until you find a correct exchanged one
							dict_designed[orig_input[index]][entry][1] += 1 #when found add 1
		if non_wt_only == True:
			count_exchanged = []
			for key in dict_designed.keys():
				count_exchange_iter = []
				for entries in dict_designed[key]:
					count_exchange_iter.append(entries[1])
				count_exchanged.append(count_exchange_iter)
			
			dict_count_only = {}
			for iteration in range(len(count_exchanged)):
				dict_count_only[letters[iteration]] = count_exchanged[iteration]
			
			heatmap_exchanged_res = pd.DataFrame(dict_count_only, index=list(dict_count_only.keys()))
			heatmap_exchanged_res.replace(0, np.nan, inplace=True)
			max_value = max(list(heatmap_exchanged_res.max()))
			normalized_heatmap = heatmap_exchanged_res.div(heatmap_exchanged_res.sum(axis=0), axis=1)
			if sort_hydrophobic == True:
				name_fig = 'heatmap_redesigned_residues_sort_hydrophobic'
				order_hydrophilic = []
				order_hydrophobic = []

				for i in polarity_aas['hydrophilic']:
					if i in aas_alphabet.keys():
						i = aas_alphabet[i]
						order_hydrophilic.append(i)

				for j in polarity_aas['hydrophobic']:
					if j in aas_alphabet.keys():
						j = aas_alphabet[j]
						order_hydrophobic.append(j)
				combined_order = sorted(order_hydrophobic)+sorted(order_hydrophilic)
				heatmap_exchanged_res = heatmap_exchanged_res[combined_order]
				heatmap_exchanged_res = heatmap_exchanged_res.reindex(combined_order)
				plt.figure(figsize=(16,9))
				ax = sns.heatmap(heatmap_exchanged_res, cmap="PuBu",linewidths = 0.5, square= True, linecolor='lightgrey', 
								vmin=0, vmax=max_value, cbar_kws={"pad": 0.01}) #rocket
				ax.xaxis.tick_top()
				ax.xaxis.set_label_position('top') 
				cbar = ax.collections[0].colorbar
				cbar.ax.tick_params(labelsize=18)
				plt.xlabel('Original residue', fontsize=20, labelpad=20)
				plt.ylabel('Designed residue', fontsize=20, labelpad=20)
				ax.collections[0].colorbar.set_label("Number of exchanges", labelpad=20)
				ax.figure.axes[-1].yaxis.label.set_size(20)
				plt.yticks(fontsize=18)
				plt.xticks(fontsize=18)
				if save =='yes':
						plt.savefig(self.figs_out+'/'+name_fig+'.png', 
												dpi=300, bbox_inches='tight')
						plt.clf()
			else:
				name_fig = 'heatmap_redesigned_residues'
				# normalized_heatmap = heatmap_exchanged_res/max_value
				plt.figure(figsize=(16,9))
				ax = sns.heatmap(heatmap_exchanged_res, cmap="PuBu",linewidths = 0.5, square= True, linecolor='lightgrey', 
												vmin=0, vmax=max_value, cbar_kws={"pad": 0.01}) #PuBu
				ax.xaxis.tick_top()
				ax.xaxis.set_label_position('top') 
				cbar = ax.collections[0].colorbar
				cbar.ax.tick_params(labelsize=18)
				plt.xlabel('Original residue', fontsize=24, labelpad=20)
				plt.ylabel('Designed residue', fontsize=24, labelpad=20)
				ax.collections[0].colorbar.set_label("Number of exchanges", labelpad=20)
				ax.figure.axes[-1].yaxis.label.set_size(22)
				plt.yticks(fontsize=20)
				plt.xticks(fontsize=20)
				# plt.show()
				if save =='yes':
					plt.savefig(self.figs_out+'/'+name_fig+'.png', 
												dpi=300, bbox_inches='tight')
					plt.clf()
		else:
			for seq in all_seqs[1:]:
				for index in range(len(seq)):
					index_redes.append(index)
					for entry in range(len(dict_designed[orig_input[index]])):
						if seq[index] == dict_designed[orig_input[index]][entry][0]: #iterate over the residues until you find a correct exchanged one
							dict_designed[orig_input[index]][entry][1] += 1 #when found add 1
			count_exchanged = []
			for key in dict_designed.keys():
				count_exchange_iter = []
				for entries in dict_designed[key]:
					count_exchange_iter.append(entries[1])
				count_exchanged.append(count_exchange_iter)
						
			dict_count_only = {}
			for iteration in range(len(count_exchanged)):
				dict_count_only[letters[iteration]] = count_exchanged[iteration]
			# return dict_count_only
			heatmap_exchanged_res = pd.DataFrame(dict_count_only, index=list(dict_count_only.keys()))
			heatmap_exchanged_res.replace(0, np.nan, inplace=True)
			max_value = max(list(heatmap_exchanged_res.max())) #if absolute values scale is needed
			normalized_heatmap = heatmap_exchanged_res.div(heatmap_exchanged_res.sum(axis=0), axis=1)
			if sort_hydrophobic == True:
				name_fig = 'heatmap_all_residues_sort_hydrophobic'
				order_hydrophilic = []
				order_hydrophobic = []

				for i in polarity_aas['hydrophilic']:
					if i in aas_alphabet.keys():
						i = aas_alphabet[i]
						order_hydrophilic.append(i)

				for j in polarity_aas['hydrophobic']:
					if j in aas_alphabet.keys():
						j = aas_alphabet[j]
						order_hydrophobic.append(j)

				combined_order = sorted(order_hydrophobic)+sorted(order_hydrophilic)
				normalized_heatmap = normalized_heatmap[combined_order]
				normalized_heatmap = normalized_heatmap.reindex(combined_order)
				plt.figure(figsize=(16,9))
				ax = sns.heatmap(normalized_heatmap, cmap="mako_r",linewidths = 0.5, square= True, linecolor='lightgrey', 
											vmin=0, vmax=1, cbar_kws={"pad": 0.01}) #PuBu
				ax.xaxis.tick_top()
				ax.xaxis.set_label_position('top') 
				cbar = ax.collections[0].colorbar
				cbar.ax.tick_params(labelsize=18)
				plt.xlabel('Original residue', fontsize=24, labelpad=20)
				plt.ylabel('Designed residue', fontsize=24, labelpad=20)
				ax.collections[0].colorbar.set_label("Frequency of exchange", labelpad=20)
				ax.figure.axes[-1].yaxis.label.set_size(22)
				plt.yticks(fontsize=20)
				plt.xticks(fontsize=20)
				if save =='yes':
					plt.savefig(self.figs_out+'/'+name_fig+'.png', 
											dpi=300, bbox_inches='tight')
					plt.clf()
			else: 
				name_fig = 'heatmap_redesigned_residues_with_wt_norm_bycol'
				plt.figure(figsize=(16,9))
				ax = sns.heatmap(normalized_heatmap, cmap="mako_r",linewidths = 0.5, square= True, linecolor='lightgrey', 
											vmin=0, vmax=1, cbar_kws={"pad": 0.01}) #PuBu
				ax.xaxis.tick_top()
				ax.xaxis.set_label_position('top') 
				cbar = ax.collections[0].colorbar
				cbar.ax.tick_params(labelsize=18)
				plt.xlabel('Original residue', fontsize=20, labelpad=20)
				plt.ylabel('Designed residue', fontsize=20, labelpad=20)
				ax.collections[0].colorbar.set_label("Frequency of exchange", labelpad=20)
				ax.figure.axes[-1].yaxis.label.set_size(20)
				plt.yticks(fontsize=18)
				plt.xticks(fontsize=18)
				if save =='yes':
					plt.savefig(self.figs_out+'/'+name_fig+'.png', 
											dpi=300, bbox_inches='tight')
					plt.clf()
		return normalized_heatmap

	def heatmap_allseq_allres(self, save='yes', sort_hydrophobic=False):
		'''
		Plots a heatmap showing the occurrences of designed amino acids at each position of the full sequence of the input and plots corresponding native amino acid identities.

		Parameters:
				- save (bool, optional): Indicates whether to save the plot. Default is 'yes'.
				- sort_hydrophobic (bool, optional): Specifies whether to sort residues by hydrophobicity. Default is False

		Returns:
				None
		'''
		orig_pssm, norm_msa_pssm = mpnn_analysis.gen_pssm(self, heatmap='yes')
		cmap_dict = {0: '#F28D90'}
		cmap = ListedColormap([cmap_dict[i] for i in range(len(cmap_dict))])
		if sort_hydrophobic == True:
			name_fig = 'fullseq_heatmap_exchanges_sorted'
			order_hydrophilic = []
			order_hydrophobic = []
			for i in polarity_aas['hydrophilic']:
				if i in aas_alphabet.keys():
					i = aas_alphabet[i]
					order_hydrophilic.append(i)
			for j in polarity_aas['hydrophobic']:
				if j in aas_alphabet.keys():
					j = aas_alphabet[j]
					order_hydrophobic.append(j)
			combined_order = sorted(order_hydrophobic)+sorted(order_hydrophilic)
			norm_msa_pssm_order = norm_msa_pssm[combined_order]
			# seq_length = len(norm_msa_pssm_order.T)
			# min_width, max_width = 10, 50
			# min_height, max_height = 3, 10
			# # Calculate width and height based on sequence length
			# width = min(max(min_width * seq_length, min_width), max_width)
			# height = min(max(min_height * seq_length, min_height), max_height)
			# plt.figure(figsize=(width, height))
			plt.figure(figsize=(50, 5)) #gotta find a better way to control the size
			ax = sns.heatmap(norm_msa_pssm_order.T, cmap='Blues_r', #crest_r
											cbar_kws={"orientation": "vertical", "pad": 0.002})
			ax.collections[0].colorbar.set_label("Frequency at position")
			ax.figure.axes[-1].yaxis.label.set_size(16)
			plt.xlabel('Residue number', fontsize=16, labelpad=20)
			plt.ylabel('Amino acid', fontsize=16, labelpad=20)#PuBu
			ax1 = sns.heatmap(orig_pssm[combined_order].T, cmap=cmap, alpha=0.7,
							cbar=False) #if orig_pssm has 0 values instead of NaN, then the background in the overlaied heatmap will appear
		else:
			name_fig = 'fullseq_heatmap_exchanges'
			print(norm_msa_pssm)
			seq_length = len(norm_msa_pssm.T)
			min_width, max_width = 10, 50
			min_height, max_height = 5, 10
			# Calculate width and height based on sequence length
			# width = min(max(min_width * seq_length, min_width), max_width)
			# height = min(max(min_height * seq_length, min_height), max_height)
			# plt.figure(figsize=(width, height))
			plt.figure(figsize=(50, 5))
			ax = sns.heatmap(norm_msa_pssm.T, cmap='Blues_r', #crest_r
											cbar_kws={"orientation": "vertical", "pad": 0.002})
			ax.collections[0].colorbar.set_label("Frequency at position")
			ax.figure.axes[-1].yaxis.label.set_size(16)
			plt.xlabel('Residue number', fontsize=16, labelpad=20)
			plt.ylabel('Amino acid', fontsize=16, labelpad=20)#PuBu
			ax1 = sns.heatmap(orig_pssm.T, cmap=cmap, alpha=0.7,
							cbar=False) #if orig_pssm has 0 values instead of NaN, then the background in the overlaied heatmap will appear
		if save =='yes':
			plt.savefig(self.figs_out+'/'+name_fig+'.png', 
									dpi=300, bbox_inches='tight')
			plt.clf()