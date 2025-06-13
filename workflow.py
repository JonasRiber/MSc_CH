### Workflow for Calling and Annotating SNPs for CH
### Author: Jonas Riber Jørgensen
### Date: 31-01-2025
### Inputs: bam files


import glob # find files and directories
import os # interact with os to find home directory '~'
import re # RegEx module
import pandas as pd
import ast #convert args for python 
from gwf import Workflow # gwf workflow module

gwf = Workflow()
ncores = 4

wd = "/path/to/working/directory"

# Reference
input_reference = '/path/to/reference_genome.fna'

# find all cfdna samples in main_data
bamfiles_cfdna = glob.glob('/path/to/data/files/*cfdna*.bam')
# find all buffycoat samples in main_data
bamfiles_bcdna = glob.glob('/path/to/data/files/*buffycoat*.bam')


# create empty dictionary. should hold patient ID as key and a dict as value
file_dict = {}

# if cfdna samples exist
# iteratively fill patientID as keys of file_dict and
# new directory containing files of type .bam
for f in bamfiles_cfdna:
	m = re.search(r"(?<=main_data/).*?(?=/)", f).group()
	if m not in file_dict.keys():
		file_dict[m] = {'bam':[]}
	file_dict[m]['bam'].append(f)

# iteratively fill buffycoat-dna bamfiles matching patient id
processed_sampleids = []
for f in bamfiles_bcdna:
	m = re.search(r"(?<=main_data/).*?(?=/)", f).group()
	sampleID = re.search(r"SAMPLE_ID_(cfdna|buffycoat)", f).group()
	if (m in file_dict.keys()) and (sampleID not in processed_sampleids):
		file_dict[m]['bam'].append(f)
		processed_sampleids.append(sampleID)




# check how many sample for a patient has ctDNA fraction
def no_ctDNA_samples(ID):
	"""
	function to load meta data and count the number of samples with fraction of ctDNA within
	returns 2 values: number of samples with ctDNA fraction and total number of samples
	"""
	patient_file = "/path/to/ctDNA_detection_file"	
	ID_col = "biobankID"
	sample_id_col = 'sampleID_MOMA'
	cfDNA_fraction_col = "TUMOR_FRACTION_spec"

	df = pd.read_csv(patient_file, delimiter=";", dtype={ID_col: str, sample_id_col: str})
	df.columns = df.columns.str.strip()
	df[ID_col] = df[ID_col].astype(str).str.strip()
	df[cfDNA_fraction_col] = pd.to_numeric(df[cfDNA_fraction_col], errors='coerce')

	#reformat ID - biobankid
	target = re.search(r'C(\d+)', ID).group(1)
	
	samp_count = 0
	ct_count = 0
	ct_sample_id = []
	for i in range(len(df)):
		x = target[1:]==df[ID_col][i]
		if x == True: #match id
			samp_count += 1
			if df[cfDNA_fraction_col][i] > 0: #if fraction above 0 then count ctdna sample
				ct_sample_id.append(df[sample_id_col][i])
				ct_count += 1
				
	return ct_count, samp_count, ct_sample_id



processed_patients = {}

# regions to run on - bedfile created by create_driver_bed.py
#regions_bedfile = "data/misc/driver_genes_and_10_percent.bed" # <-68 driver genes + 10% of genome
regions_bedfile = "data/misc/driver_genes_no_chrX.bed" # <- only 68 driver genes

#### test samples ####
pilot_patients = ['PATIENT_IDs']

samples_required = 10
patient_counter = 0
sample_counter = 0
# https://gwf.app/guide/
# gwf target for each patientID
# collect all samples per patient and combine
for patientID in file_dict:
	#check if ctDNA fraction in samples 
	ctDNA_samples = no_ctDNA_samples(patientID)	
	if patientID in pilot_patients: # for testing only include pilot patients
		continue

	
	#number of samples without ctDNA fraction
	if ctDNA_samples[1] - ctDNA_samples[0] >= samples_required: 
		#remove samples with ctDNA fraction
		bam_files = file_dict[patientID]['bam']
		file_dict[patientID]['bam'] = [bam for bam in bam_files if not any(sample_id in bam for sample_id in ctDNA_samples[2])]
		number_of_samples = len(file_dict[patientID]['bam'])
	else:
		continue

	number_of_samples = len(file_dict[patientID]['bam'])
	
	patient_counter += 1
	sample_counter += number_of_samples
	print(f'pt nr: {patient_counter}\t ID: {patientID}\t nr samples: {number_of_samples}')
	
	
	#######################
	### Call mismatches ###
	# Use pysam_mnv to call mismatches in bam files
	call_files = [] 		# store filepaths to samples for current spatientID
	for sample in file_dict[patientID]['bam']:
		sampleID = re.search(r"SAMPLE_ID_(cfdna|buffycoat)", sample).group()

		call_output_folder = f'data/{patientID}/'
		os.makedirs(call_output_folder, exist_ok=True) 	# ensure that folder exists, if not make it 
		call_output = f'data/{patientID}/{sampleID}_mismatches.tsv'

		gwf.target(
			name=f'Getting_mismatches_for_{sampleID}',
			inputs=[sample],
			outputs=[call_output],
			cores=ncores, # cores to reserve
			memory='4G', # memory to allocate
			walltime='00:30:00', # time to reserve
			account='ACCOUNT') << f"""
			python scripts/pysam_mnv_all_mm.py {sample} {call_output} {regions_bedfile}
			"""
		call_files.append(call_output)
	
	############################
	### Merge the mismatches ###
	# Use merge_call_files to merge the previous calls together
	# subtracts positions from blackslists
	merged_mis_output = f"data/2merged_samples/{patientID}_merged_mismatches.tsv" 
	gwf.target(
		name=f'Merging_samples_for_{patientID}',
		inputs=[call_files],
		outputs=[merged_mis_output],
		cores=ncores, # cores to reserve
		memory='4G', # memory to allocate
		walltime='00:30:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/merge_call_files2.sh {merged_mis_output} {patientID} {" ".join(call_files)}
		""" 
	
	#######################
	### apply blacklist ###
	blacklist_mis_output = f"data/2merged_samples/{patientID}_w_blacklist.tsv" 
	gwf.target(
		name=f'Apply_blacklists_for_{patientID}',
		inputs=[merged_mis_output],
		outputs=[blacklist_mis_output],
		cores=ncores, # cores to reserve
		memory='12G', # memory to allocate
		walltime='00:30:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/apply_blacklist.sh {merged_mis_output} > {blacklist_mis_output}
		""" 
	
	#####################
	### Collect depth ###
	# Uses the pysam_consensus_depth and the called mismatches to get the depth of all the positions
	depth_files = []
	for sample in file_dict[patientID]['bam']:
		sampleID = re.search(r"SAMPLE_ID_(cfdna|buffycoat)", sample).group()
		
		depth_output = f'data/{patientID}/{sampleID}_depth.tsv'
		gwf.target(
			name=f'Getting_depth_for_{sampleID}',
			inputs=[sample, blacklist_mis_output],
			outputs=[depth_output],
			cores=ncores, # cores to reserve
			memory='1G', # memory to allocate
			walltime='08:00:00', # time to reserve
			account='ACCOUNT') << f""" 
			python scripts/pysam_consensus_depth.py {sample} {blacklist_mis_output} {depth_output}
			"""
		depth_files.append(depth_output)
	
	#############################
	### merge the depth files ###
	# Use the merge_depth_files to get the accumulative depth of all the samples
	merged_depth_output = f"data/2merged_samples/{patientID}_merged_depth.tsv" 
	gwf.target(
		name=f'Merging_depth_files_for_{patientID}',
		inputs=[depth_files],
		outputs=[merged_depth_output],
		cores=ncores, # cores to reserve
		memory='2G', # memory to allocate
		walltime='04:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/merge_depth_files.sh {merged_depth_output} {" ".join(depth_files)}
		""" 


	##################################
	### merge mismatches and depth ###
	# use awk script to combine the mismatches and depth, followed by VAF calculation
	merged_output = f'data/3merged_calls_and_depth/{patientID}_merged_mnv_calls_depth.tsv'
	gwf.target(
		name=f'Merging_calls_and_depth_for_{patientID}',
		inputs=[merged_mis_output, merged_depth_output],
		outputs=[merged_output],
		cores=ncores, # cores to reserve
		memory='3G', # memory to allocate
		walltime='04:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/merge_calls_and_depth.sh {blacklist_mis_output} {merged_depth_output} {merged_output}
		""" 

	####################
	### call somatic ###
	# use the call_somatic_from_snps.sh to get most extreme expected alternative count
	p_value = "1e-10"
	somatic_output = f'data/4somatic_calls/{patientID}_somatic_calls_w_depth.tsv'
	gwf.target(
		name=f'Add_somatic_col_to_{patientID}',
		inputs=[merged_output],
		outputs=[somatic_output],
		cores=ncores, # cores to reserve
		memory='1G', # memory to allocate
		walltime='01:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/call_somatic_from_snps.sh {merged_output} {p_value} > {somatic_output}
		""" 
	
	
	###################################
	### R visualization per patient ###
	# Use r to make overview plots per patient
	plot_output = f'results/plots/{patientID}'
	gwf.target(
		name = f'SNP_visualization_for_{patientID}',
		inputs = [somatic_output],
		outputs = [f'{plot_output}/{patientID}_1.1VAF_distribution.png'], 	# Only checks if the first plot is changed atm.
		cores = ncores,
		memory = '8G',
		walltime = '01:00:00',
		account ='ACCOUNT') << f"""
		Rscript scripts/analysis/1visualize_calls.R {somatic_output} {plot_output} 
	 	"""	


	###########################
	### Filtering the calls ###
	#uses the Rscript filtercalls.R to set minimum depth, minimum alternative count, max depth,
	#remove multiallelic positions
	filtered_output = f'data/4somatic_calls/{patientID}_filtered_calls.tsv'
	gwf.target(
		name = f'Filtering_for_{patientID}',
		inputs = [somatic_output],
		outputs = [filtered_output],
		cores = ncores,
		memory = '4G',
		walltime = '04:00:00',
		account ='ACCOUNT') << f"""
		Rscript scripts/filtercalls.R {somatic_output} {filtered_output} 
	 	"""	

	###############################
	### exclude recurrent snps  ###
	# Use snpcalls done by Gustav as blacklist 
	exclude_output = f'data/4somatic_calls/{patientID}_exclude_calls.tsv'
	gwf.target(
		name=f'Excluding_snps_from_blacklist_for_{patientID}',
		inputs=[filtered_output],
		outputs=[exclude_output],
		cores=ncores, # cores to reserve
		memory='3G', # memory to allocate
		walltime='04:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/exclude_gustav_snpcalls.sh {filtered_output} {patientID} {exclude_output}
		""" 

	######################
	### annotate calls ###
	# use annotate_calls with ANNOVAR to annotate for regional function, gene and exonic function (reformats to a 1-based indexing)
	annotated_output = f'data/5annotated_calls/{patientID}_annotated.tsv'
	gwf.target(
		name=f'Annotating_for_{patientID}',
		inputs=[exclude_output],
		outputs=[annotated_output],
		cores=ncores, # cores to reserve
		memory='5G', # memory to allocate
		walltime='01:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/annotate_calls.sh {exclude_output} {annotated_output}
		""" 


	#################################################
	### Visualization after annotation and filter ###
	# Use r to continue overview plots per patient
	plot_output = f'results/plots/{patientID}'
	gwf.target(
		name = f'Visualization_after_filter_and_annotation_{patientID}',
		inputs = [annotated_output],
		outputs = [f'{plot_output}/{patientID}_6.1somatic_vs_germ_VAF_distribution_filtered.png'], 	# Only checks if the first plot is changed atm.
		cores = ncores,
		memory = '8G',
		walltime = '01:00:00',
		account ='ACCOUNT') << f"""
		Rscript scripts/analysis/2visualize_filtered_and_annotated.R {annotated_output} {plot_output} 
	 	"""	

	###################################
	### match variants to whitelist ###
	# Rscript to make a now column where CH status is marked
	# creates a unfiltered marked file, and one filtered for only the CH variants
	CHstatus_output = f'results/annotated_samples/'
	gwf.target(
		name = f'CH_marking_based_on_whitelist_for_{patientID}',
		inputs = [annotated_output],
		outputs = [f'{CHstatus_output}/{patientID}_CH_marked_full.tsv', f'{CHstatus_output}/{patientID}_CH_marked_filtered.tsv'], 	# Only checks if the first plot is changed atm.
		cores = ncores,
		memory = '8G',
		walltime = '01:00:00',
		account ='ACCOUNT') << f"""
		Rscript scripts/whitelist_filtering.R {annotated_output} {CHstatus_output} 
	 	"""	

	
	### fragment length

	#add patient to processed list
	processed_patients[patientID] = number_of_samples

#################################
### Add info to metadata file ###
processed_patients = repr(processed_patients)
processed_patients = '"' + processed_patients + '"'

metadata_file = "backup/data/Metadata.tsv"
gwf.target(
	name = 'Add_metadata_stats',
	inputs = [filtered_output],
	outputs = [metadata_file],
	cores = ncores,
	memory = '1G',
	walltime = '01:00:00',
	account ='ACCOUNT') << f"""
	python scripts/run_metadata_scripts.py {processed_patients}
	"""	

###########################################
### Mutational spectrum on all variants ###
# Rscript for genrating muts signature plots
all_variant_mut_plots = f'results/all_variants_mutsig/'
gwf.target(
	name = f'Plotting_mut_signature_for_all_variants',
	inputs = [CHstatus_output],
	outputs = [f'{all_variant_mut_plots}/all_variants_ratio_CHpos.png', f'{all_variant_mut_plots}/all_variants_trinucleotide_diff_w_cossim.png'], 	# Only checks if the first plot is changed atm.
	cores = ncores,
	memory = '8G',
	walltime = '02:00:00',
	account ='ACCOUNT') << f"""
	Rscript scripts/analysis/mut_signature.R {CHstatus_output} {all_variant_mut_plots} 
	"""	



### run before the below jobs work!!!! 
### Rscript scripts/create_variant_patient_table.R results/annotated_samples/ data/0_positions/
 
completed_jobs = []
### rerunning depth caller on 0_positions
for patientID in file_dict:
	#check if ctDNA fraction in samples 
	ctDNA_samples = no_ctDNA_samples(patientID)	
	#print(patientID, ctDNA_samples)
	if patientID in pilot_patients: # for testing only include pilot patients
		continue

	
	#number of samples without ctDNA fraction
	if ctDNA_samples[1] - ctDNA_samples[0] >= samples_required: 
		#remove samples with ctDNA fraction
		bam_files = file_dict[patientID]['bam']
		file_dict[patientID]['bam'] = [bam for bam in bam_files if not any(sample_id in bam for sample_id in ctDNA_samples[2])]
		number_of_samples = len(file_dict[patientID]['bam'])
	else:
		continue

	

	zero_pos_file = f"data/0_positions/{patientID}_0_positions.tsv"

	#####################################
	### Collect depth for 0_positions ###
	# Uses the pysam_consensus_depth and the called mismatches to get the depth of all the positions
	depth_files_zero_pos = []
	for sample in file_dict[patientID]['bam']:
		sampleID = re.search(r"SAMPLE_ID_(cfdna|buffycoat)", sample).group()
		
		depth_output = f'data/0_positions/single_samples/{patientID}_{sampleID}_0pos_depth.tsv'
		gwf.target(
			name=f'Getting_0pos_depth_for_{sampleID}',
			inputs=[sample],
			outputs=[depth_output],
			cores=ncores, # cores to reserve
			memory='1G', # memory to allocate
			walltime='1:00:00', # time to reserve
			account='ACCOUNT') << f"""
			python scripts/pysam_consensus_depth.py {sample} {zero_pos_file} {depth_output}
			"""
		depth_files_zero_pos.append(depth_output)
	

	#############################
	### merge the depth files ###
	# Use the merge_depth_files to get the accumulative depth of all the samples
	merged_depth_output = f"data/0_positions/merged_depths/{patientID}_merged_0pos_depth.tsv" 
	gwf.target(
		name=f'Merging_0pos_depth_files_for_{patientID}',
		inputs=[depth_output],
		outputs=[merged_depth_output],
		cores=ncores, # cores to reserve
		memory='2G', # memory to allocate
		walltime='04:00:00', # time to reserve
		account='ACCOUNT') << f"""
		bash scripts/merge_depth_files.sh {merged_depth_output} {" ".join(depth_files_zero_pos)}
		""" 
	
	completed_jobs.append(merged_depth_output)


###
print(f"Total samples used: {sample_counter}")



### Rscript scripts/variant_patient_table2.R 
# runs the variant tests

	

#####################################################


