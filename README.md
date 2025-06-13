# MSc project - Jonas Riber Jørgensen
## Clonal hematopoiesis variant calling workflow

This repository contains most of the code used for the CH variant calling and analysis.

The workflow.py file utilizes gwf to run the jobs on the GenomeDK cluster.

#### Disclaimer 
Since the data contains potentially sensitive personalized information, it is not available here.
Additionally, not all the scripts are made available since they might contain internal patient IDs or similar.
Any paths and functions that might contain or hint at sensitive information have been removed; hence, the code will not run.



### Variant calling
Runs the caller scripts in the following order:
- pysam_mnv_all_mm.py
- merge_call_files2.sh
- apply_blacklist.sh
- pysam_consensus_depth.py
- merge_depth_files.sh
- merge_calls_and_depth.sh
- call_somatic_from_snps.sh
- filtercalls.R
- exclude_gustav_snpcalls.sh
- annotate_calls.sh
- whitelist_filtering.R

- Manually run Rscript scripts/create_variant_patient_table.R results/annotated_samples/ data/0_positions/
- pysam_consensus_depth (for non mismatched positions)
- merge_depth_files.sh

- Manually run Rscript scripts/variant_patient_table2.R


### Analysis scripts
- 1visualize_calls.R
- 2visualize_filtered_and_annotated.R
- mut_signature.R
