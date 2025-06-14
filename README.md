# MSc project - Jonas Riber Jørgensen
## Clonal hematopoiesis variant calling workflow

This repository contains most of the code used for the CH variant calling and analysis.

The workflow.py file utilizes gwf to run the jobs on the GenomeDK cluster.

### Disclaimer 

Since the data contains potentially sensitive personalized information, it is not available here. 
Additionally, not all the scripts are made available since they might contain internal patient IDs or similar. 
Any paths and functions that might contain or hint at sensitive information have been removed; hence, the code will not run. 



### Variant calling
The workflow is run by being located in the directory with the workflow.py file and using:

gwf run

The workflow Runs the caller scripts in the following order:

- **pysam_mnv_all_mm.py**           (Input single samples)
- **merge_call_files2.sh**          (merges the mismatches in single samples)
- **apply_blacklist.sh**           (applies the various blacklists)
- **pysam_consensus_depth.py**      (Input single samples and mismatch positions)
- **merge_depth_files.sh**          (merges depth files)
- **merge_calls_and_depth.sh**      (merges the mismatches corosponding depth)
- **call_somatic_from_snps.sh**     (apply somatic site testing)
- **filtercalls.R**                 (apply filters)
- **exclude_gustav_snpcalls.sh**    (exclude additional germline snps)
- **annotate_calls.sh**             (Use ANNOVAR to annotate variants)
- **whitelist_filtering.R**         (Match variants to driver list)

- Manually run Rscript scripts/create_variant_patient_table.R results/annotated_samples/ data/0_positions/
  (Combine all the positions into one fill)
  
- **pysam_consensus_depth**     (collect depth for all positions including non-mismatched ones - input single samples)
- **merge_depth_files.sh**      (merge depth samples)

- Manually run Rscript scripts/variant_patient_table2.R
  (Create variant patient matrix)


### Analysis scripts
Some of the included analysis scripts used to generate some of the plots in the report

- **1visualize_calls.R**                      (creates summary plots right after the mismatch caller)
- **2visualize_filtered_and_annotated.R**     (creates summary plots before the filtering script)
- **dataset_overview.R**                      (dataset timeline plot)
- **mut_signature.R**                         (Creates the mutational spectra from the called output)
- **pheat_map.R**                             (Creates the heatmap plots)
- **pval_analysis.R**                           (pvalue analysis and adjusted pvals)
- **stepwise_depth.R**                        (Summary depth distributions through out the workflow, and additional bimodal investigation)
- **variant_analysis.R**                        (variant count plots)  
