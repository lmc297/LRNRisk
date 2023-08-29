# LC's og LRNrisk prototype
#python lrnrisk_prototype.py --gtdb gtdb_example/gtdbtk.bac120.summary.tsv --virulence pima_example/pima_output/vfdb_core.fasta.tabular --amr pima_example/pima_output/AMR_pima_md_2023_02_02.fasta.tabular --blacklist lrnrisk_prototype_blacklist.tsv --vf_distribution vfdb_dist_beta_80p80q.tsv --amr_distribution pima_dist_beta_80p80q.tsv --output test
# G's updated script, with LC's edits
#python lrn_risk.py --gtdb gtdb_example/gtdbtk.bac120.summary.tsv --virulence_factors_file pima_example/pima_output/vfdb_core.fasta.tabular --amr_determinants_file pima_example/pima_output/AMR_pima_md_2023_02_02.fasta.tabular --blacklist_file lrnrisk_prototype_blacklist.tsv --vf_distribution_file vfdb_dist_beta_80p80q.tsv --amr_distribution_file pima_dist_beta_80p80q.tsv --vfdb_output_file test_vfdb --amr_output_file test_amr --blacklist_output_file test_blacklist
# Test GTDB results from non-Bacilli
python lrn_risk.py --gtdb gtdb_example/gtdbtk.bac120.summary.francisella.tsv --virulence_factors_file pima_example/pima_output/vfdb_core.fasta.francicella.tabular --amr_determinants_file pima_example/pima_output/AMR_pima_md_2023_02_02.fasta.tabular --blacklist_file lrnrisk_prototype_blacklist.tsv --vf_distribution_file vfdb_dist_beta_80p80q.tsv --amr_distribution_file pima_dist_beta_80p80q.tsv --vfdb_output_file test_vfdb --amr_output_file test_amr --blacklist_output_file test_blacklist
