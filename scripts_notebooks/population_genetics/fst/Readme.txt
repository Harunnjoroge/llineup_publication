# In this analysis, we load genotype data that is downloaded from Malariagen using the [./fst/extract whole genome genotype.ipynb ]
# We load a csv file containing permutated control_phase or Location data, the output of the [./fst/llineup_random_permuation.r]
# Fst between the different cohorts are calculated and saved as rds. These are visualised using the [./fst/calculate_Fst_window_p_value.r]