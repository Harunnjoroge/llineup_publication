import pickle
import allel
import pandas as pd
import numpy as np
import itertools
import random
from scipy import stats
import os
from re import *
#import seaborn as sns

# Write a function to calculate burrows' delta from a pair of loci
def burrows_delta(calls, allele_i, allele_j, unbiased = True):
	allele_num = calls.shape[1] * 2
	# This only works if the alleles can only be 0 or 1
	p = ((np.sum(calls[0]) * allele_i) + (allele_num - np.sum(calls[0])) * (1 - allele_i)) / allele_num
	q = ((np.sum(calls[1]) * allele_j) + (allele_num - np.sum(calls[1])) * (1 - allele_j)) / allele_num
	hom_i = 2*allele_i
	hom_j = 2*allele_j
	het = 1
	double_homozygotes = np.sum(np.logical_and(calls[0] == hom_i, calls[1] == hom_j))
	hom_i__het_j = np.sum(np.logical_and(calls[0] == hom_i, calls[1] == het))
	het_i__hom_j = np.sum(np.logical_and(calls[0] == het, calls[1] == hom_j))
	het__het = np.sum(np.logical_and(calls[0] == het, calls[1] == het))
	delta = (2*double_homozygotes +
	         hom_i__het_j +
	         het_i__hom_j +
	         het__het / 2
	        ) / calls.shape[1] - 2*p*q
	if unbiased:
		delta = delta*calls.shape[1] / (calls.shape[1] - 1)
	return(delta)

# Write a function to calculate r_hat_delta:
def Weir_r_hat_delta(calls, allele_i, allele_j):
	allele_num = calls.shape[1] * 2
	# This only works if the alleles can only be 0 or 1
	p = ((np.sum(calls[0]) * allele_i) + (allele_num - np.sum(calls[0])) * (1 - allele_i)) / allele_num
	q = ((np.sum(calls[1]) * allele_j) + (allele_num - np.sum(calls[1])) * (1 - allele_j)) / allele_num
	hom_i = 2*allele_i
	hom_j = 2*allele_j
	h_i = np.sum(calls[0] == hom_i) / calls.shape[1]
	h_j = np.sum(calls[1] == hom_j) / calls.shape[1]
	delta_hat = burrows_delta(calls, allele_i, allele_j)
	r_hat_delta = delta_hat / np.sqrt( (p*(1-p) + (h_i - p**2)) * (q*(1-q) + (h_j - q**2)) )
	return(r_hat_delta)

# For a given pair of loci, calculate the mean of the r values for all pairs of alleles
def mean_r_squared_delta(calls):
	r_00 = Weir_r_hat_delta(calls, 0, 0) ** 2
	r_01 = Weir_r_hat_delta(calls, 0, 1) ** 2
	r_10 = Weir_r_hat_delta(calls, 1, 0) ** 2
	r_11 = Weir_r_hat_delta(calls, 1, 1) ** 2
	return(np.mean([r_00, r_01, r_10, r_11]))

# Calculate the expected r^2 due only to sampling (see Table 1 in Weir 1979), based on sample size S
def expected_r_squared(S):
	if S < 30:
		E = 0.0018 + 0.907/S + 4.44/(S**2)
	else:
		E = 1/S + 3.19/(S**2)
	return(E)

# For a set of genotype calls, calculate the mean r_squared, adjusted to the expectation. In case there are
# more locus permutations than are computationally feasible, we can set the number of allele pairs to randomly
# choose (N_subset). If N_subset = None, use all permutations
def overall_r_squared(call_set, N_subset = 1000, mindist = 1, positions = None):
	# If we want the minimum distance between markers to be greater than 1, then we need to know the positions of the
	# markers
	if mindist > 1 and positions is None:
		raise Exception('If mindist > 1, a list of marker positions is needed.')
	# What is the number of possible permutations?
	possible_permutations = int(call_set.shape[0] * (call_set.shape[0] - 1) / 2)
	# If the number of possible permutations is smaller than N_subset, or if N_subset is None, use all permutations
	if N_subset is None or N_subset >= possible_permutations:
		# Create an iterator for the permutations
		permutations = itertools.combinations(range(call_set.shape[0]), r = 2)
		# Set up the vector of r values
		r = np.array([])
		for p in permutations:
			if mindist > 1:
				this_distance = positions[p[1]] - positions[p[0]]
				if this_distance < mindist:
					continue
			these_calls = call_set[p, :]
			r = np.append(r, mean_r_squared_delta(these_calls))
	# Otherwise, use a random subset of permutations
	else :
		# the code below calculated all permutations and then chose a random subset of them. I think that for the
		# number of locus pairs that we are thinking of using, it will be a lot quicker to just sample them randomly
		# and kick out any repeats
#		# Create an iterator for the permutations
#		permutations = itertools.combinations(range(call_set.shape[0]), r = 2)
#		# Choose which permutations we will use
#		select_permutations = np.sort(random.sample(range(possible_permutations), N_subset))
#		# The islice function will skip a number of steps ahead in the iterator, so we need to change this array of
#		# choices into the number of steps to skip from the previous one
#		select_steps = select_permutations - np.concatenate([[0], select_permutations[:-1] + 1])
#		# Set up the vector of r values
#		r = np.array([])
#		for i, s in enumerate(select_steps):
#			p = next(itertools.islice(permutations, s, s+1))
#			if mindist > 1:
#				this_distance = positions[p[1] - positions[p[0]]]
#				if this_distance < mindist:
#					continue
#			these_calls = call_set[p, :]
#			r = np.append(r, mean_r_squared_delta(these_calls))
		# We choose pairs of loci at random. I'm going to be lazy here and make it so we never use the same locus
		# twice. This will cause problems if N_subset is ever larger than the call_set / 2, but I don't think that
		# this is realistically going to happen.
		gamete1 = np.array(random.sample(range(call_set.shape[0]), N_subset))
		unused = np.setdiff1d(np.array(range(call_set.shape[0])), gamete1)
		gamete2 = np.random.choice(unused, N_subset, replace = False)
		permutations = zip(gamete1, gamete2)
		r = np.array([])
		for p in permutations:
			if mindist > 1:
				this_distance = positions[p[1] - positions[p[0]]]
				if this_distance < mindist:
					continue
			these_calls = call_set[p, :]
			r = np.append(r, mean_r_squared_delta(these_calls))
	return(r)

# Write a function to calculate r_squared using variants from two different chromosomes.
# N_subset is a bit misleading here. If you set it to None, the function will take all pairwise combinations of SNPs
# between the two chromosomes. But if you set N_subset to a value, it will just take that many pairs. So setting
# N_subset to the number of SNPs is not the same as setting it to None.
def twochrom_r_squared(call_set1, call_set2, N_subset=None):
	if N_subset is None or N_subset > min([call_set1.shape[0], call_set2.shape[0]]):
		# Create an iterator for the permutations
		permutations = itertools.product(range(call_set1.shape[0]), range(call_set2.shape[0]))
		print('Using all pairwise combinations.')
	else:
		gamete1 = np.array(random.sample(range(call_set1.shape[0]), N_subset))
		gamete2 = np.array(random.sample(range(call_set2.shape[0]), N_subset))
		permutations = zip(gamete1, gamete2)
		print('Using a subset of loci with one pair per selected locus.')
	# Set up the vector of r values
	r = np.array([])
	for p in permutations:
		these_calls = np.vstack([call_set1[p[0], :], call_set2[p[1], :]])
		r = np.append(r, mean_r_squared_delta(these_calls))
	return(r)

# Write a function to calculate effective population size from the r_squared values and sample size
def Weir_Ne(r_squared_vector, S, c = None):
	expected = expected_r_squared(S)
	observed = np.mean(r_squared_vector)
	adjusted = observed - expected
	# If no recombination rate was given, use the equations from Waples & Do 2008
	if c is None:
		if S < 30:
			Ne = (0.308 + np.sqrt(0.308**2 - 2.08 * adjusted)) / (2*adjusted)
		else:
			Ne = ((1/3) + np.sqrt((1/9) - 2.76 * adjusted)) / (2*adjusted)
	# Otherwise, use the equations from Appendix 1 of Hollenbeck et al 2016
	else:
		gamma = ((1-c)**2 + c**2)/(2*c*(2-c))
		Ne = gamma / adjusted
	if Ne < 0:
		return np.Inf
	else:
		return(Ne)

# Write a function that will take genotypes and calculate the r squared between the two chromosomes
def process_genotypes(genotypes, callsets, chromsize, sample_size = None, all_pairwise_thresh = 0):
	genotypes_1 = []
	genotypes_2 = []
	rs = []
	if sample_size is None:
		sample_size = genotypes[0].shape[1]
		for i, g in enumerate(genotypes):
			print(str(i))
			# remove singletons
			singleton_filter = np.logical_and(np.sum(g, 1) > 1, np.sum(g, 1) < (2*sample_size - 1))
			g = g[singleton_filter, :]
			pos = callsets[i]['variants/POS'][singleton_filter]
			# Let's split the table into two, one for each scaffold
			scaff1 = pos < chromsize
			g1 = g[scaff1, :]
			g2 = g[np.invert(scaff1), :]
			genotypes_1 += [g1]
			genotypes_2 += [g2]
			num_loci = min([g1.shape[0], g2.shape[0]])
			if num_loci < all_pairwise_thresh:
				rs += [twochrom_r_squared(g1, g2)]
			else:
				rs += [twochrom_r_squared(g1, g2, num_loci)]
	else:
		for i, gg in enumerate(genotypes):
			print(str(i))
			g = gg[:, :sample_size]
			# remove singletons and monomorphic loci
			singleton_filter = np.logical_and(np.sum(g, 1) > 1, np.sum(g, 1) < (2*sample_size - 1))
			g = g[singleton_filter, :]
			pos = callsets[i]['variants/POS'][singleton_filter]
			# Let's split the table into two, one for each scaffold
			scaff1 = pos < chromsize
			g1 = g[scaff1, :]
			g2 = g[np.invert(scaff1), :]
			genotypes_1 += [g1]
			genotypes_2 += [g2]
			num_loci = min([g1.shape[0], g2.shape[0]])
			if num_loci < all_pairwise_thresh:
				rs += [twochrom_r_squared(g1, g2)]
			else:
				rs += [twochrom_r_squared(g1, g2, num_loci)]
	return(genotypes_1, genotypes_2, rs)

# Write a function to perform a mann-whitney U test on rs values and report the overall findings
def MannWhit(rs1, rs2, verbose = True, sc1 = 'scenario 1', sc2 = 'scenario 2', pop_size = 'N'):
	u_tests = np.array([stats.mannwhitneyu(x, y).pvalue for (x, y) in zip(rs1, rs2)])
	if verbose:
		if sc2 == None:
			print('A significant difference in r_squared between different simulations with a starting population of ' + str(pop_size) + ' and a ' + sc1 + ' scenario was found in ' + str(np.sum(u_tests < 0.05)) + ' out of ' + str(len(u_tests)) + ' simulations.')
		else:
			print('A significant difference in r_squared between ' + sc1 + ' and ' + sc2 + ' scenarios with a starting population of ' + str(pop_size) + ' was found in ' + str(np.sum(u_tests < 0.05)) + ' out of ' + str(len(u_tests)) + ' simulations.')
	return(u_tests)

# Load the files
allfiles = os.listdir()
vcf_files = list(filter(lambda x: '.vcf' in x, allfiles))
popsizes = [100,1000,10000]
for popsize in popsizes:
	print(f'\nProcessing files for popsize {popsize}')
	crash_files = list(filter(lambda x: f'_crash25_N{popsize}-' in x, vcf_files))
	nocrash_files = list(filter(lambda x: f'_nocrash_N{popsize}_' in x, vcf_files))

	print('\tLoading vcf')
	callsets_crash = [allel.read_vcf(f) for f in crash_files]
	callsets_nocrash = [allel.read_vcf(f) for f in nocrash_files]

	print('\tConverting call data')
	genotypes_crash = [np.sum(c['calldata/GT'], 2) for c in callsets_crash]
	genotypes_nocrash = [np.sum(c['calldata/GT'], 2) for c in callsets_nocrash]

	print('\tProcessing genotypes')
	genotypes_crash_1, genotypes_crash_2, rs_crash = process_genotypes(genotypes_crash, callsets_crash, 50000000)
	genotypes_nocrash_1, genotypes_nocrash_2, rs_nocrash = process_genotypes(genotypes_nocrash, callsets_nocrash, 50000000)

	# Save the results
	print('\tPickling results')
	with open(f'crash_N{popsize}.pickle', 'wb') as f :
		pickle.dump((genotypes_crash_1, genotypes_crash_2, rs_crash), f)
	with open(f'nocrash_N{popsize}.pickle', 'wb') as f :
		pickle.dump((genotypes_nocrash_1, genotypes_nocrash_2, rs_nocrash), f)

	print('\tRunning Mann Whitney tests')
	u_tests = MannWhit(rs_crash, rs_nocrash, sc1 = 'nocrash', sc2 = '10x crash', pop_size = popsize)

	# What is the variation in actual population size estimate?
	sample_size = genotypes_crash[0].shape[1]
	W_crash = np.array([Weir_Ne(x, sample_size, c = 0.5) for x in rs_crash])
	W_nocrash = np.array([Weir_Ne(x, sample_size, c = 0.5) for x in rs_nocrash])
	W_values = pd.DataFrame({f'W_crash_n{sample_size}': W_crash, f'W_nocrash_n{sample_size}': W_nocrash})

	# Smaller sample size
	smaller_sample_sizes = [x for x in [750,500,200,100,75,50,20] if x < sample_size]
	for ss in smaller_sample_sizes:
		print(f'\n\tTrying out a sample size of {ss}:')
		genotypes_crash_1_ss, genotypes_crash_2_ss, rs_crash_ss = process_genotypes(genotypes_crash, callsets_crash, 50000000, ss)
		genotypes_nocrash_1_ss, genotypes_nocrash_2_ss, rs_nocrash_ss = process_genotypes(genotypes_nocrash, callsets_nocrash, 50000000, ss)

		u_tests_ss = MannWhit(rs_crash_ss, rs_nocrash_ss, sc1 = 'nocrash', sc2 = '10x crash', pop_size = popsize)

		# What is the variation in actual population size estimate?
		W_crash_ss = np.array([Weir_Ne(x, ss, c = 0.5) for x in rs_crash_ss])
		W_nocrash_ss = np.array([Weir_Ne(x, ss, c = 0.5) for x in rs_nocrash_ss])
		W_values[f'W_crash_n{ss}'] = W_crash_ss
		W_values[f'W_nocrash_n{ss}'] = W_nocrash_ss

		with open(f'crash_N{popsize}_n{ss}.pickle', 'wb') as f :
			pickle.dump((genotypes_crash_1_ss, genotypes_crash_2_ss, rs_crash_ss), f)
		with open(f'nocrash_N{popsize}_n{ss}.pickle', 'wb') as f :
			pickle.dump((genotypes_nocrash_1_ss, genotypes_nocrash_2_ss, rs_nocrash_ss), f)

	W_values.to_csv(f'W_values_N{popsize}.csv', sep = '\t')
