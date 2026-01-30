from sys import stdout
import numpy as np
import msprime

# Recombination rate
rho = 1e-8
rates = [rho, 0.5, rho, 0]

# Mutation rate
mu = 1e-8

# Generations since crash
crash_t = 25
 
# Sample size (number of haploptypes, so number of individuals will be half of this)
Nh = 10000

# Chromosome sizes
chrom_size = 50000000
positions = [0, chrom_size-1, chrom_size, 2*chrom_size]
recomb = msprime.RecombinationMap(positions = positions, rates = rates, num_loci = chrom_size*2)

# Number of replicates
num_rep = 5

def simulate_and_save(currentN, recomb, demog, n_rep, note = '', should_save = True):
	if demog is None:
		print('\nSimulating a population of size ' + str(currentN) + ' with two scaffolds of size ' + str(int(recomb.get_length()/2)) + '.')
		tree_sequences = msprime.simulate(Nh, Ne = currentN, recombination_map = recomb, num_replicates = n_rep)
		filename_root = 'population_2chrom_g' + str(int(recomb.get_length()/2)) + '_nocrash_N' + str(currentN) + '_r' + str(rho) + '_mu' + str(mu) + note + '_'
	else:
		print('\nSimulating a population of size ' + str(demog[0].initial_size) + ' crashing to ' + str(currentN) + ' with two scaffolds of size ' + str(int(recomb.get_length()/2)) + '.')
		tree_sequences = msprime.simulate(Nh, Ne = currentN, recombination_map = recomb, demographic_events = demog, num_replicates = n_rep)
		filename_root = 'population_2chrom_g' + str(int(recomb.get_length()/2)) + '_crash' + str(crash_t) + '_N' + str(demog[0].initial_size) + '-' + str(currentN) + '_r' + str(rho) + '_mu' + str(mu) + note + '_'
	mts = []
	g = []
	for i, ts in enumerate(tree_sequences):
		tree_number = str(i).zfill(len(str(num_rep-1)))
		print('\tCreating tree number ' + tree_number)
		stdout.flush()
		mutated_tree_sequence = msprime.mutate(ts, rate = mu)
		gen = mutated_tree_sequence.genotype_matrix()
		mts += [mutated_tree_sequence]
		g += [gen]
		if should_save:
			with open(filename_root + tree_number + '.vcf', 'w') as vcf_file:
				# Need the "position transform", because without it I got an error about a SNP ending up having position 0
				mutated_tree_sequence.write_vcf(vcf_file, 2, position_transform = lambda x: 1 + np.round(x))
	return(mts, g)

# Population sizes and demographics
startNs = [100000]
crashNs = [int(round(x/10, 3)) for x in startNs]

for j in range(len(startNs)):
	startN = startNs[j]
	crashN = crashNs[j]
	print('Starting population ' + str(startN))
	stdout.flush()
	demog = [msprime.PopulationParametersChange(time = crash_t, initial_size = startN, growth_rate = 0)]
	dd = msprime.DemographyDebugger(Ne = crashN, demographic_events = demog)
	print('\nDemographic history for crash (' + str(startN) + ' to ' + str(crashN) + ').')
	dd.print_history()

	# Simulate the populations before the crash
	tree_sequences_nocrash = simulate_and_save(startN, recomb, None, num_rep, '_beforecrash')
	# And after the crash
	#tree_sequences_crash = simulate_and_save(crashN, recomb, demog, num_rep, '_crash')


