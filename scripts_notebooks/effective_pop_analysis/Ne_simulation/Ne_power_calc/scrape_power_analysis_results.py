# I was silly and wrote the Ne_estimation...py files such that they didn't save
# the results of the power of the r-squared comparisons to a convenient table. So
# this script is needed to really clumsily scrape the log files and pull out the
# results.

import re
import pandas as pd

def scrape(fn, max_sample_size):
    # I told you this code was going to be clumsy
    results = {'N': [],
               'n': [],
               'successes': [],
               'trials': []}
    with open(fn) as f:
        for line in f:
            # It just gets worse and worse doesn't it
            if re.search('Processing files for popsize', line):
                sample_size = max_sample_size
            if re.search('^A significant difference in r_squared', line):
                N = re.findall('\\d+',
                               re.findall('starting population of \\d+', line)[0]
                              )[0]
                successes = re.findall('\\d+',
                                       re.findall('was found in \\d+', line)[0]
                                       )[0]
                trials = re.findall('\\d+',
                                      re.findall('out of \\d+ simulations', line)[0]
                                     )[0]
                results['N'] += [N]
                results['n'] += [sample_size]
                results['successes'] += [successes]
                results['trials'] += [trials]
            if re.search('Trying out a sample size', line):
                sample_size = re.findall('\\d+', line)[0]
    output_df = pd.DataFrame(results)
    return(output_df)

# For the wf simulations, we simulated 250 samples (500 haplotypes)
# for the coalescent simulations, we simulated 500 samples (1000 haplotypes)
wf_simulations = scrape('wf_simulations/Ne_estimation.log', 250)
coal_main_simulations = scrape('coalescent_simulations/Ne_estimation.log', 500)
coal_N10000_simulations = scrape('coalescent_simulations/Ne_estimation_continued_10000.log', 500)
coal_N100000_simulations = scrape('coalescent_simulations/Ne_estimation_continued_100000.log', 500)
coal_simulations = pd.concat([coal_main_simulations,
                              coal_N10000_simulations,
                              coal_N100000_simulations])

wf_simulations.insert(loc = 0, column = 'simulation_type', value = 'wf')
coal_simulations.insert(loc = 0, column = 'simulation_type', value = 'coalescent')

all_simulations = pd.concat([wf_simulations, coal_simulations])

all_simulations.to_csv('power_calculations.csv', sep = '\t', index = False)
