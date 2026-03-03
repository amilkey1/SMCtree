# SMCtree
SMC for divergence-time estimation.

General settings to be specified in smctree.conf file:\
`rnseed` = pseudorandom number seed (ex. 123)\
`nthreads` = number of threads if using multithreading (default value = 1, ex. 1)\
`startmode`= simulate or run SMC (options are sim or smc, ex. smc)\
`subset` = a string defining a parittion subset (ex. first:1-1000)\
`nparticles` = number of particles (ex. 10000)\
`filename` = name of file containing data (ex. sim.nex)\
`ngroups` = number of subgroups, total number of particles = nparticles * ngroups (default is 1.0, ex. 5)\
`save_every` = save one out of every n trees (default is 1.0, ex. 100)\
`run_on_empty` = run without likelihood (default is false, ex. false)

Parameter settings:\
`lambda` = per lineage speciation rate assumed for the Yule model, mean lambda if estimating lambda (ex. 5.0)\
`root_age` = root age for birth death model (ex. 1.0)\
`root_age_min` = min root age for birth death model if estimated (ex. 1.0)\
`root_age_max` = max root age for birth death model if estimated (ex. 2.0)\
`mu` = per lineage extinction rate assumed for the birth death model, mean mu if estimating mu (ex. 0.3)\
`clock_rate` = clock rate if fixed, mean clock rate if estimating clock rate (ex. 1.0)\
`savememory` = save memory by recalculating all partials each round (default is false, ex. true)\
`model` = model to use for likelihood calculations (default is JC, ex. HKY, current implementation only allows JC or HKY)\
`verbosity` = level of output (default is 1, ex. 2)\
`kappa`= value of kappa transition rate ratio divided by transversion rate ratio (default is 1.0, ex. 5.0)\
`base_frequencies` = base frequencies (default is all 0.25, ex. 0.3, 0.3, 0.2, 0.2)\
`plus_G` = use gamma rate heterogeneity (default is false, ex. true)\
`gamma_rate_var` = rate variance for +G model (default is 1.0, only used when plus_G = true, ex. 1.0)\
`plus_I` = use proportion of invariable sites rate heterogeneity model (default is false, ex. true)\
`pinvar` = proportion of invariable sites (default is 0.0, only used when plus_I = true, ex. 0.8)\
`relative_rates` = relative rates by locus (default is all equal, ex. 1.8, 2.2)\
`estimate_lambda` = estimate birth rate (default is false, ex. true)\
`estimate_mu` = estimate death rate (default is false, ex. false)\
`estimate_root_age` = estimate root age (default is false, ex. true)\
`estimate_clock_age` = estimate clock rate (default is false, ex. true)

Fossil settings:\
`fossil` = a string defining a fossil (ex. Ursus_abstrusus 1.8–5.3 4.3 (4.3 is time, 1.8-5.3 is prior range)\
`taxset` = a string defining a taxon set (ex. Ursinae: Helarctos_malayanus Melursus_ursinus Ursus_abstrusus Ursus_americanus Ursus_arctos Ursus_maritimus Ursus_spelaeus Ursus_thibetanus)

Simulation settings:\
`simfnprefix` = prefix of files in which to save simulated trees and data if startmode is 'sim' (ex. sim)\
`simntaxa` = number of taxa to simulate if startmode is 'sim' (ex. 12)\
`simlambda` = true speciation rate for simulating tree under the Yule model if startmode is 'sim' (ex. 5.0)\
`simmu` = true extinction rate for simulating tree under the constant-rates Birth-Death model if startmode is 'sim' (ex. 0.1)\
`simrho` = true extant taxon sampling rate for simulating tree under the constant-rates Birth-Death model if startmode is 'sim' (ex. 1.0, note SMC cannot currently sample this parameter)\
`simrootage` = true root age for simulating tree under the constant-rates Birth-Death model if startmode is 'sim' (ex. 1.0)\
`simclockrate` = true clock rate for simulating tree under constant-rates Birth-Death model if startmode is 'sim' (ex. 1.0)\
`simfossil` = a string defining a fossil for simulations only (ex. Ursus_abstrusus 1.8–1.8 1.8 (1.8 is minimum clade depth))

Validation settings - do not use for analyses:\
`ruv`  = run ruv analysis (default false, ex. true)\
`coverage` = run coverage analysis (default false, ex. true)\
`sim_dir` = location of directory with simulation results (default . , ex. ../sim)

