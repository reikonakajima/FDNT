* all.train: training data based on all FDNTRing runs up to Oct 28 2013
* net.out:   3.10.10.1 network trained for multiplicative bias on this
* test.sm:   sm macro to plot data and model at arbitrary parameter combinations
* with _c:   same thing, but for additive bias (parametrized as c=<e2>/psf_e2)
* make_annz_input[_c].py:  read a directory of Reiko's .dat files and turn them into lines that ANNz can read (with 1/SN^2, e, exp(-r), corresponding bias & uncertainty on that bias in it); should run this on all available directories and paste output into one big file (this is how to make all.train)


Train with something like ../src/annz_train arch.3.10.10.1.net all_c.train all_c.train net_c.out 1

Call Supermongo and run

macro read test.sm
varrad_c net_c.out 30 0.3
