from kpop.all import *

# Create 3 ancestral populations
num_ind = 100
num_loci = 104
popA = Population.random(num_ind, num_loci, id='A')
popB = Population.random(num_ind, num_loci, id='B')
popC = Population.random(num_ind, num_loci, id='C')

# Create a parental multi-population and some mixtures
parental = popA + popB + popC
popAB = parental.new_admixed_population([0.5, 0.5, 0.0], 50, id='AB')
popAAB = parental.new_admixed_population([0.75, 0.25, 0.0], 50, id='AAB')
popABC = parental.new_admixed_population([0.25, 0.25, 0.5], 50, id='ABC')
child = popAB + popAAB + popABC

# Re-calculate frequencies for parental populations
popA.freqs = popA.empirical_freqs()
popB.freqs = popB.empirical_freqs()
popC.freqs = popC.empirical_freqs()

# Infer admixture coefficients from children using different classifiers
t0 = time()
child_maxlike = parental.admixture_classify(child, classifier='maxlike')
dt1 = time() - t0

t0 = time()
child_admix = parental.admixture_classify(child, classifier='admixture')
dt2 = time() - t0


# Now we compute the total accumulated error
def err1(pop):
    mixtures = [[0.5, 0.5, 0.0], [0.75, 0.25, 0.0], [0.25, 0.25, 0.5]]
    mixtures = np.array(mixtures)
    errors = []
    for sub, orig, prob in zip(pop.populations, child.populations, mixtures):
        cum_error = 0
        for ind in sub:
            cum_error += abs(prob - ind.admixture_q.encode([0, 1, 2])).sum()
        errors.append(cum_error)

    return np.array(errors)

err1 = err1(child_maxlike)
err2 = err1(child_admix)
print(err1, err1.sum(), dt1)
print(err2, err2.sum(), dt2)

# Plot _results
child_maxlike.plot.admixture()
child_admix.plot.admixture()
