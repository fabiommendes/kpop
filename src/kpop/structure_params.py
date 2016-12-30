MAINPARAMS_DEFAULTS = dict(
    #
    # Basic Program Parameters
    #

    # (int) number of populations assumed
    maxpops=2,

    # (int) length of burnin period
    burnin=10000,

    # (int) number of MCMC reps after burnin
    numreps=20000,


    #
    # Input/Output files
    #

    # (str) name of input data file
    infile='<infile>',

    # (str) name of output data file
    outfile='/dev/stdout',


    #
    # Data file format
    #

    # (int) number of diploid individuals in data file
    numinds='<N:int>',

    # (int) number of loci in data file
    numloci='<J:int>',

    # (int) ploidy of data
    ploidy=2,

    # (int) value given to missing genotype data
    missing=-9,

    # (B) store data for individuals in a single line
    onerowperind=1,

    # (B) Input file contains individual labels
    label=1,

    # (B) Input file contains a population identifier
    popdata=0,

    # (B) Input file contains a flag which says whether to use popinfo when
    # USEPOPINFO==1
    popflag=0,

    # (B) Input file contains a location identifier
    locdata=0,

    # (B) Input file contains phenotype information
    phenotype=0,

    # (int) Number of additional columns of data before the genotype data start.
    extracols=0,

    # (B) data file contains row of marker names
    markernames=0,

    # (B) data file contains dominant markers (eg AFLPs) and a row to indicate
    # which alleles are recessive
    recessivealleles=0,

    # (B) data file contains row of map distances between loci
    mapdistances=0,


    #
    # Advanced data file options
    #

    # (B) Data are in correct phase (relevant for linkage model only)
    phased=0,

    # (B) the data for each individual contains a line indicating phase
    # (linkage model)
    phaseinfo=0,

    # (B) the phase info follows a Markov model.
    markovphase=0,

    # (int) for use in some analyses of polyploid data
    notambiguous=-999
)


EXTRAPARAMS_DEFAULTS = dict(
    {'lambda': 1.0},
    noadmix=0,
    linkage=0,
    usepopinfo=0,
    locprior=0,
    freqscorr=1,
    onefst=0,
    inferalpha=1,
    popalphas=0,
    alpha=1.0,
    inferlambda=0,
    popspecificlambda=0,

    # priors
    fpriormean=0.01,
    fpriorsd=0.05,
    unifprioralpha=1,
    alphamax=10.0,
    alphapriora=1.0,
    alphapriorb=2.0,
    log10rmin=-4.0,
    log10rmax=1.0,
    log10rpropsd=0.1,
    log10rstart=-2.0,

    # using prior population info (usepopinfo)
    gensback=2,
    migrprior=0.01,
    pfrompopflagonly=0,

    # locprior model for using location information
    locispop=0,
    locpriorinit=1.0,
    maxlocprior=20.0,

    # output options
    printnet=1,
    printlambda=1,
    printqsum=1,
    sitebysite=0,
    printqhat=0,
    updatefreq=100,
    printlikes=0,
    intermedsave=0,
    echodata=1,

    # next 3 are for collecting distribution of q:
    ancestdist=0,
    numboxes=1000,
    ancestpint=0.90,

    # miscellaneous
    computeprob=1,
    admburnin=500,
    alphapropsd=0.025,
    startatpopinfo=0,
    randomize=1,
    seed=2245,
    metrofreq=10,
    reporthitrate=0,
)
