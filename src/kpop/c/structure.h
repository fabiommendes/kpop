
/*Mechanics of the program, should not need to be changed. */
#define UNDERFLO   1e-100 /*DBL_MIN*1.1*/	/*Get worried about numbers that are smaller than this */
#define STOREFREQ  1		/*frequency of recording the likelihood */
#define UNASSIGNED -9		/*missing Data are placed in this population. Should be <0 */
#define STRLEN     200		/* max length of string used to store data file name */
#define GENELEN 15              /* max length of gene name */
#define LABELLEN   12		/* max length of string used to store individ label */
#define VERBOSE    0		/* (B) Print run-time details to screen */
#define ALLELLEN   10		/* integer data in data file are read in as strings 
				   of this length */
#define MAXINDS    3		/* print 2*MAXINDS-1 individuals when printing data to screen */
#define BANNERFREQ 10		/*freq of printing banner for update stuff */
#define LAMBDAMAX  10           /*max value of lambda*/
#define LAMBDAPROPSD 0.3          /*SD for proposal for lambda*/

/*return appropriate positions in each array... */
#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))	/*Geno */
#define TransPos(loc,pos) ((MAXALLELES)*(loc)+(pos))	/*Translation */
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))	/* Z */
#define SiteBySiteSumPos(ind,line,loc,pop) ((ind)*(LINES)*(NUMLOCI)*(MAXPOPS)+(line)*(NUMLOCI)*(MAXPOPS)+(loc)*(MAXPOPS)+(pop))
	/* P */
#define DiploidSiteBySiteSumPos(ind,pop2,loc,pop) ((ind)*(MAXPOPS)*(NUMLOCI)*(MAXPOPS)+(pop2)*(NUMLOCI)*(MAXPOPS)+(loc)*(MAXPOPS)+(pop))
	
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))	/* Q */
#define RTransitProbPos(loc,line,pop) ((loc)*MAXPOPS*LINES+(line)*MAXPOPS+(pop))
#define DiploidRTransitProbPos(loc,pop,pop2) ((loc)*MAXPOPS*MAXPOPS+(pop)*MAXPOPS+(pop2))
#define EpsPos(loc,allele) ((loc)*(MAXALLELES)+(allele))	/* Epsilon */
#define NumAFromPopPos(pop,allele) ((pop)*(MAXALLELES)+(allele))	/* NumAFromPop */
	/*AncestDist */
#define AncestDistPos(ind,pop,box) ((ind)*(MAXPOPS)*(NUMBOXES)+(pop)*(NUMBOXES)+(box))
#define PostProbsPos(pop,gen) ((pop)*(GENSBACK+1)+(gen))	/*PostProbs-Q update with pop */
	/*UsePopProbs */
#define UsePPrPos(ind,pop,gen) ((ind)*(MAXPOPS)*(GENSBACK+1)+(pop)*(GENSBACK+1)+(gen))
#define PhasePos(ind,loc)   ((ind)*(NUMLOCI)+(loc))
#define SquarePos(pop,pop2) ((pop)*(MAXPOPS)+(pop2))
#define MarkernamePos(loc,pos) ((loc)*(GENELEN)+pos)

#define ClustSizePos(pop) (NUMLOCATIONS + (pop))

#define AlphaPos(loc, pop) ((loc)==NUMLOCATIONS ? (pop) : (MAXPOPS*((loc)+1)+(pop)))

#define LocPriorPos(loc, pop) (1 + (((loc)==NUMLOCATIONS) ? (pop) : (MAXPOPS*((loc)+1)+(pop))))


/*NumAlleles[loc] */
/*Correls[loc] */

/*GLOBAL VARIABLES WHICH ARE SET IN THE mainparams AND extraparams FILES;
   SEE EXPLANATIONS THERE */

struct IND {				/*Non-genotype info for each individual */
  char Label[LABELLEN];
  int Population;
  int PopFlag;
  int Phenotype;
  int Location;
  int myloc;
};

/*when adding new parameters, need to make changes in 4 places: 
   (1) here, (2) extraparams, (3) SetValue (params.c), (4) PreSetValues
   (params.c) */


/*Data File */
char DATAFILE[STRLEN + 1];
char OUTFILE[STRLEN + 1];
int NUMINDS, NUMLOCI, MISSING, LABEL, POPDATA, LINES;
int POPFLAG, PHENOTYPE, EXTRACOLS;
int ONEROWPERIND; 
int RECESSIVEALLELES;
int NOTAMBIGUOUS;
/*Program Parameters */
int MAXPOPS, BURNIN, NUMREPS;
/*Program options */
int USEPOPINFO, INFERALPHA, INFERLAMBDA,POPSPECIFICLAMBDA, PFROMPOPFLAGONLY;
int POPALPHAS, COMPUTEPROB, NOADMIX, ADMBURNIN;
int MAPDISTANCES, MARKERNAMES, LINKAGE, PHASED,PHASEINFO,MARKOVPHASE;
		   
/*Output Options */
int UPDATEFREQ, PRINTLIKES, INTERMEDSAVE, PRINTKLD, PRINTNET, PRINTLAMBDA;
int ECHODATA, ANCESTDIST, NUMBOXES, GENSBACK;
double ANCESTPINT;
int PRINTQHAT;
int PRINTQSUM;
int SITEBYSITE;
/*Priors */
int FREQSCORR, UNIFPRIORALPHA, ONEFST;
double MIGRPRIOR, ALPHA, FPRIORMEAN,FPRIORSD, LAMBDA, MIGRPRIOR;
double ALPHAMAX, ALPHAPRIORA, ALPHAPRIORB;
/*Miscellaneous */
int STARTATPOPINFO;
double ALPHAPROPSD, LOG10RPROPSD,LOG10RMIN,LOG10RMAX,LOG10RSTART;

/*Relatively obscure options */
int RANDOMIZE, METROFREQ, REPORTHITRATE, STARTATPOPINFO;

/*Parameters defined in STRATparams */
int NUMSIMSTATS, NUMPHENS, POOLFREQ, LOCUSxONLY, MISSINGPHENO, PHENOTYPECOL;
double EMERROR;

/*LocPrior parameters (melissa added 7/12/07)*/
int LOCPRIOR, NUMLOCATIONS;
int UPDATELOCPRIOR;
double LOCPRIORINIT, MAXLOCPRIOR;
int LOCDATA, LOCISPOP;
double LOCPRIORSTEP;
int PRINTLOCPRIOR;

/* seed parameter, melissa added 1/8/08 */
int SEED;

/* options hidden by user but still used by program*/
int NOQS, POSTERIOR,PICTUREFILE,INDIVIDUALR;




/*Needed for program (not input) */
int MAXALLELES;			/*maximum number of alleles at any locus */
int ORIGMISSING;		/*value of missing alleles may be changed when the
				   data are translated; this is the original value */
int ORIGMISSINGPHENO;           /*ditto for phenotypes in STRAT*/
int NOALPHA;			/*Alpha not used because prior pop model is used instead */

/*function declarations------ */
void Kill ();			/*abort program */

void Welcome (FILE * file);
int FindMaxAlleles (int *Geno, int *Recessive);

