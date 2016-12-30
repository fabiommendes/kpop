

/*=======================================================

  STRUCTURE.C

  Program for inferring population structure using multilocus
  genotype data.

  Code written by Daniel Falush, Melissa Hubisz, and Jonathan Pritchard

  See additional details in README file.

  =========================================================*/
#define VERSION "2.3.4 (Jul 2012)"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "mymath.h"
#include "ran.h"
#include "params.h"
#include "datain.h"
#include "output.h"

void InitializeZ (int *Geno, struct IND *Individual, int *Z);
void UpdateQAdmixture (double *Q, int *Z, double *Alpha, struct IND *Individual);

/*========================================================
  ==========================================================*/
void
Welcome (FILE * file)
{
  fprintf (file, "\n\n");
  fprintf (file, "----------------------------------------------------\n");
  fprintf (file, "STRUCTURE by Pritchard, Stephens and Donnelly (2000)\n");
  fprintf (file, "     and Falush, Stephens and Pritchard (2003)\n");
  fprintf (file, "       Code by Pritchard, Falush and Hubisz\n");
  fprintf (file, "             Version %s\n", VERSION);
  fprintf (file, "----------------------------------------------------\n");

  fprintf (file, "\n\n");
  fflush(file);
}
/*-------------------------------------------------------*/
void Kill ()                            /*exit program */
{
  printf ("\nExiting the program due to error(s) listed above.\n\n");
  exit (1);
}
/*---------------------------------------*/
void CheckParamCombinations()
{
  if ((LINES>2) && (USEPOPINFO==1))
  {
    printf(" USEPOPINFO does not work for ploidy >2\n");
    Kill();
  }
  if ((PHASEINFO) && (LINES!=2))
  {
    printf("phase info is only applicable to diploid data!! \n");
    Kill();
  }
  if (LINES ==2 && PHASEINFO==1 && PHASED ==0 && MARKOVPHASE==-1 && LINKAGE)
  {
    printf("You need to specify a phase model using the parameter MARKOVPHASE!! \n");
    Kill();
  }
  if (LINKAGE && !MAPDISTANCES)
  {
    printf("Map distance information is required in order to run the linkage model. \n");
    Kill();
  }

  if ((LINKAGE) && (!PHASED) && (LINES!=2))
  {
    printf("unphased data only permitted for diploids!! \n");
    Kill();
  }

  if ((LINKAGE) && (NOADMIX))
  {
    printf("please choose the LINKAGE option or the NOADMIX option, but not both \n");
    Kill();
  }

  if ((INFERLAMBDA) && (FREQSCORR))
  {
    printf("Warning!, choosing both INFERLAMBDA and FREQSCORR parameters may leave the computer with insufficient information to estimate either parameter accurately \n");
  }


  if ((LINKAGE==0) && (SITEBYSITE))
  {
    printf("SITEBYSITE is not currently implemented for no-linkage model\n");
    SITEBYSITE=0;
  }


  if (((NOADMIX) || (LINKAGE)) && ADMBURNIN >= BURNIN)
  {
    printf("The ADMBURNIN should finish before the BURNIN!! \n");
    Kill();
  }

  if ((PFROMPOPFLAGONLY) && (!(POPFLAG)))
  {
    printf("PFROMPOPFLAGONLY can only be turned on when the data file contains POPFLAG data\n");
    Kill();
  }
  if (LINES>2 && RECESSIVEALLELES && MISSING==NOTAMBIGUOUS)
  {
    printf("The code for missingdata (MISSING) should be set differently to the code (NOTAMBIGUOUS) for loci whose genotypes are known");
    Kill();
  }
  if (LOCPRIOR && LINKAGE) {
    printf("LOCPRIOR model is not compatible with linkage model\n");
    Kill();
  }
  if (LOCPRIOR && USEPOPINFO) {
    printf("LOCPRIOR model is not compatible with USEPOPINFO option\n");
    Kill();
  }
  /*  if (RANDOMIZE && SEED!=-1) {
      printf("Warning: Seed from parameter file will not be used as RANDOMIZE is set to 1.  SEED in output file will be random seed drawn from time\n");
      }*/  /* modified by JKP */

  if (RANDOMIZE)
    printf("Note: RANDOMIZE is set to 1. The random number generator will be initialized using the system clock, ignoring any specified value of SEED.\n");

}
/*---------------------------------------*/
/*void FreeAll (int *Geno, double *Mapdistance, char *Markername,
              struct IND *Individual, int *Translation,
              int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP, double *Epsilon,
              double *Fst, int *NumLociPop, double *PSum, double *QSum,
              double *SiteBySiteSum,
              double *FstSum, int *AncestDist, double *UsePopProbs, double *R,
              double *sumR, double *varR, double *LocPrior, double *sumLocPrior)
*/
void FreeAll(double *Mapdistance, double *Phase, int *Phasemodel, double *lambda, double *sumlambda,
             char *Markername, int *Geno, int* PreGeno, int* Recessive, struct IND *Individual,
             int *Translation, int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP,
             double *R, double *sumR, double *varR, double *Epsilon, double *SumEpsilon, double *Fst,
             double *FstSum, int *NumLociPop, double *PSum, double *QSum, double *SiteBySiteSum,
             int *AncestDist, double *UsePopProbs, double *LocPrior, double *sumLocPrior,
             double *Alpha, double *sumAlpha, double *sumIndLikes, double *indLikesNorm)
{
  /** these variables are calloc'd in main and freed in the same order */

  free (Mapdistance);
  
  free(Phase);
  if (LINES==2 && PHASED == 0) {
    free(Phasemodel);
  }
  


  free(lambda);
  free(sumlambda);

  
  free (Markername);
  free (Geno);

  if (RECESSIVEALLELES) {
    free(PreGeno);
    free(Recessive);
  }


  free (Individual);
  free (Translation);
  free (NumAlleles);
  
  free (Z);
  free (Z1);
  
  free (Q);
  free (P);
  
  free (LogP);
  free (R);
  free (sumR);
  free (varR);

  free (Epsilon);
  
  if (FREQSCORR) {
    free(SumEpsilon);
  }
  
  
  free (Fst);
  free (FstSum);
  free (NumLociPop);
  
  free (PSum);
  free (QSum);

  if ( SITEBYSITE) {
    free (SiteBySiteSum);
  }

  if (ANCESTDIST) {
    free (AncestDist);
  }

  if (USEPOPINFO)  {
    free (UsePopProbs);
  }

  if (LOCPRIOR) {
    free(LocPrior);
    free(sumLocPrior);
  }

  
  free(Alpha);
  free(sumAlpha);

  free(sumIndLikes);
  free(indLikesNorm);
  
}
/*---------------------------------------*/
void
PrintTranslation (int *Translation, int *NumAlleles)
{
  int loc, allele;
  printf ("Translation matrix:\n");
  for (loc = 0; loc < NUMLOCI; loc++)
  {
    for (allele = 0; allele < NumAlleles[loc]; allele++)
      printf ("%2d ", Translation[TransPos (loc, allele)]);
    printf ("\n");
  }
  printf ("\n");
  for (loc = 0; loc < NUMLOCI; loc++)
    printf ("%d ", NumAlleles[loc]);
  printf ("\n");


}

/*-------------------------------------------*/
void
InitFromGeogPop (int *Geno, struct IND *Individual, int *Z, int verbose)
{
  /*initialize the population of origin of each allele. These are
    started in their given populations (when there is population info).
    It is assumed that the populations are numbered 1..K.  If the
    population identifier is out of range, the allele is assigned to a
    random population.  These is used to check that the MCMC is not
    producing unexpected results because it couldn't find the mode. */

  int ind, line, loc;
  int poperrors = 0;
  int pop;

  if (!(POPDATA)) {
    InitializeZ (Geno, Individual, Z);
    if (verbose) {
      printf ("Starting from a random configuration because POP=0\n");
    }
  } else {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (line = 0; line < LINES; line++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
          if (Geno[GenPos (ind, line, loc)] == MISSING) {
            Z[ZPos (ind, line, loc)] = UNASSIGNED;
          } else {
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
              poperrors++;
            }
          }
        }
      }
    }
    if (verbose) {
      printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
    }
    if ((verbose) && (poperrors)) {
      printf ("WARNING: unable to initialize %d individuals to the predefined\n", poperrors);
      printf ("populations because their population names were not in the range\n");
      printf ("{1..%d}.  These individuals were initialized at random.\n",MAXPOPS);
    }
  }
}
/*---------------------------------------*/
void
InitializeZ (int *Geno, struct IND *Individual, int *Z)
{
  /*initialize the population of origin of each allele. I pick these at
    random because this seems to produce better behaviour than starting
    out with everybody in one pop, for instance.  I also set missing Data
    to the unassigned pop from the start.  STARTATPOPINFO indicates that
    individuals should be started at their given populations. It is
    assumed that the populations are numbered 1..K.  If the population
    identifier is out of range, the allele is assigned to a random
    population. */

  int ind, line, loc;
  int allele;
  int pop;

  for (ind = 0; ind < NUMINDS; ind++) {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)]; /*missing data */
        if (allele == MISSING) {
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
        } else {  /*------data present-----------*/
          if ((STARTATPOPINFO) && (POPDATA)) {    /*use given pops as initial Z */
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
            }
          } else {          /*initial Z random */
            Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
          }
        }
        /*printf("%d ",Z[ZPos(ind,line,loc)]); */
      }
      /*printf("\n"); */
    }
  }

  if ((STARTATPOPINFO) && (POPDATA)) {
    printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
  }
}
/*---------------------------------------*/
void
InitFreqPriors (double *Epsilon, double *Fst, int *Geno, int *NumAlleles)
{
  int ind, line, loc, allele,pop;
  int value;
  int *Count /*[MAXALLELES] */ ;        /*stores number of copies of each allele
                                          at present locus */
  int total;

  if (!(FREQSCORR)) {           /*allele frequencies uncorrelated across populations */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] = LAMBDA;
      }
    }
    for (pop = 0; pop< MAXPOPS; pop++) {
      Fst[pop] = 0.5;
    }
  } else {                      /*correlated allele frequencies------------------- */
    Count = calloc (MAXALLELES, sizeof (int));
    if (Count == NULL) {
      printf ("Error in assigning memory, InitFreqPriors\n");
      Kill ();
    }

    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Count[allele] = 0;
      }
      total = 0;
      for (ind = 0; ind < NUMINDS; ind++) {
        for (line = 0; line < LINES; line++) {
          value = Geno[GenPos (ind, line, loc)];
          if (value != MISSING) {
            total++;
            if ((value < 0) || (value >= NumAlleles[loc])) {
              printf ("WARNING: expected allele value, InitFreqPriors: loc %d, allele %d\n", loc, value);
            } else {
              Count[value]++;
            }
          }
        }
      }

      /*Start Epsilon at (mean sample frequencies) */
      /* add lambda to all counts to ensure positive frequencies
       * for recessive model etc */
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] =
            (((double) LAMBDA +
              (double) Count[allele]) / ((double) NumAlleles[loc] *
                                         (double) LAMBDA +
                                         (double) total));
        /* ((double) Count[allele] / total); */
      }
      /*printf("\n"); */
    }


    for (pop= 0; pop < MAXPOPS; pop++) {
      Fst[pop] = FPRIORMEAN;    /*start Fst at the prior mean */
    }
    if (Count != NULL) {
      free (Count);
    }
  }                             /*end, correlated allele frequencies------------- */

}
/*---------------------------------------*/


void
CheckPopPriors (struct IND *Individual)
    /*This function is called when USEPOPINFO==1 to check that the
      prior information on populations is ok */
{
  int ind;
  int numnopop;

  if ((MIGRPRIOR < 0.0) || (MIGRPRIOR > 1.0)) {
    printf ("MIGRPRIOR (which is currently set to %1.3f) must be in the range [0.0, 1.0]\n", MIGRPRIOR);
    Kill ();
  }

  if (!(POPDATA)) {
    printf ("Can't apply USEPOPINFO because no POPDATA in input data file\n");
    Kill ();
  }

  if (!(POPFLAG)) {              /*if no popflag, then assume that everybody should be used */
    for (ind = 0; ind < NUMINDS; ind++) {
      Individual[ind].PopFlag = 1;
    }
  }

  /*Check that the given population is within range for all individuals, and if
    not, then turn popflag off for that individual. */
  for (ind = 0; ind < NUMINDS; ind++) {
    if (Individual[ind].PopFlag) {
      if ((Individual[ind].Population < 1) || (Individual[ind].Population > MAXPOPS)) {
        printf ("Warning: population prior for individual %d is %d, which is not\n",
                ind + 1, Individual[ind].Population);
        printf ("  in the range 1..%d.  Population prior for this individual will\n",
                MAXPOPS);
        printf ("  be ignored\n");
        Individual[ind].PopFlag = 0;
      }
    }
  }

  if ((INFERALPHA) && (!(NOADMIX))) {     /*check whether alpha is needed at all */
    numnopop = 0;
    for (ind = 0; ind < NUMINDS; ind++) {
      if (Individual[ind].PopFlag == 0) {
        numnopop++;
      }
    }
    if (numnopop == 0) {
      NOALPHA = 1;
      INFERALPHA = 0;
    }
  }
}


/*GetNumLocations: Melissa added 7/12/07.  Sets the variable NUMLOCATIONS and also setse
  all the individuals so that ind[i].myloc is in (0...NUMLOCATIONS).  ind[i].loc is unchanged
  to whatever the input file indicates*/
void GetNumLocations (struct IND *ind) {
  int maxloc=0, i, j, *freq, pos;
  for (i=0; i<NUMINDS; i++) {
    /* for now we're not dealing with unknown location */
    if (ind[i].Location < 0) {
      printf("individual %s has negative location!  locations should be >= 0\n", ind[i].Label);
      Kill();
    }
    if (ind[i].Location > maxloc) {
      maxloc = ind[i].Location;
    }
  }

  freq = malloc((maxloc+1)*sizeof(int));
  for (i=0; i<=maxloc; i++) {
    freq[i]=0;
  }
  for (i=0; i<NUMINDS; i++) {
    freq[ind[i].Location]++;
  }

  pos=0;
  for (i=0; i<=maxloc; i++) {
    if (freq[i]==0) {
      continue;
    }
    for (j=0; j<NUMINDS; j++) {
      if (ind[j].Location==i) {
        ind[j].myloc = pos;
      }
    }
    pos++;
  }
  free(freq);
  NUMLOCATIONS = pos;
}


/*---------------------------------------*/
void
InitializeSums (double *PSum, double *QSum, double *SiteBySiteSum, double *FstSum,
                int *NumAlleles, int *AncestDist, double *UsePopProbs,double *SumEpsilon)
    /*initialize arrays which store sums of parameters */
{
  int loc, pop, allele, ind, box, gen, line;

  for (loc = 0; loc < NUMLOCI; loc++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        PSum[PPos (loc, pop, allele)] = 0.0;
      }
    }
  }

  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      QSum[QPos (ind, pop)] = 0.0;
    }
  }

  if (SITEBYSITE) {
    if ((LINKAGE) && (!PHASED)) {
      for (ind = 0; ind < NUMINDS; ind++) {
        for (loc = 0; loc < MAXPOPS; loc++) {
          for (line = 0; line < LINES; line++) {
            for (pop = 0; pop < MAXPOPS; pop++) {
              SiteBySiteSum[DiploidSiteBySiteSumPos (ind, line, loc, pop)] = 0.0;
            }
          }
        }
      }
    } else {
      for (ind = 0; ind < NUMINDS; ind++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
          for (line = 0; line < LINES; line++) {
            for (pop = 0; pop < MAXPOPS; pop++) {
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] = 0.0;
            }
          }
        }
      }
    }
  }

  if (ANCESTDIST) {
    for (box = 0; box < NUMBOXES; box++) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (ind = 0; ind < NUMINDS; ind++) {
          AncestDist[AncestDistPos (ind, pop, box)] = 0;
        }
      }
    }
  }

  if (USEPOPINFO) {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (gen = 0; gen <= GENSBACK; gen++) {
          UsePopProbs[UsePPrPos (ind, pop, gen)] = 0.0;
        }
      }
    }
  }

  if (FREQSCORR) {
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        SumEpsilon[ EpsPos(loc,allele)] = 0.0;
      }
    }
    for (pop = 0; pop < MAXPOPS; pop++) {
      FstSum[pop] = 0.0;
    }
  }
}


/*-----------------------------------------------------*/
void
InitializeR (double *R, double *sumR, double *varR)
{
  int ind;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (LOG10RSTART> LOG10RMIN && LOG10RSTART<LOG10RMAX) {
      R[ind]=exp(LOG10RSTART*2.302585092994046);
    } else {
      R[ind] = exp((LOG10RMIN+LOG10RMAX)*1.15129);
    }
  }
  sumR[ind] = 0.0;
  varR[ind] = 0.0;
}
/*---------------------------------------*/
void InitializeGeno (int *Geno, int *PreGeno)
{
  int ind, line, loc;
  for (ind = 0; ind < NUMINDS; ind++) {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        Geno[GenPos (ind, line, loc)] = PreGeno[GenPos (ind, line, loc)];
      }
    }
  }
}

/*---------------------------------------*/
void Initialization (int *Geno, int *PreGeno,
                     struct IND *Individual, int *Translation,
                     int *NumAlleles, int *Z, int *Z1, double *Epsilon,
                     double *SumEpsilon,
                     double *Fst,double *PSum, double *Q, double *QSum,
                     double *SiteBySiteSum, double *FstSum,
                     int *AncestDist, double *UsePopProbs, double *Alpha,
                     double *sumAlpha, double *sumR, double *varR,
                     double *sumlikes, double *sumsqlikes,
                     int *savefreq, double *R, double *lambda, double *sumlambda,
                     double *Phase, int *Recessive, double *LocPrior,
                     double *sumLocPrior, int LocPriorLen,
                     double *sumIndLikes, double *indlike_norm)

    /*
      This function is in charge of initializing the data arrays and other
      parameters to appropriate values */
{
  int pop, ind, loc, i;

  NOALPHA = 0;   /*initialize the alphas*/
  for (pop=0; pop<MAXPOPS; pop++) {
    Alpha[pop] = ALPHA;
    sumAlpha[pop] = 0.0;
    lambda[pop]=LAMBDA;
    sumlambda[pop]=0.0;
  }

  if (LOCPRIOR) {
    if (NOADMIX==0) {
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        for (pop=0; pop<MAXPOPS; pop++) {
          Alpha[AlphaPos(loc, pop)] = ALPHA;
          sumAlpha[AlphaPos(loc, pop)] = 0.0;
        }
      }
    }
    for (i=0; i<LocPriorLen; i++) {
      sumLocPrior[i]=0.0;
      LocPrior[i] = 1.0/(double)MAXPOPS;
    }
    LocPrior[0] = LOCPRIORINIT;
  }

  *sumlikes = 0.0;
  *sumsqlikes = 0.0;

  if (INTERMEDSAVE > 0) {
    *savefreq = (int) NUMREPS / (INTERMEDSAVE + 1);
    /*need to worry about saving one time too many due to integer truncation */
    if (((*savefreq) * (INTERMEDSAVE + 1)) < NUMREPS) {
      (*savefreq)++;
    }
  } else {
    *savefreq = 0;
  }

  if (LINKAGE) {
    InitializeR (R, sumR, varR);
  }

  if (RECESSIVEALLELES) {
    CountAlleles (PreGeno, NumAlleles, Translation, Recessive);
    InitializeGeno (Geno, PreGeno);
  } else {
    CountAlleles (Geno, NumAlleles, Translation, Recessive);  /*recode alleles to {0,..,1-k} */
  }

  if (STARTATPOPINFO) {
    InitFromGeogPop (Geno, Individual, Z, 1);    /*set Z to geog origin */
  } else {
    InitializeZ (Geno, Individual, Z);  /*set Z to random initial values */
  }

  InitFreqPriors (Epsilon, Fst, Geno, NumAlleles);      /*set priors on allele freqs */
  InitializeSums (PSum, QSum, SiteBySiteSum, FstSum, NumAlleles, AncestDist, UsePopProbs,SumEpsilon);

  for (ind=0; ind<NUMINDS; ind++) {
    for (pop=0; pop<MAXPOPS; pop++) {
      Q[QPos(ind, pop)] = 1.0/MAXPOPS;
    }
  }

  if (USEPOPINFO) {
    CheckPopPriors (Individual);
    InitFromGeogPop (Geno, Individual, Z, 0);
  }

  /* bug -- corrected by Daniel April 2004
     if ((LINKAGE) && (!(PHASED)))
     for (ind = 0; ind < NUMINDS; ind++)
     for (loc = 0; loc < NUMLOCI; loc++)
     Phase[PhasePos (ind, loc)] = 0.5;
  */
  if ((LINKAGE) && (!(PHASED))&&(! (PHASEINFO))) {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        Phase[PhasePos (ind, loc)] = 0.5;
      }
    }
  }

  for (ind=0; ind<NUMINDS; ind++) {
    sumIndLikes[ind] = indlike_norm[ind] = 0.0;
  }
}


/*-----------------------------------------*/
double LogProbQ (double *Q, double onealpha, struct IND *Individual)
{
  /*return log prob of q given alpha [for single alpha in all populations].
    See notes 5/13/99 */
  double sum;
  double runningtotal;
  int ind, pop;
  int numinds = 0;              /*this is the number of individuals without pop. info */
  double sqrtunder = sqrt (UNDERFLO);

  sum = 0.0;

  runningtotal = 1.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */
      numinds++;
      for (pop = 0; pop < MAXPOPS; pop++) {
        /*being more careful with underflow caused by very small values of Q */
        if (Q[QPos (ind, pop)] > sqrtunder) {             /*0-values lead to underflow */
          runningtotal *= Q[QPos (ind, pop)];
        } else {
          runningtotal *= sqrtunder;
        }

        if (runningtotal < sqrtunder)  {  /*this is to avoid having to take logs all the time */
          if (runningtotal == 0.0) {
            printf ("*");
          }
          sum += (onealpha - 1.0) * log (runningtotal);
          runningtotal = 1.0;
        }
      }
    }
  }

  sum += (onealpha - 1.0) * log (runningtotal);
  sum += (mylgamma (MAXPOPS * onealpha) - MAXPOPS * mylgamma (onealpha)) * numinds;

  /*printf("%1.2e ",sum);
    printf("%d ",numinds); */
  return (sum);
}


/*-----------------------------------------*/
double LogProbQonepop (double *Q, double popalpha, double sumalphas, struct IND *Individual,int pop)
{
  /*return log prob of q given alpha--for one element of q ONLY.  This version is for
    updates where there is one alpha for each population.  Everything cancels out of
    the M-H ratio except one term of the gamma function, top and bottom, and the
    relevant product in the q's */

  double sum;
  double runningtotal;
  int ind;
  int numinds = 0;              /*this is the number of individuals without pop. info */
  double sqrtunder = sqrt (UNDERFLO);
  sum = 0.0;
  runningtotal = 1.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */

      numinds++;
      /*being more careful with underflow caused by very small values of Q */
      /*0-values lead to underflow */

      if (Q[QPos (ind, pop)] > sqrtunder) {
        runningtotal *= Q[QPos (ind, pop)];
      } else {
        runningtotal *= sqrtunder;
      }

      if (runningtotal < sqrtunder) {    /*this is to avoid having to take logs all the time */
        if (runningtotal == 0.0) {
          printf ("*");
        }
        sum += (popalpha - 1.0) * log (runningtotal);
        runningtotal = 1.0;
      }
    }
  }

  sum += (popalpha - 1.0) * log (runningtotal);
  sum += (mylgamma (sumalphas) - mylgamma (popalpha)) * numinds;
  return (sum);
}


/*-----------------------------------------*/
double AlphaPriorDiff (double newalpha, double oldalpha)
{
  /*returns log diff in priors for the alpha, assuming a gamma prior on alpha
    See notes 7/29/99 */
  return ((ALPHAPRIORA - 1) * log (newalpha / oldalpha) + (oldalpha - newalpha) / ALPHAPRIORB);
}


/*-----------------------------------------*/
void UpdateAlpha (double *Q, double *Alpha, struct IND *Individual, int rep)
{
  /* Produce new *Alpha using metropolis step.  There are two cases
     here: either there is the same alpha for all populations, or we do a
     separate Metropolis update for the alpha for each population.*/

  double newalpha;
  /*  double logoldprob;
   *  double lognewprob; */
  double unifrv;
  double threshold;
  double logprobdiff = 0;
  double sumalphas;
  int pop, numalphas,i;

  if (!((NOADMIX) && ((rep >= ADMBURNIN) || (rep > BURNIN)))) {
    /*don't update alpha in these cases*/
    if (POPALPHAS) {
      numalphas = MAXPOPS;
    }
    else numalphas = 1;
    for (pop = 0; pop < numalphas; pop++) {
      newalpha = RNormal (Alpha[pop], ALPHAPROPSD); /*generate proposal alpha */

      /*reject immed. if out of range*/
      if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {
        if (!(UNIFPRIORALPHA)) {
          logprobdiff = AlphaPriorDiff (newalpha, Alpha[pop]);
        }
        /*compute probabilities */
        if (POPALPHAS)  { /*different alphas in each population*/
          sumalphas = 0.0;  /*need to send in sum of alphas*/
          for (i=0; i<MAXPOPS; i++)  {
            sumalphas += Alpha[i];
          }

          /*compute probabilities for M-H ratio*/
          logprobdiff -= LogProbQonepop (Q, Alpha[pop], sumalphas,Individual,pop);
          sumalphas += newalpha - Alpha[pop];
          logprobdiff += LogProbQonepop (Q, newalpha, sumalphas, Individual,pop);
        } else  {  /*same alpha for all populations*/
          logprobdiff += LogProbQ (Q, newalpha, Individual);
          logprobdiff -= LogProbQ (Q, Alpha[pop], Individual);
        }

        /*accept new alpha with min of 1 and exp(logprobdiff) */
        threshold = exp (logprobdiff);
        unifrv = rnd ();

        /*printf("%d %.3f %.3f %.4f     ",pop,Alpha[pop],newalpha,threshold);
          if (pop==MAXPOPS-1) printf("\n");*/

        if (unifrv < threshold) {
          Alpha[pop] = newalpha;

          if (!(POPALPHAS)) { /*if same alpha in all populations*/
            for (pop = 1; pop < MAXPOPS; pop++) {
              Alpha[pop] = newalpha;
            }
          }
        }
      }
    }
  }
}

/*--------------------------------------------*/

/* returns log(exp(a)+exp(b)) without causing overflow issues */
double logsumexp(double a, double b)
{
  if (a-b > 100) {
    return a;
  }
  if (b-a > 100) {
    return b;
  }
  return a + log(1.0 + exp(b-a));
}

/* returns log Pr(Q|LocPrior) for a subset of individuals at a location */
double LogProbQ_LocPrior_loc(double *Q, double *Alpha, struct IND *Individual, int loc)
{
  double sumalpha=0.0, sumgammaalpha=0.0, like=0.0;
  int ind, pop, numind=0;

  for (ind=0; ind<NUMINDS; ind++) {
    if (Individual[ind].myloc!=loc) {
      continue;
    }
    for (pop=0; pop<MAXPOPS; pop++) {
      like += (Alpha[pop]-1.0)*log(Q[QPos(ind, pop)]);
    }
    numind++;
  }

  for (pop=0; pop<MAXPOPS; pop++) {
    sumalpha += Alpha[pop];
    sumgammaalpha += mylgamma(Alpha[pop]);
  }
  like += (mylgamma(sumalpha) - sumgammaalpha)*(double)numind;

  return like;
}


/* updates Alpha under LocPrior model */
void UpdateAlphaLocPrior(double *Q, double *Alpha, double *LocPrior,
                         struct IND *Individual)
{
  double diff, newalpha, oldalpha, lprobQ, globalpha, new_lprobQ;
  int pop, loc, pos;

  /* first update global alpha */
  for (pop=0; pop < MAXPOPS; pop++) {
    oldalpha = Alpha[pop];
    newalpha = RNormal(oldalpha, ALPHAPROPSD);
    if (newalpha >= ALPHAMAX || newalpha <= 0.0) {
      continue;
    }
    diff = 0.0;
    for (loc=0; loc<NUMLOCATIONS; loc++) {
      diff += (newalpha-oldalpha)*LocPrior[0]*log(Alpha[AlphaPos(loc,pop)]) - mylgamma(newalpha*LocPrior[0]) + mylgamma(oldalpha*LocPrior[0]) + (newalpha-oldalpha)*LocPrior[0]*log(LocPrior[0]);
    }

    if (diff > 0.0 || RandomReal(0,1) < exp(diff)) {
      Alpha[pop] = newalpha;
    }
  }

  /* now update location-specific alphas */
  for (loc=0; loc<NUMLOCATIONS; loc++) {
    pos = AlphaPos(loc, 0);
    lprobQ = LogProbQ_LocPrior_loc(Q, &Alpha[pos], Individual, loc);
    for (pop=0; pop<MAXPOPS; pop++) {
      globalpha = Alpha[pop];
      oldalpha = Alpha[pos+pop];
      newalpha = RNormal(oldalpha, ALPHAPROPSD);

      if (newalpha <= 0.0) {
        continue;
      }
      Alpha[pos+pop] = newalpha;
      new_lprobQ = LogProbQ_LocPrior_loc(Q, &Alpha[pos], Individual, loc);
      diff = (globalpha*LocPrior[0]-1.0)*log(newalpha/oldalpha) - LocPrior[0]*(newalpha-oldalpha) + new_lprobQ - lprobQ;
      if (diff >= 0.0 || RandomReal(0,1) < exp(diff)) {
        lprobQ = new_lprobQ;
      } else {
        Alpha[pos+pop] = oldalpha;
      }
    }
  }
}


/*--------------------------------------------*/

void UpdatePopLambda (double *LogP, double *lambda, int *NumAlleles)
    /*updates a lambda for each population*/
{
  double new;
  double sum;
  double sumlogp;
  int loc, pop, allele;

  for (pop=0;pop<MAXPOPS;pop++) {
    new = RNormal (lambda[pop], LAMBDAPROPSD); /*proposal*/

    if ((new > 0.0) && (new < LAMBDAMAX)) {
      sum = 0.0;
      for (loc=0; loc < NUMLOCI; loc++) {   /*compute log of likelihood ratio*/
        if (NumAlleles[loc] > 1) {
          /*norm constants*/
          sum +=  mylgamma((double) NumAlleles[loc]*new);
          sum -=  mylgamma((double) NumAlleles[loc]*lambda[pop]);
          sum +=  (double) NumAlleles[loc] * mylgamma(lambda[pop]);
          sum -=  (double) NumAlleles[loc] * mylgamma(new);

          /*printf("%d %1.3f ----- ",loc,sum);*/

          sumlogp = 0.0;
          for (allele=0; allele<NumAlleles[loc]; allele++) {
            sumlogp += LogP[PPos(loc,pop,allele)];
          }
          sum += (new - lambda[pop])*sumlogp;
          /*printf("%1.3f\n",sum);*/
        }
      }

      if (rnd() < exp(sum)) {
        lambda[pop]=new;
      }
    }
  }
}


/*--------------------------------------------*/
void UpdateLambda (double *LogP,double *Epsilon, double *lambda, int *NumAlleles)
    /*updates single value of lambda.  If FREQSCORR is turned on, this is based on the
      ancestral frequencies (Epsilon); otherwise it is based on ALL
      the population frequencies, P.  Uniform prior for lambda assumed */
{
  double new;
  double sum;
  double sumlogp;
  int loc, pop, allele,stoppop;

  new = RNormal (lambda[0], LAMBDAPROPSD); /*proposal*/

  if ((new > 0.0) && (new < LAMBDAMAX)) {
    if (FREQSCORR) {
      stoppop=1;
    } else {
      stoppop = MAXPOPS;
    }

    sum = 0.0;
    for (loc=0; loc < NUMLOCI; loc++) {  /*compute log of likelihood ratio*/
      if (NumAlleles[loc] > 1) {
        /*norm constants*/
        sum += (double) stoppop * mylgamma((double) NumAlleles[loc]*new);
        sum -= (double) stoppop * mylgamma((double) NumAlleles[loc]*lambda[0]);
        sum += (double) stoppop * (double) NumAlleles[loc] * mylgamma(lambda[0]);
        sum -= (double) stoppop * (double) NumAlleles[loc] * mylgamma(new);
        /*printf("%d %1.3f ----- ",loc,sum);*/

        sumlogp = 0.0;
        for (pop=0; pop<stoppop; pop++) {
          for (allele=0; allele<NumAlleles[loc]; allele++) {
            if (FREQSCORR) {
              sumlogp += log(Epsilon[EpsPos(loc,allele)]);
            } else {
              sumlogp += LogP[PPos(loc,pop,allele)];
            }
          }
        }
        sum += (new - lambda[0])*sumlogp;
        /*printf("%1.3f\n",sum);*/
      }
    }

    if (rnd() < exp(sum)) {
      for (pop=0;pop<MAXPOPS;pop++) {
        lambda[pop]=new;
      }
    }
  }
}


/*============================================*/
double
FPriorDiff (double newf, double oldf)
{
  /*returns log diff in priors for the correlation, f. See notes 5/14/99, and 7/15/99 */

  return ((FPRIORMEAN*FPRIORMEAN/(FPRIORSD*FPRIORSD) - 1) * log (newf / oldf) + (oldf - newf) *FPRIORMEAN/(FPRIORSD*FPRIORSD));

}


/*-----------------------------------------*/
double
FlikeFreqs (double f, double *Epsilon, double *LogP, int *NumAlleles, int pop)
{
  /*returns the log probability of the allele frequencies (for a particular pop)
    given the prior and the z (number of each allele in each population).
    Here f is the value of Fst for that locus */

  /* If numalleles=1 this seems to be ok
   * here passes Epsilon into mylgamma. Does this cause problems if epsilon very small? */
  int allele;
  double sum;
  int loc;
  double frac = (1.0-f)/f;

  sum = NUMLOCI*mylgamma(frac);
  for (loc=0; loc<NUMLOCI; loc++) {
    for (allele=0; allele < NumAlleles[loc]; allele++) {
      sum += frac*Epsilon[EpsPos (loc, allele)]*LogP[PPos(loc,pop,allele)];
      sum -= mylgamma( frac*Epsilon[EpsPos (loc, allele)]);
    }
    if (NumAlleles[loc]==0) {
      sum -=mylgamma(frac); /* should not be counting sites with all missing data */
    }
  }
  return sum;
}


/*-----------------------------------------*/
void
UpdateFst (double *Epsilon, double *Fst,
           double *LogP, int *NumAlleles)
    /*update the correlation factor, Fst, for each population*/
{

  double newf,oldf;
  double logprobdiff;
  double unifrv;
  double threshold;
  int pop1,pop2;
  int numpops1, numpops2;

  /*------Update f ()----See notebook, 5/14/99-----------*/

  /*There are two models: either there is a different F for each population,
    in which case we go through the entire loop K times; otherwise there
    is a single F, in which case we sum the likelihood ratio across populations.*/

  /*control the outer loop*/
  if (ONEFST) {
    numpops1 = 1;
  } else {
    numpops1 = MAXPOPS;
  }

  for (pop1 = 0; pop1 < numpops1; pop1++) {
    /*generate proposal f */
    oldf = Fst[pop1];
    newf = RNormal (oldf, FPRIORSD);

    /*reject if propopal < 0 or greater than 1 */
    if (newf > 0.0 && newf<1.0) {
      /*compute prior ratio */
      logprobdiff = FPriorDiff (newf, oldf);

      /*compute log likelihood diff */
      if (ONEFST) {
        numpops2 = MAXPOPS;
      } else {
        numpops2 = pop1+1;
      }
      for (pop2 = pop1; pop2 < numpops2; pop2++){
        logprobdiff += FlikeFreqs (newf, Epsilon, LogP, NumAlleles, pop2);
        logprobdiff -= FlikeFreqs (oldf, Epsilon, LogP, NumAlleles, pop2);
      }

      /*decide whether to accept, and then update*/

      if (logprobdiff >= 0.0) {   /*accept new f */
        for (pop2 = pop1; pop2 < numpops2; pop2++) {
          Fst[pop2] = newf;
        }
      } else {                 /*accept new parameter with prob p */
        threshold = exp (logprobdiff);
        unifrv = rnd ();
        if (unifrv < threshold) {
          for (pop2 = pop1; pop2 < numpops2; pop2++) {
            Fst[pop2] = newf;
          }
        }
      }
    }
  }
}



/*===============================================*/
void
GetNumFromPop (int *NumAFromPop, int *Geno, int *Z, int loc, int numalleles,struct IND *Individual)
{
  /*Fill in the number of each allele from each pop */
  int ind, line, pop, allele;
  /* int genpos; */
  int allelevalue;
  int popvalue;

  for (pop = 0; pop < MAXPOPS; pop++) {
    for (allele = 0; allele < numalleles; allele++) {
      NumAFromPop[NumAFromPopPos (pop, allele)] = 0;
    }
  }

  if (PFROMPOPFLAGONLY) {     /*this option uses only individuals with POPFLAG=1 to update P*/
    for (ind = 0; ind < NUMINDS; ind++) {
      if (Individual[ind].PopFlag == 1) {    /*individual must have popflag turned on*/
        for (line = 0; line < LINES; line++) {
          popvalue = Z[ZPos (ind, line, loc)];
          allelevalue = Geno[GenPos (ind, line, loc)];

          if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
            NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
          }
        }
      }
    }
  } else {       /*standard update--use everybody to update P */
    for (ind = 0; ind < NUMINDS; ind++) {
      for (line = 0; line < LINES; line++) {
        popvalue = Z[ZPos (ind, line, loc)];
        allelevalue = Geno[GenPos (ind, line, loc)];

        if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
          NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
        }
      }
    }
  }
}


/*------------------------------------------*/
void
GetNumLociPop (int *NumLociPop, int *Z, int ind)
{
  /*Fill in the number of alleles that each individual has from each pop */
  int loc, line, pop;
  int from;

  for (pop = 0; pop < MAXPOPS; pop++) {
    NumLociPop[pop] = 0;
  }

  for (loc = 0; loc < NUMLOCI; loc++) {
    for (line = 0; line < LINES; line++) {
      from = Z[ZPos (ind, line, loc)];
      if (from != UNASSIGNED)
        NumLociPop[from]++;
    }
  }
}


/*------------------------------------------*/
void IndependenceUpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double Lambda)
    /*this is the alternative update to the one below, proposed by Graham */
{
  int loc, pop, allele;
  /*  double difference; */
  double Sum;
  double frac;
  double *trialepsilon,*parameters;

  trialepsilon = calloc (MAXALLELES, sizeof (double));
  parameters=calloc(MAXALLELES,sizeof(double));

  if (trialepsilon == NULL || parameters == NULL) {
    printf ("warning: unable to allocate memory in UpdateEpsilon\n");
    Kill ();
  }

  for (loc = 0; loc < NUMLOCI; loc++) {
    if (NumAlleles[loc] > 1) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        parameters[allele] = Lambda;
        for (pop = 0; pop < MAXPOPS; pop++) {
          parameters[allele] +=
              (1.0 - Fst[pop]) * P[PPos (loc, pop, allele)] / Fst[pop];
        }
      }

      RDirichlet (parameters, NumAlleles[loc], trialepsilon);
      Sum = 0.0;
      /* compute Hastings (transition) ratio and prior ratio*/
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Sum +=
            (parameters[allele] - Lambda) * (log(Epsilon[EpsPos (loc, allele)]) -
                                             log(trialepsilon[allele]));
      }

      /*compute likelihood ratio*/

      for (pop = 0; pop < MAXPOPS; pop++) {
        frac = (1.0 - Fst[pop]) / Fst[pop];
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Sum +=
              mylgamma (frac * Epsilon[EpsPos (loc, allele)]);
          Sum -=
              mylgamma (frac * trialepsilon[allele]);
          Sum +=
              frac * (trialepsilon[allele] - Epsilon[EpsPos (loc, allele)])
              * LogP[PPos (loc,pop,allele)];
        }
      }

      if (rnd () < exp (Sum)) {
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Epsilon[EpsPos (loc, allele)] = trialepsilon[allele];
        }
      }
    }
  }
  free (trialepsilon);
  free (parameters);
}



/*------------------------------------------*/
void
UpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double lambda)
    /*update the ancestral allele freq vector Epsilon.  This is done
      by picking 2 alleles at each locus, and changing their frequencies.
      See notes May 30; June 20, 2001*/

{
  int loc,pop,allele1,allele2;
  double difference,invsqrtnuminds;
  double sum;
  double frac;
  /* double ratio; */

  /*here we choose between two different updates that we believe have different mixing
    properties, especially for small lambda. The independence update uses a
    Dirichlet prior independent of current epsilon while the update below uses a small normal jump */
  if (rnd()<0.5) {
    IndependenceUpdateEpsilon(P,LogP, Epsilon, Fst,NumAlleles, lambda);
  } else {
    /*this sets the range from which the proposal is drawn*/
    invsqrtnuminds=pow((double)NUMINDS,-0.5);

    for (loc=0;loc<NUMLOCI;loc++) {
      if (NumAlleles[loc]>1) {
        allele1=RandomInteger(0,NumAlleles[loc]-1);

        do {
          allele2=RandomInteger(0,NumAlleles[loc]-1);
        } while (allele1==allele2);

        difference=RandomReal(0,invsqrtnuminds);

        /*check that the proposals are in range*/
        if ((Epsilon[EpsPos(loc,allele1)]+difference<1.0) &&
            (Epsilon[EpsPos(loc,allele2)]-difference>0.0)) {

          sum=0.0;
          for (pop=0; pop<MAXPOPS; pop++) { /*compute likelihood ratio*/
            frac = (1.0-Fst[pop])/Fst[pop];

            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele1)]);
            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele2)]);
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele1)]+difference));
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele2)]-difference));

            sum += frac*difference*LogP[PPos (loc, pop, allele1)];
            sum -= frac*difference*LogP[PPos (loc, pop, allele2)];
          }

          if (lambda != 1.0) {              /*compute prior ratio*/
            /*TEMP: log added by JKP 6/30/03 as I think this was previously
              an error.  Now doing testing */
            sum += log(pow( (Epsilon[EpsPos (loc, allele1)] + difference)*
                            (Epsilon[EpsPos (loc, allele2)] - difference)/
                            (Epsilon[EpsPos (loc, allele1)])/
                            (Epsilon[EpsPos (loc, allele2)]), (double) lambda-1.0));
          }

          /*if (loc==3)
            {
            printf("%1.3f->%1.3f   %1.3f->%1.3f     ",Epsilon[EpsPos(loc,0)],
            Epsilon[EpsPos(loc,0)]
            +(allele1==0)*difference-(allele1==1)*difference,
            Epsilon[EpsPos(loc,1)],
            Epsilon[EpsPos(loc,1)]
            +(allele2==0)*difference-(allele2==1)*difference);
            printf("%1.3f %1.3f     MH=%1.5f\n",
            P[PPos (loc, 0, 0)],
            P[PPos (loc, 1, 0)],
            exp(sum));
            }*/


          if (rnd() < exp(sum)) {
            Epsilon[EpsPos(loc,allele1)]+=difference;
            Epsilon[EpsPos(loc,allele2)]-=difference;
          }
        }
      }
    }
  }
}


/*------------------------------------------*/
void UpdateP (double *P, double *LogP, double *Epsilon, double *Fst,
              int *NumAlleles, int *Geno, int *Z, double *lambda, struct IND *Individual)
    /*Simulate new allele frequencies from Dirichlet distribution */
{
  int loc, pop, allele;
  double *Parameters;           /*[MAXALLS] **Parameters of posterior on P */
  int *NumAFromPop;             /*[MAXPOPS][MAXALLS] **number of each allele from each pop */

  Parameters = calloc (MAXALLELES, sizeof (double));
  NumAFromPop = calloc (MAXPOPS * MAXALLELES, sizeof (int));
  
  if ((Parameters == NULL) || (NumAFromPop == NULL)) {
    printf ("WARNING: unable to allocate array space in UpdateP\n");
    Kill ();
  }

  for (loc = 0; loc < NUMLOCI; loc++) {
    /*count number of each allele from each pop */
    GetNumFromPop (NumAFromPop, Geno, Z, loc, NumAlleles[loc], Individual);
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        if (FREQSCORR) {
          Parameters[allele] = Epsilon[EpsPos (loc, allele)]
              *(1.0- Fst[pop])/Fst[pop]
              + NumAFromPop[NumAFromPopPos (pop, allele)];
        } else {
          Parameters[allele] = lambda[pop]
              + NumAFromPop[NumAFromPopPos (pop, allele)];
        }
      }
      /*return a value of P simulated from the posterior Di(Parameters) */
      LogRDirichlet (Parameters, NumAlleles[loc], P + PPos (loc, pop, 0),LogP +PPos(loc,pop,0));

      /*need to worry about underflow in UpdateEpsilon due to
        allele frequencies being set to zero---hence previously used the
        following hack, however now pass LogP instead

        for (allele=0;allele<NumAlleles[loc];allele++) if
        (P[PPos(loc,pop,allele)]<1E-20) {
        P[PPos(loc,pop,allele)]=1E-20;

        for (pop=0; pop<MAXPOPS; pop++)
        {
        printf(" loc =%d pop= %d fst=%f ",loc,pop,Fst[pop]);
        for (allele=0;allele<NumAlleles[loc];allele++)
        printf (" Epsilon= %.5E P= %.5E Parameters=  %.5E Num= %d",
        Epsilon[EpsPos(loc,allele)],P[PPos(loc,pop,allele)],
        Parameters[allele],NumAFromPop[NumAFromPopPos (pop, allele)]);
        printf("\n");
        }
        } */
    }
  }

  free (Parameters);
  free (NumAFromPop);
}


/*----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate locprior*/
void UpdateQMetro (int *Geno, int *PreGeno, double *Q, double *P,
                   double *Alpha, int rep, struct IND *Individual, int *Recessive)
    /*The goal with this function is to improve mixing when alpha is
      small.  The problem is that in that situation the Q's can't move
      efficiently.  What I do here is simply pick a proposal Q from the
      prior, and try to move there.  I move from q0->q1 with prob: Min(1,
      [P(X|q1,P)/P(X|q0,P)]
      --- Notes 28 July 99 */

    /* note that this code overlaps substantially with updateqmetrorecombine (the
     * version used for the linkage model).  Both versions should be updated in tandem. */

{
  double *PriorQ1;    /*[MAXPOPS]; */
  double *CurrentQ;             /*[MAXPOPS]; */
  double *TestQ;                /*[MAXPOPS]; */
  int pop;
  double logdiff;
  double randomnum;
  int ind;
  int numhits = 0;
  /*  int i, ok; */

  /*  PriorQ1 = calloc (MAXPOPS, sizeof (double)); */
  CurrentQ = calloc (MAXPOPS, sizeof (double));
  TestQ = calloc (MAXPOPS, sizeof (double));

  if ((CurrentQ == NULL) || (TestQ == NULL)) {
    printf ("WARNING: error in assigning memory in function UpdateQMetro\n");
    Kill ();
  }

  /*  for (pop = 0; pop < MAXPOPS; pop++)
      PriorQ1[pop] = Alpha[pop];
  */

  PriorQ1 = Alpha;

  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */

      /*-------compute/record newq and oldq----------*/
      /*ok = 1;
        do
        { */

      if (LOCPRIOR) {
        PriorQ1 = &Alpha[AlphaPos(Individual[ind].myloc, 0)];
      }

      RDirichlet (PriorQ1, MAXPOPS, TestQ);     /*return TestQ, sampled from the prior */
      /*  for (i=0;i<MAXPOPS;i++)
          if (TestQ[i]==0) { ok=0; break;}
          }
          while (ok==0); */


      if (rep == 0) {             /*If this is the first rep, Q will not be initialized */
        for (pop = 0; pop < MAXPOPS; pop++) {
          Q[QPos (ind, pop)] = (double) 1 / MAXPOPS;
        }
      }

      for (pop = 0; pop < MAXPOPS; pop++) {
        CurrentQ[pop] = Q[QPos (ind, pop)];
      }

      /*-------Do metropolis test of newq-------*/

      logdiff = 0.0;
      /* logdiff += log(TestQ[pop]) - log(CurrentQ[pop]);
         logdiff = logdiff*(alpha-1.0); removed prior prob bit */

      logdiff += CalcLikeInd (Geno, PreGeno, TestQ, P, ind, Recessive);  /*likelihood bit */
      logdiff -= CalcLikeInd (Geno, PreGeno, CurrentQ, P, ind, Recessive);

      randomnum = RandomReal (0.0, 1.0);
      if (randomnum < exp (logdiff)) {    /*accept */
        for (pop = 0; pop < MAXPOPS; pop++) {
          Q[QPos (ind, pop)] = TestQ[pop];
        }
        numhits++;
      }
      /*for (pop=0;pop<MAXPOPS; pop++)
        printf("%1.3f %1.3f    ",CurrentQ[pop],TestQ[pop]);
        if (randomnum < exp(logdiff) ) printf("  Accepted ");
        printf("\n"); */
    }
  }
  if (REPORTHITRATE) {           /*does this once every UPDATEFREQ reps */
    if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {
      printf ("Acceptance rate in UpdateQMetro %1.3f\n",
              (double) numhits / NUMINDS);
    }
  }

  /*  free (PriorQ1); */
  free (CurrentQ);
  free (TestQ);
}


/*----------------------------------------*/
void
UpdateQNoAdmix (int *Geno, double *Q, double *P, struct IND *Individual, double *LocPrior)
    /*Assign each individual to exactly one population according to the
      conditional probabilities. */
{

  int ind, line, loc, pop;
  int allele;
  double *ProbsVector;          /*[MAXPOPS] */
  double sumlogs, sum;
  double runningtotal;
  double max=0.0, prob;
  int pickedpop;

  ProbsVector = calloc (MAXPOPS, sizeof (double));

  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */
      for (pop = 0; pop < MAXPOPS; pop++) {      /*make a vector of log probs for each pop */
        sumlogs = 0;

        /*Melissa added 7/12/07*/
        if (LOCPRIOR) {
          prob = LocPrior[LocPriorPos(Individual[ind].myloc, pop)];
        } else {
          prob = 1.0;
        }

        runningtotal = prob;
        for (line = 0; line < LINES; line++) {
          for (loc = 0; loc < NUMLOCI; loc++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
              runningtotal *= P[PPos (loc, pop, allele)];
              if (runningtotal < UNDERFLO) {      /*this is to avoid having to
                                                  take logs all the time */
                sumlogs += log (runningtotal);
                runningtotal = 1;
              }
            }
          }
        }
        ProbsVector[pop] = sumlogs + log (runningtotal);
        if (pop==0 || ProbsVector[pop] > max) {
          max = ProbsVector[pop];
        }
      }

      sum = 0.0;
      for (pop=0; pop < MAXPOPS; pop++) {
        sum += (ProbsVector[pop] = exp(ProbsVector[pop]-max));
      }

      for (pop = 0; pop < MAXPOPS; pop++) {
        Q[QPos (ind, pop)] = 0.0;
      }

      pickedpop = PickAnOption (MAXPOPS, sum, ProbsVector);
      Q[QPos (ind, pickedpop)] = 1.0;

    }
  }
  
  free (ProbsVector);
}


/*-----------------------------------------*/
void UpdateQAdmixture (double *Q, int *Z, double *Alpha, struct IND *Individual)
    /*update Q: proportion of ancest of each ind in each pop. */
{
  int *NumLociPop;              /*[MAXPOPS] -- number of alleles from each pop (by ind) */
  double *Parameters;           /*[MAXPOPS] -- dirichlet parameters of posterior on Q */
  int ind, pop, loc;
  double *usealpha=NULL;

  NumLociPop = calloc (MAXPOPS, sizeof (int));
  Parameters = calloc (MAXPOPS, sizeof (double));
  if ((NumLociPop == NULL) || (Parameters == NULL)) {
    printf ("WARNING: unable to allocate array space in UpdateQAdmixture\n");
    Kill ();
  }

  if (LOCPRIOR==0) {
    usealpha=Alpha;
  }
  
  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */
      GetNumLociPop (NumLociPop, Z, ind);
      if (LOCPRIOR) {
        loc = Individual[ind].myloc;
        usealpha = &Alpha[AlphaPos(loc, 0)];
      }

      for (pop = 0; pop < MAXPOPS; pop++) {       /*compute dirichlet parameters */
        Parameters[pop] = usealpha[pop] + NumLociPop[pop];
      }      
      RDirichlet (Parameters, MAXPOPS, Q + QPos (ind, 0));      /*generate dirichlet RVs */
    }
  }
  
  free (NumLociPop);
  free (Parameters);
}



/*-----------------------------------------*/
void
UpdateQWithPops (int *Geno, double *Q, double *P, int *Z, double *Alpha,
                 int rep, struct IND *Individual, double *UsePopProbs)
    /*This version of UpdateQ is used when the USEPOPINFO flag is
      turned on, indicating that the prior information about population
      information should be taken into account.  It assumes that individuals
      are normally from the given populations, but are occasionally
      migrants, or children, grandchildren, etc of migrants.  Individuals
      can also be marked as unknown, in which case the usual method for
      estimating Q is used. Individuals who lack prior pop info are treated
      by the regular Q updates */
{
  double *Input, *Output;       /*both [MAXPOPS] */
  int ind, pop, loc, allele1, allele2, homepop;
  int pop1;
  int gen;
  double sumlogs;
  double sumsofar;
  double runningtotal;
  double *PostProbs;            /*[MAXPOPS][GENSBACK+1] */
  double *Priors;               /*[GENSBACK+1] */
  double p;
  double weights;
  double sumweights;
  double sumprobs;
  double rannum;
  double sqrtunder = sqrt (UNDERFLO);
  double like=0.0;              /*likelihood of current genotype */
  double maxprob;


  /*compute prior probs on each option */
  /*the correct-as-assigned option has weight MIGRPRIOR; the remaining
    probability is split among the other options */

  Input = calloc (MAXPOPS, sizeof (double));
  Output = calloc (MAXPOPS, sizeof (double));
  PostProbs = calloc (MAXPOPS * (GENSBACK + 1), sizeof (double));
  Priors = calloc (GENSBACK + 1, sizeof (double));
  if ((Input == NULL) || (Output == NULL) || (PostProbs == NULL) || (Priors == NULL)) {
    printf ("Warning!! Error in assigning memory in UpdateQWithPops\n");
    Kill ();
  }

  weights = 0.5;
  sumweights = 0.0;
  for (gen = 0; gen < GENSBACK + 1; gen++) {
    weights *= 2;               /*sequence goes 1, 2, 4, 8,.... */
    Priors[gen] = weights;
    sumweights += weights * (MAXPOPS - 1);
  }
  
  for (gen = 0; gen < GENSBACK + 1; gen++) {
    Priors[gen] = MIGRPRIOR * Priors[gen] / sumweights;
  }


  /*now compute ancestry */
  
  for (ind = 0; ind < NUMINDS; ind++) {
    if (Individual[ind].PopFlag) {        /*use prior population information */
      homepop = Individual[ind].Population - 1;
      for (gen = 0; gen < GENSBACK + 1; gen++) {
        if (gen>0) {
          p = exp ((gen - 1) * log (0.5));
        } else {
          p=1.0;
        }

        for (pop = 0; pop < MAXPOPS; pop++) {
          sumlogs = 0;
          runningtotal = 1;

          for (loc = 0; loc < NUMLOCI; loc++) {
            allele1 = Geno[GenPos (ind, 0, loc)];
            if (LINES>1) {
              allele2 = Geno[GenPos (ind, 1, loc)];
            } else {
              allele2=MISSING;
            }
            
            if ((allele1 != MISSING) && (allele2 != MISSING)) {
              /* no missing alleles */
              if (gen == 0) {             /*both alleles from same population */
                like = P[PPos (loc, pop, allele1)] * P[PPos (loc, pop, allele2)];
              } else if (gen > 0) {
                /*at most one allele from the outside population */
                like =
                    0.5 * p * P[PPos (loc, pop, allele1)] * P[PPos (loc, homepop, allele2)] +
                    0.5 * p * P[PPos (loc, homepop, allele1)] * P[PPos (loc, pop, allele2)] +
                    (1.0 - p) * P[PPos (loc, homepop, allele1)] * P[PPos (loc, homepop, allele2)];
              }
            } else if (allele1 != MISSING) {       /*1 allele missing */
              like = p * P[PPos (loc, pop, allele1)]
                  + ((double) 1.0 - p) * P[PPos (loc, homepop, allele1)];
            } else if (allele2 != MISSING) {
              like = p * P[PPos (loc, homepop, allele2)]
                  + ((double) 1.0 - p) * P[PPos (loc, homepop, allele2)];
            } else {
              like = 1.0;       /*both alleles missing */
            }

            if (like > sqrtunder) {      /*0-values lead to underflow */
              runningtotal *= like;
            } else {
              runningtotal *= sqrtunder;
            }
            
            if (runningtotal < sqrtunder) {      /*this is to avoid having to
                                                  take logs all the time */
              sumlogs += log (runningtotal);
              runningtotal = 1;
            }
          }

          /*note: this next bit of code is changed from
            previously, where I took off some standard
            normalizing constant, and put a real probability
            (not log prob) into PostProbs */
          sumlogs = sumlogs + log (runningtotal);
          PostProbs[PostProbsPos (pop, gen)] = sumlogs;

        }
      }

      /*find maximum log prob */
      maxprob = PostProbs[PostProbsPos (0, 0)];
      for (gen = 0; gen < GENSBACK + 1; gen++) {
        for (pop = 0; pop < MAXPOPS; pop++) {
          if (PostProbs[PostProbsPos (pop, gen)] > maxprob) {
            maxprob = PostProbs[PostProbsPos (pop, gen)];
          }
        }
      }

      /*subtract off minimum to avoid underflow-type problems,and
        exponentiate to turn into regular probabilities */
      for (gen = 0; gen < GENSBACK + 1; gen++) {
        for (pop = 0; pop < MAXPOPS; pop++) {
          PostProbs[PostProbsPos (pop, gen)] =
              exp (PostProbs[PostProbsPos (pop, gen)] - maxprob);
        }
      }
      /*compute posterior probabilities (unnormalized) */
      sumprobs = 0.0;
      for (gen = 0; gen < GENSBACK + 1; gen++) {
        for (pop = 0; pop < MAXPOPS; pop++) {
          if ((pop == homepop) && (gen == 0)) {   /* non-migrant case */
            PostProbs[PostProbsPos (pop, gen)] *= (1 - MIGRPRIOR);
            sumprobs += PostProbs[PostProbsPos (pop, gen)];
          } else if ((pop == homepop) && (gen != 0)) {     /* dealt with elsewhere */
            PostProbs[PostProbsPos (pop, gen)] = 0.0;
          } else {
            /*migrant case */
            PostProbs[PostProbsPos (pop, gen)] =
                PostProbs[PostProbsPos (pop, gen)] * Priors[gen];
            sumprobs += PostProbs[PostProbsPos (pop, gen)];
          }
          /*printf("%2.8f  ",Infer1Probs[pop][gen]); */
        }
        /*printf("\n"); */
      }
      /*printf("%d %f\n\n",homepop+1,sumprobs); */

      if (rep >= BURNIN) {
        for (gen = 0; gen < GENSBACK + 1; gen++) {
          for (pop = 0; pop < MAXPOPS; pop++) {
            UsePopProbs[UsePPrPos (ind, pop, gen)]
                += PostProbs[PostProbsPos (pop, gen)] / sumprobs;
          }
        }
      }

      /*now pick an option from the posterior... */
      rannum = RandomReal (0.0, sumprobs);
      sumsofar = 0.0;
      for (gen = 0; gen < GENSBACK + 1; gen++) {
        for (pop = 0; pop < MAXPOPS; pop++) {
          sumsofar += PostProbs[PostProbsPos (pop, gen)];
          if (rannum < sumsofar) {        /*set ancest to this guy */
            for (pop1 = 0; pop1 < MAXPOPS; pop1++) {
              Q[QPos (ind, pop1)] = 0;
            }
            p = exp (gen * log (0.5));  /*amount of migrant ancestry */
            Q[QPos (ind, pop)] = p;
            Q[QPos (ind, homepop)] += 1 - p;
            /*printf("P: %d, G: %d\n",pop+1,gen);
              printf("%1.3f %1.3f %1.3f\n",
              Ancest[ind][0],Ancest[ind][1],Ancest[ind][2]); */

            /*set remaining ancestry to home pop */
            gen = GENSBACK + 1;
            pop = MAXPOPS;      /*breaking out of the loop */
          }
        }
      /*printf("\n"); */
      }
    }
  }
  
  free (Input);
  free (Output);
  free (PostProbs);
  free (Priors);
}



/*-----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate LocPriors*/
void UpdateQ (int *Geno, int *PreGeno, double *Q, double *P, int *Z, double *Alpha, int rep,
              struct IND *Individual, double *UsePopProbs, int *Recessive, double *LocPrior)
    /*update Q: proportion of ancest of each ind in each pop. Three
      main options here.  One is update Q with admixture, one is update
      Q with no admixture, and one is to use prior population info.
      Even when USEPOPINFO is turned on, the other Q updates are still
      called, in order to deal with any individuals who lack population
      info.  The other functions are written to ignore individuals with
      pop data when USEPOPINFO is on. */
{

  if (USEPOPINFO) {               /*update with prior population information */
    UpdateQWithPops (Geno, Q, P, Z, Alpha, rep, Individual, UsePopProbs);
  }

  if (NOADMIX) {                  /*no admixture model */
    /* don't use ADMIBURNIN with LOCPRIOR models */
    if ((rep > ADMBURNIN) || (rep > BURNIN) || LOCPRIOR) {
      UpdateQNoAdmix (Geno, Q, P, Individual, LocPrior);
    } else {
      UpdateQAdmixture (Q, Z, Alpha, Individual);       /*initial burnin */
    }
  } else {
    /*admixture model */
    if (METROFREQ > 0 && rep%METROFREQ==0) {
      UpdateQMetro (Geno, PreGeno, Q, P, Alpha, rep, Individual, Recessive);
    } else {
      UpdateQAdmixture (Q, Z, Alpha, Individual);
    }
  }
  
}




/*------------------------------------------*/
/*Melissa added 7/12/07*/
void UpdateLocPriorNoAdmix(double *Q, double *LocPrior,
                           struct IND *Individual) {
  int i, j;
  double **dcount, newpp, diff;
  double logdiff, e1, e2, g1, g2, *eta;
  int loc, pop1, pop2;

  /* get counts in each location, cluster */
  dcount = malloc(NUMLOCATIONS*sizeof(double*));
  for (i=0; i<NUMLOCATIONS; i++) {
    dcount[i] = malloc(MAXPOPS*sizeof(double));
    for (j=0; j<MAXPOPS; j++) {
      dcount[i][j] = 0.0;
    }
  }
  for (i=0; i<NUMINDS; i++) {
    for (j=0; j<MAXPOPS; j++) {
      dcount[Individual[i].myloc][j] += (double)Q[QPos(i, j)];
    }
  }

  /* first update r (LocPrior[0]) */
  eta = &LocPrior[1];
  newpp = RandomReal(LocPrior[0]-LOCPRIORSTEP, LocPrior[0]+LOCPRIORSTEP);
  if (newpp > 0.0 && newpp < MAXLOCPRIOR) {
    logdiff = mylgamma(newpp) - mylgamma(LocPrior[0]);
    for (i=0; i<MAXPOPS; i++) {
      logdiff += mylgamma(LocPrior[0]*eta[i]) - mylgamma(newpp*eta[i]);
    }

    logdiff *= NUMLOCATIONS;
    for (i=0; i<NUMLOCATIONS; i++) {
      for (j=0; j<MAXPOPS; j++) {
        logdiff += (newpp-LocPrior[0])*eta[j]*log(LocPrior[LocPriorPos(i,j)]);
      }
    }

    if (logdiff >= 0.0 || RandomReal(0,1) < exp(logdiff)) {
      LocPrior[0] = newpp;
    }
  }

  /* now update eta */
  pop1 = RandomInteger(0, MAXPOPS-1);
  pop2 = RandomInteger(0, MAXPOPS-2);
  if (pop2>=pop1) {
    pop2++;
  }
  
  diff = RandomReal(0, ALPHAPROPSD);
  e1 = eta[pop1]+diff;
  e2 = eta[pop2]-diff;
  if (e1 < 1.0 && e2 >0.0) {
    /* don't need to consider prior on alpha if using Dirichlet(1,1,..,1) */
    logdiff = NUMLOCATIONS*(mylgamma(LocPrior[0]*eta[pop1]) + mylgamma(LocPrior[0]*eta[pop2])
                            -mylgamma(LocPrior[0]*e1)
                            -mylgamma(LocPrior[0]*e2));
    for (i=0; i<NUMLOCATIONS; i++) {
      logdiff += LocPrior[0]*((e1-eta[pop1])*log(LocPrior[LocPriorPos(i, pop1)]) +
                              (e2-eta[pop2])*log(LocPrior[LocPriorPos(i, pop2)]));
    }

    if (logdiff >= 0.0 || RandomReal(0, 1) < exp(logdiff)) {
      eta[pop1] = e1;
      eta[pop2] = e2;
    }
  }

  /* now update gamma */
  for (loc=0; loc<NUMLOCATIONS; loc++) {
    pop1 = RandomInteger(0, MAXPOPS-1);
    pop2 = RandomInteger(0, MAXPOPS-2);
    if (pop2 >= pop1) {
      pop2++;
    }
    diff = RandomReal(0, ALPHAPROPSD);
    g1 = LocPrior[LocPriorPos(loc, pop1)]+diff;
    g2 = LocPrior[LocPriorPos(loc, pop2)]-diff;
    if (g1 < 1.0 && g2 > 0.0) {
      logdiff =
          (LocPrior[0]*eta[pop1]-1.0)*(log(g1)-log(LocPrior[LocPriorPos(loc, pop1)])) +
          (LocPrior[0]*eta[pop2]-1.0)*(log(g2)-log(LocPrior[LocPriorPos(loc, pop2)])) +
          dcount[loc][pop1]*(log(g1)-log(LocPrior[LocPriorPos(loc, pop1)])) +
          dcount[loc][pop2]*(log(g2)-log(LocPrior[LocPriorPos(loc, pop2)]));
      if (logdiff >= 0.0 || RandomReal(0,1)  < exp(logdiff)) {
        LocPrior[LocPriorPos(loc, pop1)] = g1;
        LocPrior[LocPriorPos(loc, pop2)] = g2;
      }
    }
  }
  for (i=0; i<NUMLOCATIONS; i++) {
    free(dcount[i]);
  }
  free(dcount);
}


void UpdateLocPriorAdmix(double *Q, double *LocPrior, double *Alpha, struct IND *Individual) {
  double diff;
  double currPP, newPP, PPdiff, globalpha, localpha, alpharat;
  int k, loc;

  diff = 0.0;
  currPP = LocPrior[0];
  newPP = RandomReal(currPP - LOCPRIORSTEP, currPP + LOCPRIORSTEP);
  if (newPP > 0.0 && newPP < MAXLOCPRIOR) {
    PPdiff = newPP - currPP;
    diff = 0.0;
    for (k=0; k<MAXPOPS; k++) {
      globalpha = Alpha[k];
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        localpha = Alpha[AlphaPos(loc, k)];
        alpharat = localpha/globalpha;
        diff += globalpha*PPdiff*log(localpha) -
                PPdiff*localpha -
                mylgamma(globalpha*newPP) +
                mylgamma(globalpha*currPP) +
                globalpha*newPP*log(newPP) -
                globalpha*currPP*log(currPP);
      }
    }
    if (diff > 0.0 || RandomReal(0, 1) < exp(diff)) {
      LocPrior[0] = newPP;
    }
  }
}


void UpdateLocPrior(double *Q, double *LocPrior, double *Alpha, struct IND *Individual) {
  if (NOADMIX) {
    UpdateLocPriorNoAdmix(Q, LocPrior, Individual);
  } else {
    UpdateLocPriorAdmix(Q, LocPrior, Alpha, Individual);
  }
}

/*------------------------------------------*/
void UpdateGeno (int *PreGeno, int *Geno, double *P, int *Z,
                 int *Recessive, int *NumAlleles, double *Q)
    /* this function updates the imputed genotypes when the genotypes are
     * ambiguous due to recessive markers or inbreeding.  */
{

  int ind, loc, allele, line, dom, toggle,rejectioncount,allelecount,notmissingcount;
  int *AlleleUsed=NULL, *AllelePresent=NULL;
  double *AlleleProbs, Sum;
  static int RejectionThreshold = 1000000;
  /*  int pop;
   *  double temp, Sum1, Sum2; */

  if (LINES == 2) {
    AlleleProbs = calloc (4, sizeof (double));
  } else {
    AlleleProbs = calloc (MAXALLELES, sizeof (double));
    AlleleUsed = calloc (MAXALLELES, sizeof (int));
    AllelePresent = calloc (MAXALLELES, sizeof (int));
  }

  if (LINES == 2) {  /* special update for diploids; general case by rejection sampling */
    for (ind = 0; ind < NUMINDS; ind++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        if (PreGeno[GenPos (ind, 0, loc)] != MISSING
            && PreGeno[GenPos (ind, 1, loc)] != MISSING) {
          if (PreGeno[GenPos (ind, 0, loc)] ==
              PreGeno[GenPos (ind, 1, loc)]) {
            for (dom = 0; dom < 4; dom++) {
              AlleleProbs[dom] = 0.0;
            }

            if (RECESSIVEALLELES
                && Recessive[loc] != MISSING    /* bug fixed 05072007 */
                && PreGeno[GenPos (ind, 0, loc)] != Recessive[loc]) {
              AlleleProbs[0] =
                  P[PPos(loc, Z[ZPos(ind, 0, loc)], Recessive[loc])] *
                  P[PPos(loc, Z[ZPos(ind, 1, loc)], PreGeno[GenPos(ind, 0, loc)])];
              AlleleProbs[1] =
                  P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                  P[PPos(loc, Z[ZPos (ind, 1, loc)], Recessive[loc])];
            }

            AlleleProbs[2] =
                P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                P[PPos(loc, Z[ZPos (ind, 1, loc)], PreGeno[GenPos (ind, 1, loc)])];

            Sum = AlleleProbs[0] + AlleleProbs[1] + AlleleProbs[2];
            dom = PickAnOption (3, Sum, AlleleProbs);

            if (dom == 0) {
              Geno[GenPos (ind, 0, loc)] = Recessive[loc];
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 1) {
              Geno[GenPos (ind, 1, loc)] = Recessive[loc];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 2) {
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else {
              printf ("surely some mistake in UpdateGeno\n");
              Kill(); /* modify by William - kill is not a standard ANSI C function */
            }
          }
        }
      }
    }
  } else { /* general n-ploid case.  Select genotype by rejection sampling */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (ind = 0; ind < NUMINDS; ind++) {
        if (Recessive[loc]==NOTAMBIGUOUS) {
          for (line=0;line<LINES;line++) {
            Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
          }
        } else {
          allelecount=0;
          notmissingcount=0;
          
          for (allele = 0; allele < NumAlleles[loc]; allele++) {
            AllelePresent[allele] = 0;
          }
          
          for (line = 0; line < LINES; line++) {
            if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
              notmissingcount+=1;
              if (AllelePresent[PreGeno[GenPos (ind, line, loc)]]==0) {
                AllelePresent[PreGeno[GenPos (ind, line, loc)]] = 1;
                allelecount+=1;
              }
            }
          }
          
          if (allelecount==notmissingcount) {  /* if number of alleles equal to number of slots then nothing to do */
            for (line=0;line<LINES;line++) {
              Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
            }
          } else {
            toggle = 1;
            rejectioncount=0;
            while (toggle && rejectioncount < RejectionThreshold) {
              rejectioncount+=1;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                AlleleUsed[allele] = 0;
              }
              
              for (line = 0; line < LINES; line++) {
                if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    if (AllelePresent[allele] || allele == Recessive[loc]) {
                      AlleleProbs[allele] =
                          P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                    } else {
                      AlleleProbs[allele] = 0.0;
                    }
                  }
                  Sum = 0.0;
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    Sum += AlleleProbs[allele];
                  }
                  dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                  Geno[GenPos (ind, line, loc)] = dom;
                  AlleleUsed[dom] = 1;
                }
              }
              
              toggle = 0;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AlleleUsed[allele] == 0 && AllelePresent[allele] == 1) {
                  toggle = 1;
                }
              }
            }
            
            if (toggle==1) {
              /* rejection method failed, set a lower rejection threshold for future iterations */
              if (RejectionThreshold > 100) {
                RejectionThreshold /= 2;
                if (RejectionThreshold <100) {
                  RejectionThreshold = 100;
                }
              } 
              /*
                printf("sorry, STRUCTURE has tried to simulate alleles for individual %d at locus %d 1,000,000 times by a rejection method and failed", ind+1, loc+1);
                Kill();
              */
              line=0;

              /*  first fill in one copy of each present allele. */
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AllelePresent[allele] == 1) {
                  Geno[GenPos(ind,line,loc)]= allele;
                  line++;
                }
              }

              /* now fill in remaining nonmissing slots. */
              for (line=allelecount; line<notmissingcount; line++) {
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  if (AllelePresent[allele] || allele == Recessive[loc]) {
                    AlleleProbs[allele] = P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                  } else {
                    AlleleProbs[allele] = 0.0;
                  }
                }

                Sum = 0.0;
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  Sum += AlleleProbs[allele];
                }
                dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                Geno[GenPos (ind, line, loc)] = dom;
                AlleleUsed[dom] = 1;
              }

              /* now fill in missing slots */
              for (line=notmissingcount;line<LINES;line++) {
                Geno[GenPos (ind, line, loc)]=MISSING;
              }
            }
          }
        }
      }
    }
  }
  free (AlleleProbs);
  if (LINES != 2) {
    free (AlleleUsed);
    free (AllelePresent);
  }
}


/*-----------------------------------------*/
void
UpdateZ (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,int rep)
    /*update Z: population origin of each allele */
{
  int ind, line, loc, pop;
  double *Cutoffs, *Cutoffs2 /*[MAXPOPS] */ ;
  /*Cutoffs contains unnormalized probabilities of
    an allele coming from each population */
  /*Cutoffs2 shows the same thing, independent of the individual's other alleles*/
  double sum=0.0, sum2=0.0;
  int allele;

  Cutoffs = calloc (MAXPOPS, sizeof (double));
  Cutoffs2 = calloc (MAXPOPS, sizeof (double));
  if (Cutoffs == NULL || Cutoffs2 == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZ\n");
    Kill ();
  }

  for (ind = 0; ind < NUMINDS; ind++) {  /*go through all alleles in sample */
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)];

        if (allele == MISSING) {   /*Missing Data */
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
          if (SITEBYSITE && rep+1>BURNIN) {
            sum=0.0;
            sum2=1.0;
            for (pop=0;pop<MAXPOPS;pop++) {
              Cutoffs[pop]=Q[QPos(ind,pop)];
              sum+=Cutoffs[pop];
            }
          }
        } else {
          /*Data present */
          sum = 0.0;    /*compute prob of each allele being from each pop */
          sum2 = 0.0;
          for (pop = 0; pop < MAXPOPS; pop++) {
            Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
            sum += Cutoffs[pop];
          }

          if (SITEBYSITE && rep+1>BURNIN && !POSTERIOR) {
            for (pop=0;pop<MAXPOPS;pop++) {
              Cutoffs2[pop] = P[PPos (loc, pop, allele)];
              sum2 += Cutoffs2[pop];
            }
          }
          Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs);
        }
        
        if (SITEBYSITE && rep+1>BURNIN) {
          for (pop=0;pop<MAXPOPS;pop++) {
            if (POSTERIOR) {
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs[pop]/sum;
            } else {
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs2[pop]/sum2;
            }
          }
        }
      }
    }
  }

  free (Cutoffs);
  free (Cutoffs2);
}



/*---------------------------------------*/
double
Forward (int *Z, double *IndividualQ, double *P, int *Geno, double Rec, int ind,
         double *RTransitProb, double *Mapdistance, double *Phase,int *Phasemodel)
{
  long pop, pop2,  loc, line;
  double loglikelihood, temp, problinked, tempP00, tempP01,
      tempP10, tempP11,asum;
  double *sum1,*sum2;
  double sqrtunder = sqrt (UNDERFLO);
  sum1=calloc(MAXPOPS,sizeof(double));
  sum2=calloc(MAXPOPS,sizeof(double));
  if (sum1==NULL || sum2==NULL) {
    printf("WARNING: unable to allocate array space in Forward\n");
    Kill();
  }
  
  loglikelihood = 0.0;
  if (PHASED) {
    for (line = 0; line < LINES; line++) {
      /* set RTransitProb for the first locus, for each population */
      for (pop = 0; pop < MAXPOPS; pop++) {
        RTransitProb[RTransitProbPos (0, line, pop)] = IndividualQ[pop];
        if (Geno[GenPos (ind, line, 0)] != MISSING) {
          RTransitProb[RTransitProbPos (0, line, pop)] *= P[PPos (0, pop, Geno[GenPos (ind, line, 0)])];
        }
      }
      /* rest of the loci */
      for (loc = 1; loc < NUMLOCI; loc++) {
        if (Mapdistance[loc] < 0) {
          problinked=0;  /*this is the code for unlinked loci*/
        } else {
          problinked = exp (-Rec * Mapdistance[loc]);
        }
        
        /*these steps are time critical. Could save further time by storing values of exp(-Rec*mapdistance) */
        asum=0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          asum+=RTransitProb[RTransitProbPos(loc-1,line,pop)];
        }
        
        for (pop = 0; pop < MAXPOPS; pop++) {
          if (Geno[GenPos (ind, line, loc)] == MISSING) {
            RTransitProb[RTransitProbPos (loc, line, pop)] =
                (1.0-problinked)*IndividualQ[pop]*asum + problinked* RTransitProb[RTransitProbPos (loc - 1, line, pop)];
        /*necessary to incorporate recombination even if the data is missing */
          } else {
            RTransitProb[RTransitProbPos (loc, line, pop)] =
                ((1.0-problinked)*IndividualQ[pop]*asum + problinked * RTransitProb[RTransitProbPos (loc - 1, line, pop)])
                * P[PPos (loc, pop, Geno[GenPos (ind, line, loc)])];
          }
        }
        
        /*prevent underflows */
        temp = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          temp += RTransitProb[RTransitProbPos (loc, line, pop)];
        }
        if (temp < sqrtunder) {
          loglikelihood += log (sqrtunder);
          for (pop = 0; pop < MAXPOPS; pop++) {
            RTransitProb[RTransitProbPos (loc, line, pop)] /= sqrtunder;
          }
        }
      }

      temp = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        temp += RTransitProb[RTransitProbPos (NUMLOCI - 1, line, pop)];
      }
      loglikelihood += log (temp);
    }                           /* do the next LINE */
  } else {

    /* set RTransitProb for the first locus, for each population */
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        if (Geno[GenPos (ind, 0, 0)] != MISSING) {
          tempP00 = P[PPos (0, pop, Geno[GenPos (ind, 0, 0)])];
          tempP01 = P[PPos (0, pop2, Geno[GenPos (ind, 0, 0)])];
        } else {
          tempP00 = 1.0;
          tempP01 = 1.0;
        }
        
        if (Geno[GenPos (ind, 1, 0)] != MISSING) {
          tempP10 = P[PPos (0, pop, Geno[GenPos (ind, 1, 0)])];
          tempP11 = P[PPos (0, pop2, Geno[GenPos (ind, 1, 0)])];
        } else {
          tempP10 = 1.0;
          tempP11 = 1.0;
        }
        
        if (Phasemodel[ind]==1) {
          RTransitProb[DiploidRTransitProbPos (0, pop, pop2)] =
              IndividualQ[pop] * IndividualQ[pop2] *
              (Phase[PhasePos (ind, 0)] * tempP00 * tempP11 + (1.0 - Phase[PhasePos (ind, 0)]) * tempP01 * tempP10);
        } else {
          RTransitProb[DiploidRTransitProbPos (0, pop, pop2)] = IndividualQ[pop] * IndividualQ[pop2] * tempP00 * tempP11;
        }
      }
    }

    /* rest of the loci */
    for (loc = 1; loc < NUMLOCI; loc++) {
      if (Mapdistance[loc] < 0) {
        problinked=0;  /*this is the code for unlinked loci*/
      } else {
        problinked = exp (-Rec * Mapdistance[loc]);
      }

      /*these steps are time critical. Could save further time by storing values of exp(-Rec*mapdistance) */
      /* first digit of pop indicates whether it refers to  locus loc (0) or locus loc-1 (loc)
         if MARKOVPHASE=0, second indicates whether it refers to maternal (0) or paternal (1)
         otherwise second indicates whether it is first (0) or second(1) locus  */
      for (pop=0;pop<MAXPOPS;pop++) {
        sum1[pop]=0.0;
        sum2[pop]=0.0;
      }
      
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
          sum1[pop]+=RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)];
          sum2[pop2]+=RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)];
        }
      }
      
      asum=0.0;
      for (pop=0;pop<MAXPOPS;pop++) {
        asum+=sum1[pop];
      }
      
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        for (pop = 0; pop < MAXPOPS; pop++) {

          /* first digit of tempP indicates whether it refers to first (0) or second (1) listed allele copy.
             if Phasemodel=1
             Second digit indicates whether it refers to maternal (0) or paternal (1) strands
             otherwise
             second digit indicates whether it refers to first (0) or second (1) strands */

          if (Geno[GenPos (ind, 0, loc)] != MISSING) {
            tempP00 = P[PPos (loc, pop, Geno[GenPos (ind, 0, loc)])];
            tempP01 = P[PPos (loc, pop2, Geno[GenPos (ind, 0, loc)])];
          } else {
            tempP00 = 1.0;
            tempP01 = 1.0;
          }
          
          if (Geno[GenPos (ind, 1, loc)] != MISSING) {
            tempP10 = P[PPos (loc, pop, Geno[GenPos (ind, 1, loc)])];
            tempP11 = P[PPos (loc, pop2, Geno[GenPos (ind, 1, loc)])];
          } else {
            tempP10 = 1.0;
            tempP11 = 1.0;
          }

          /*note that for markov model (Phasemodel==0), phase information starts at locus 1  */
          if (Phasemodel[ind]==1) {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] = (
                problinked*problinked*RTransitProb[DiploidRTransitProbPos (loc - 1, pop, pop2)]
                + (1.0-problinked)*(1.0-problinked)*IndividualQ[pop]*IndividualQ[pop2]*asum
                +problinked*(1.0-problinked)*(IndividualQ[pop2]*sum1[pop]+IndividualQ[pop]*sum2[pop2])
                                                                      )* (Phase[PhasePos (ind, loc)] * tempP00 * tempP11 + (1.0 - Phase[PhasePos (ind, loc)]) * tempP01 * tempP10);
          } else {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] =  tempP00* tempP11*(
                problinked*problinked*
                (Phase[PhasePos(ind,loc)]* RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)]+
                 (1.0-Phase[PhasePos(ind,loc)])*RTransitProb[DiploidRTransitProbPos(loc-1,pop2,pop)])
                +(1.0-problinked)*(1.0-problinked)*
                IndividualQ[pop]*IndividualQ[pop2]*asum
                +problinked*(1.0-problinked)*
                (Phase[PhasePos (ind, loc)]*(IndividualQ[pop2]*sum1[pop]+IndividualQ[pop]*sum2[pop2])
                 + (1.0-Phase[PhasePos (ind, loc)])*(IndividualQ[pop]*sum1[pop2]+IndividualQ[pop2]*sum2[pop]))
                                                                                        );
          }
        }
      }
      
      /*prevent underflows */
      temp = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
          temp += RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
        }
      }

      if (temp < sqrtunder) {
        loglikelihood += log (sqrtunder);
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] /= sqrtunder;
          }
        }
      }
    }

    /* return the log of the sum of the RTransitProbs from the last locus, this time not summed over lines */
    temp = 0.0;
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        temp += RTransitProb[DiploidRTransitProbPos (NUMLOCI - 1, pop, pop2)];
      }
    }
    loglikelihood += log (temp);
  }
  
  free(sum1);
  free(sum2);
  return loglikelihood;
}



/*-----------------------------------------------------------*/
void
Backward (int *Z, double *SiteBySiteSum, double *IndividualQ, double Rec, int ind,
          double *Mapdistance, double *RTransitProb, int rep, int *Z1,
          double *Phase, double *P, int *Geno,int *Phasemodel)
{
  int loc, pop, line, pop2, answer;
  double sum, sum2, problinked;
  double *Cutoffs, *Cutoffs2;
  double *SquareCutoffs;
  double temp00,temp01,temp10,temp11;
  /*  double temp; */

  /*added 1 in next two lines because the minimum allowable size is 2 (see last
    bit of this function where these start to refer to chromosome strands).*/
  Cutoffs = calloc (MAXPOPS+1, sizeof (double));
  Cutoffs2 = calloc (MAXPOPS+1, sizeof (double));
  SquareCutoffs = calloc (MAXPOPS * MAXPOPS, sizeof (double));
  if (Cutoffs == NULL || SquareCutoffs == NULL || Cutoffs2 == NULL) {
    printf ("WARNING: unable to allocate array space in Backwards\n");
    Kill ();
  }
  
  if (PHASED) {
    for (line = 0; line < LINES; line++)
    {
      /*NUMLOCI-1th locus */
      sum = 0.0;
      sum2 = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++)
      {
        Cutoffs[pop] = RTransitProb[RTransitProbPos (NUMLOCI - 1, line, pop)];
        sum += Cutoffs[pop];
        Cutoffs2[pop] = P[PPos (NUMLOCI-1, pop, Geno[GenPos (ind, line, NUMLOCI-1)])];
        sum2 += Cutoffs2[pop];
      }

      Z[ZPos (ind, line, NUMLOCI - 1)] = PickAnOption (MAXPOPS, sum, Cutoffs);

      if (rep + 1 > BURNIN && SITEBYSITE) {
        for (pop = 0; pop < MAXPOPS; pop++)
          if (POSTERIOR)
            SiteBySiteSum[SiteBySiteSumPos (ind, line, NUMLOCI - 1, pop)] += Cutoffs[pop] / sum;
          else
            SiteBySiteSum[SiteBySiteSumPos (ind, line, NUMLOCI - 1, pop)] += Cutoffs2[pop] / sum2;
      }
      for (loc = NUMLOCI - 2; loc > -1; loc = loc - 1)
      {
        sum = 0.0;
        sum2 = 0.0;
        if (Mapdistance[loc+1] < 0) problinked=0;  /*this is the code for unlinked loci--JP*/
        else problinked = exp (-Rec * Mapdistance[loc + 1]);
        /*Same time saving approach as in Forward */
        /* here temp has bitten the dust */
        for (pop = 0; pop < MAXPOPS; pop++)
        {
          if (pop == Z[ZPos (ind, line, loc + 1)])
            Cutoffs[pop] = RTransitProb[RTransitProbPos (loc, line, pop)] * (problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, line, loc + 1)]]);
          else
            Cutoffs[pop] = RTransitProb[RTransitProbPos (loc, line, pop)] * (1.0 - problinked) * IndividualQ[Z[ZPos (ind, line, loc + 1)]];
          sum += Cutoffs[pop];
          Cutoffs2[pop] = P[PPos (loc, pop, Geno[GenPos (ind, line, loc)])];
          sum2 += Cutoffs2[pop];
        }
        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++)
            if (POSTERIOR)
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs[pop] / sum;
            else
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs2[pop] / sum2;
        }
        Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs);
      }
    }
  } else {
    if (Phasemodel[ind]==1) {
      /*have treated first loci and others together to avoid excessive repetition */
      for (loc = NUMLOCI - 1; loc > -1; loc = loc - 1) {

        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            SquareCutoffs[SquarePos (pop, pop2)] = RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
          }
        }

        if (loc < NUMLOCI - 1) {
          if (Mapdistance[loc+1] < 0) {
            problinked=0;  /*this is the code for unlinked loci*/
          } else { problinked = exp (-Rec * Mapdistance[loc + 1]);
          }
          
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              if (pop == Z1[ZPos (ind, 0, loc + 1)]) {
                SquareCutoffs[SquarePos (pop, pop2)] *= problinked + (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 0, loc + 1)]];
              } else {
                SquareCutoffs[SquarePos (pop, pop2)] *= (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 0, loc + 1)]];
              }

              if (pop2 == Z1[ZPos (ind, 1, loc + 1)]) {
                SquareCutoffs[SquarePos (pop, pop2)] *= problinked + (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 1, loc + 1)]];
              } else {
                SquareCutoffs[SquarePos (pop, pop2)] *= (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 1, loc + 1)]];
              }
            }
          }
        }

        sum = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            sum += SquareCutoffs[SquarePos (pop, pop2)];
          }
        }

        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              SiteBySiteSum[DiploidSiteBySiteSumPos (ind, pop2, loc, pop)] += SquareCutoffs[SquarePos (pop, pop2)] / sum;
            }
          }
        }

        answer = PickAnOption (MAXPOPS * MAXPOPS, sum, SquareCutoffs);
        Z1[ZPos (ind, 0, loc)] = answer / MAXPOPS;
        Z1[ZPos (ind, 1, loc)] = answer - MAXPOPS * (int) (answer / MAXPOPS);
      }
      /* we have determined the populations for maternal and paternal strands Z1 */
      /*Now we work out the populations for the first and second loci Z */
      /*Note that meaning of Zpos changes from (ind,pop,loc) to (ind,line,loc). */
      for (loc = 0; loc < NUMLOCI; loc++) {
        Cutoffs[0] = (Phase[PhasePos(ind,loc)])
            * P[PPos (loc, Z1[ZPos (ind, 1, loc)], Geno[GenPos (ind, 1, loc)])]
            * P[PPos (loc, Z1[ZPos (ind, 0, loc)], Geno[GenPos (ind, 0, loc)])];
        Cutoffs[1] = (1.0 - Phase[PhasePos(ind,loc)])
            * P[PPos (loc, Z1[ZPos (ind, 0, loc)], Geno[GenPos (ind, 1, loc)])]
            * P[PPos (loc, Z1[ZPos (ind, 1, loc)], Geno[GenPos (ind, 0, loc)])];
        sum = Cutoffs[0] + Cutoffs[1];
        answer = PickAnOption (2, sum, Cutoffs);
        if (answer == 0) {
          Z[ZPos (ind, 0, loc)] = Z1[ZPos (ind, 0, loc)];
          Z[ZPos (ind, 1, loc)] = Z1[ZPos (ind, 1, loc)];
        } else {
          Z[ZPos (ind, 0, loc)] = Z1[ZPos (ind, 1, loc)];
          Z[ZPos (ind, 1, loc)] = Z1[ZPos (ind, 0, loc)];
        }
      }
    } else {
      /*have treated first loci and others together to avoid excessive repetition */
      for (loc = NUMLOCI - 1; loc > -1; loc = loc - 1) {

        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            SquareCutoffs[SquarePos (pop, pop2)] = RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
          }
        }

        if (loc < NUMLOCI - 1) {
          if (Mapdistance[loc+1] < 0) {
            problinked=0;  /*this is the code for unlinked loci*/
          } else {
            problinked = exp (-Rec * Mapdistance[loc + 1]);
          }
          
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              if (pop == Z[ZPos (ind, 0, loc + 1)]) {
                temp00= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              } else {
                temp00= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              }

              if (pop2 == Z[ZPos (ind, 1, loc + 1)]) {
                temp11= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              } else {
                temp11= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              }

              if (pop2 == Z[ZPos (ind, 0, loc + 1)]) {
                temp01= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              } else {
                temp01= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              }

              if (pop == Z[ZPos (ind, 1, loc + 1)]) {
                temp10= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              } else {
                temp10= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              }
              SquareCutoffs[SquarePos(pop,pop2)]*=temp00*temp11*Phase[PhasePos(ind,loc+1)]+temp10*temp01*(1.0-Phase[PhasePos(ind,loc+1)]);
            }
          }
        }

        sum = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            sum += SquareCutoffs[SquarePos (pop, pop2)];
          }
        }

        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              SiteBySiteSum[ DiploidSiteBySiteSumPos (ind, pop2, loc, pop)] += SquareCutoffs[SquarePos (pop, pop2)] / sum;
            }
          }
        }

        answer = PickAnOption (MAXPOPS * MAXPOPS, sum, SquareCutoffs);
        Z[ZPos (ind, 0, loc)] = answer / MAXPOPS;
        Z[ZPos (ind, 1, loc)] = answer - MAXPOPS * (int) (answer / MAXPOPS);
      }
    }
  }

  free (Cutoffs);
  free (Cutoffs2);
  free (SquareCutoffs);
}



/*----------------------------------------*/
double
UpdateZandSingleR (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,
                   double *R, double *Mapdistance, int rep, double *Phase,
                   int *Z1,int *Phasemodel, double *sumIndLikes,
                   double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR, logtrialR,currentloglikelihood, trialloglikelihood, indlike;
  /* long loc; */
  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandSingleR\n");
    Kill ();
  }
  
  currentloglikelihood = 0.0;
  trialloglikelihood = 0.0;
  logtrialR = RNormal(log(R[0])/2.30259,LOG10RPROPSD);
  if (logtrialR<LOG10RMIN) {
    logtrialR=2*LOG10RMIN-logtrialR;
  }

  if (logtrialR>LOG10RMAX) {
    logtrialR=2*LOG10RMAX-logtrialR;
  }
  
  trialR=exp(2.30259*logtrialR);
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }
    indlike = Forward(Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                      Mapdistance, Phase,Phasemodel);
    currentloglikelihood += indlike;
    if (sumIndLikes!=NULL) {
      sumIndLikes[ind] += exp(indlike-indlike_norm[ind]);
    }

    Backward(Z, SiteBySiteSum, IndividualQ, R[ind], ind, Mapdistance,
             RTransitProb, rep, Z1, Phase, P, Geno,Phasemodel);

    trialloglikelihood += Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb,
                                   Mapdistance, Phase,Phasemodel);

  }
  /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[0],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
  if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {  /*Accept */
    R[0] = trialR;
    /*currentloglikelihood=trialloglikelihood;  commented out by JKP--see email from Daniel 9 Dec 02*/
  }
  
  for (ind = 0; ind < NUMINDS; ind++) {
    R[ind] = R[0];
  }
  
  free (RTransitProb);
  free (IndividualQ);

  return currentloglikelihood;
}



/*----------------------------------------*/
double
UpdateZandR (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,
             double *R, double *Mapdistance, int rep, double *Phase, int *Z1,int *Phasemodel, double *sumindlike, double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR,logtrialR, currentloglikelihood, trialloglikelihood,sumlikelihood;
  /*  long loc; */


  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandR\n");
    Kill ();
  }
  
  sumlikelihood=0.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }

    currentloglikelihood = Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                                    Mapdistance, Phase,Phasemodel);
    Backward (Z, SiteBySiteSum, IndividualQ, R[ind], ind, Mapdistance, RTransitProb,
              rep, Z1, Phase, P, Geno,Phasemodel);
    
    logtrialR = RNormal(log(R[ind])/2.30259,LOG10RPROPSD);
    if (logtrialR<LOG10RMIN) {
      logtrialR=2*LOG10RMIN-logtrialR;
    }
    if (logtrialR>LOG10RMAX) {
      logtrialR=2*LOG10RMAX-logtrialR;
    }
    trialR=exp(2.30259*logtrialR);
    trialloglikelihood = Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb, Mapdistance, Phase,Phasemodel);
    /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[ind],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
    if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {        /*Accept */
      R[ind] = trialR;
      sumlikelihood+=trialloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(trialloglikelihood-indlike_norm[ind]);
      }
    } else {
      sumlikelihood+=currentloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(currentloglikelihood-indlike_norm[ind]);
      }
    }
  }
  free (RTransitProb);
  free (IndividualQ);

  return sumlikelihood;
}



/*----------------------------------------*/
void
UpdateQMetroRecombine (int *Geno, double *Q, int *Z, double *P,
                       double *Alpha, int rep, struct IND *Individual,
                       double *Mapdistance, double *R, double *Phase,int *Phasemodel)
    /* This function does the same job as UpdateQMetro in the case
       when there is recombination.
       the code is very similar, but a new routine is nevertheless a good idea */
{
  double *PriorQ;               /*[MAXPOPS]; */
  double *CurrentQ;             /*[MAXPOPS]; */
  double *TestQ;                /*[MAXPOPS]; */
  int pop;
  double logdiff;
  double randomnum;
  int ind;
  int numhits = 0;
  /*  int i, ok; */
  double *RTransitProb;
  double *IndividualQ;

  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  PriorQ = calloc (MAXPOPS, sizeof (double));
  CurrentQ = calloc (MAXPOPS, sizeof (double));
  TestQ = calloc (MAXPOPS, sizeof (double));
  IndividualQ = calloc (MAXPOPS, sizeof (double));

  if ((PriorQ == NULL) || (CurrentQ == NULL) || (TestQ == NULL)) {
    printf ("WARNING: error in assigning memory in function UpdateQMetro\n");
    Kill ();
  }

  for (pop = 0; pop < MAXPOPS; pop++) {
    PriorQ[pop] = Alpha[pop];
  }

  for (ind = 0; ind < NUMINDS; ind++) {
    if ((USEPOPINFO) && (Individual[ind].PopFlag)) {/*set Q for inds with prior info*/
      for (pop=0;pop<MAXPOPS;pop++) {
        if (pop==Individual[ind].Population-1) {
          Q[QPos(ind,pop)]=1.0;
        } else {
          Q[QPos(ind,pop)]=0.0;
        }
      }
    } else {
      /* ie don't use individuals for whom prior pop info is used */

      /*-------compute/record newq and oldq----------*/
      /*ok = 1;
        do
        { */
      RDirichlet (PriorQ, MAXPOPS, TestQ);      /*return TestQ, sampled from the prior */
      /*  for (i=0;i<MAXPOPS;i++)
          if (TestQ[i]==0) { ok=0; break;}
          }
          while (ok==0); */


      if (rep == 0) {             /*If this is the first rep, Q will not be initialized */
        for (pop = 0; pop < MAXPOPS; pop++) {
          Q[QPos (ind, pop)] = (double) 1 / MAXPOPS;
        }
      }

      for (pop = 0; pop < MAXPOPS; pop++) {
        CurrentQ[pop] = Q[QPos (ind, pop)];
      }

      /*-------Do metropolis test of newq-------*/

      logdiff = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        IndividualQ[pop] = TestQ[pop];
      }
      logdiff += Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb, Mapdistance, Phase,Phasemodel);

      for (pop = 0; pop < MAXPOPS; pop++) {
        IndividualQ[pop] = CurrentQ[pop];
      }
      logdiff -= Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb, Mapdistance, Phase,Phasemodel);

      randomnum = RandomReal (0.0, 1.0);
      if (randomnum < exp (logdiff)) {    /*accept */
        for (pop = 0; pop < MAXPOPS; pop++) {
          Q[QPos (ind, pop)] = TestQ[pop];
        }
        numhits++;
      }
      /*for (pop=0;pop<MAXPOPS; pop++)
        printf("%1.3f %1.3f    ",CurrentQ[pop],TestQ[pop]);
        if (randomnum < exp(logdiff) ) printf("  Accepted ");
        printf("\n"); */

    }
  }
  if (REPORTHITRATE) {            /*does this once every UPDATEFREQ reps */
    if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {
      printf ("Acceptance rate in UpdateQMetro %1.3f\n",
              (double) numhits / NUMINDS);
    }
  }

  free (PriorQ);
  free (CurrentQ);
  free (TestQ);
  free (IndividualQ);
  free (RTransitProb);
}



/*---------------------------------------*/
/*void
  UpdateQRecombine (int *Geno, double *Q, int *Z, double *P,
  double *Alpha, int rep, struct IND *Individual,
  double *Mapdistance, double *R, double *Phase)
{
  int pop, ind;

  if (NOQS) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      Q[QPos (ind, pop)] = 1.0 / (double) MAXPOPS;
    }
  } else {
    UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                           Individual, Mapdistance, R, Phase,Phasemodel);
  }

}*/
/*=============MAIN======================================*/

int main (int argc, char *argv[])
{

  /*data--------- */
  int *Geno;                    /*NUMINDSxLINES: genotypes */
  double *R;                    /*NUMINDS */
  double *Mapdistance;          /*NUMLOCI */
  double *Phase;                /*NUMLOCI*NUMINDS */
  int *Phasemodel=NULL;         /*NUMINDS */
  char *Markername;             /*GENELEN*NUMLOCI */

  struct IND *Individual;       /*NUMINDS: records for each individual */
  int *Translation;             /*NUMLOCIxMAXALLELES: value of each coded allele */
  int *NumAlleles;              /*NUMLOCI: number of alleles at each locus */

  /* only used for recessive or inbreeding models: */
  int *PreGeno=NULL;           /*NUMINDSxLINESxNUMLOCI; diploid genotype if recessive alleles */
  int *Recessive=NULL;         /*NUMLOCI recessive allele at each locus, or -1 if there is none */


  /*Basic parameters */
  int *Z;                       /*NUMINDSx2xNUMLOCI: Z=pop of origin for each allele */
  int *Z1;
  double *Q;                    /*NUMINDSxMAXPOPS:  Q=ancestry of individuals */
  double *P;                    /*NUMLOCIxMAXPOPSxMAXALLELES: P=population allele freqs */
  double *LogP;                 /*NUMLOCIxMAXPOPSxMAXALLELES: log of P, used to prevent underflow */
  double *Epsilon;              /*NUMLOCIxMAXALLELES: Dirichlet parameter for allele
                                  frequencies. This is either LAMBDA (if uncorrelated), or
                                  ancestral allele freqs if they are correlated */
  double *Fst;          /*MAXPOPS: Factor multiplied by epsilon under the Fst model */
  double *Alpha;                /*MAXPOPS: Dirichlet parameter for degree of admixture.
                                  Start this at ALPHA, and possibly change
                                  (if INFERALPHA==1) */
  double *lambda;                /*Dirichlet prior parameter for allele frequencies;
                                   start this at LAMBDA, and update if INFERLAMBDA*/
  double *sumlambda;
  /*Summaries */
  int    *NumLociPop;           /*NUMINDSxMAXPOPS: Number of alleles from each pop (by ind) */
  double *PSum;                 /*NUMLOCIxMAXPOPSxMAXALLELES: sum of AlFreqs */
  double *QSum;                 /*NUMINDSxMAXPOPS:  sum of Ancestries */
  double *SiteBySiteSum=NULL;
  double *FstSum;               /*MAXPOPS:  Sum of Fst */
  double *SumEpsilon=NULL;      /*NUMLOCIxMAXALLELES: sum of ancestral allele freqs*/
  double *sumAlpha;              /*MAXPOPS*/
  double *sumR;                 /*NUMINDS */
  double *varR;                 /*NUMINDS */
  double recomblikelihood=0.0;
  double like;                  /*current likelihood value */
  double sumlikes;              /*sum of likelihood values */
  double sumsqlikes;            /*sum of squared likelihoods */


  /*Melissa added 7/12/07 for calculating DIC*/
  double *sumIndLikes, *indLikesNorm;

  int    *AncestDist=NULL;      /*NUMINDS*MAXPOPS*NUMBOXES histogram of Q values */
  double *UsePopProbs=NULL;     /*NUMINDS*MAXPOPS*(GENSBACK+1) This is used when the
                                  population info is used.  It stores the probability that an
                                  individual has each of a specified set of ancestry amounts */
  /*loop variables-------------- */
  int rep;                      /*MCMC iterations so far */
  int savefreq;                 /*frequency of saving to file */
  int ind;

  /*Melissa's new variables added 7/12/07 to use priors based on sampling location*/
  double *LocPrior=NULL, *sumLocPrior=NULL, LocPriorLen=0;


  /*=====Code for getting started=============================*/

  Welcome (stdout);             /*welcome */
  GetParams (0,argc,argv);      /*read in parameter values */

  CheckParamCombinations();     /*check that some parameter combinations are valid*/

  Mapdistance = calloc (NUMLOCI, sizeof (double));
  Phase = calloc (NUMLOCI * NUMINDS, sizeof (double));


  if (LINES ==2 && PHASED ==0) {
    Phasemodel=calloc(NUMINDS,sizeof(int));
    for (ind=0;ind<NUMINDS;ind++) {
      if (MARKOVPHASE) {
        Phasemodel[ind]=0;
      } else {
        Phasemodel[ind]=1;
      }
    }
  }
  
  lambda=calloc(MAXPOPS, sizeof (double));
  sumlambda=calloc(MAXPOPS, sizeof (double));

  Markername = calloc (GENELEN*NUMLOCI, sizeof (char));
  Geno = calloc (LINES * NUMLOCI * NUMINDS, sizeof (int));
  if (RECESSIVEALLELES) {
    PreGeno = calloc (LINES * NUMLOCI * NUMINDS, sizeof (int));
    Recessive = calloc (NUMLOCI, sizeof (int));
    if (PreGeno == NULL || Recessive == NULL) {
      printf ("Error (3) in assigning memory\n");Kill ();
    }
  }

  Individual = calloc (NUMINDS, sizeof (struct IND));
  if (Geno == NULL || Individual == NULL || Mapdistance == NULL || Markername == NULL) {
    printf ("Error in assigning memory (not enough space?)\n");
    Kill ();
  }
  Randomize(RANDOMIZE, &SEED);

  /*read in data file */
  if (RECESSIVEALLELES) {
    ReadInputFile(PreGeno, Mapdistance, Markername, Individual, Phase, Recessive);
  } else {
    ReadInputFile (Geno, Mapdistance, Markername, Individual, Phase, Recessive);
  }

  if (RECESSIVEALLELES) {
    MAXALLELES = FindMaxAlleles (PreGeno, Recessive);
  } else {
    MAXALLELES = FindMaxAlleles (Geno, Recessive);
  }


  /*=============set aside memory space=====================*/
  Translation = calloc (NUMLOCI * MAXALLELES, sizeof (int));
  NumAlleles = calloc (NUMLOCI, sizeof (int));
  Z = calloc (NUMINDS * LINES * NUMLOCI, sizeof (int));
  Z1 = calloc (NUMINDS * LINES * NUMLOCI, sizeof (int));
  Q = calloc (NUMINDS * MAXPOPS, sizeof (double));
  P = calloc (NUMLOCI * MAXPOPS * MAXALLELES, sizeof (double));
  LogP = calloc(NUMLOCI * MAXPOPS * MAXALLELES, sizeof(double));
  R = calloc (NUMINDS, sizeof (double));
  sumR = calloc (NUMINDS, sizeof (double));
  varR = calloc (NUMINDS, sizeof (double));
  Epsilon = calloc (NUMLOCI * MAXALLELES, sizeof (double));
  if (FREQSCORR) {
    SumEpsilon = calloc (NUMLOCI * MAXALLELES, sizeof (double));
  }
  Fst = calloc (MAXPOPS, sizeof (double));
  FstSum = calloc (MAXPOPS, sizeof (double));
  NumLociPop = calloc (NUMINDS * MAXPOPS, sizeof (int));
  PSum = calloc (NUMLOCI * MAXPOPS * MAXALLELES, sizeof (double));
  QSum = calloc (NUMINDS * MAXPOPS, sizeof (double));
  if (SITEBYSITE) {
    if (LINKAGE && !PHASED) {
      SiteBySiteSum = calloc (NUMINDS * MAXPOPS * NUMLOCI * MAXPOPS, sizeof (double));
    } else {
      SiteBySiteSum = calloc (NUMINDS * LINES * NUMLOCI * MAXPOPS, sizeof (double));
    }
  }


  if (ANCESTDIST) {
    AncestDist = calloc (NUMINDS * MAXPOPS * NUMBOXES, sizeof (int));
  }
  if (USEPOPINFO) {
    UsePopProbs = calloc (NUMINDS * MAXPOPS * (GENSBACK + 1), sizeof (double));
  }

  /*Melissa added 7/12/07*/
  if (LOCDATA>0 || LOCISPOP) {
    GetNumLocations(Individual);
  }

  /*Allocate the LocPrior vector.
    For no-admixture, it contains r, and the vectors nu and gamma.
    For admixture, it contains gamma.  The alphas_locals are stored with alpha global*/
  if (LOCPRIOR) {
    if (NOADMIX) {
      LocPriorLen = 1+MAXPOPS*(NUMLOCATIONS+1);
    } else {
      LocPriorLen=1;
    }
    LocPrior = malloc(LocPriorLen*sizeof(double));
    sumLocPrior = malloc(LocPriorLen*sizeof(double));
  }
  
  if (LOCPRIOR && NOADMIX==0) {
    Alpha = malloc(MAXPOPS*(NUMLOCATIONS+1)*sizeof(double));
    sumAlpha = malloc(MAXPOPS*(NUMLOCATIONS+1)*sizeof(double));
  } else {
    Alpha = calloc(MAXPOPS, sizeof (double));
    sumAlpha = calloc(MAXPOPS, sizeof (double));
  }

  /* this is for DIC */
  sumIndLikes = malloc(NUMINDS*sizeof(double));
  indLikesNorm = malloc(NUMINDS*sizeof(double));

  if ((Translation == NULL) || (NumAlleles == NULL) || (Z == NULL) || (Z1 == NULL) || (Q == NULL) ||
      (P == NULL) || (LogP==NULL) || (R == NULL) || (sumR == NULL) || (varR == NULL) || (Epsilon == NULL) ||
      (Fst == NULL) || (NumLociPop == NULL) ||
      (PSum == NULL) || (QSum == NULL) || ( SITEBYSITE && (SiteBySiteSum == NULL)) || (FstSum == NULL) ||
      ((ANCESTDIST) && (AncestDist == NULL)) ||
      ((USEPOPINFO) && (UsePopProbs == NULL))||(Alpha == NULL)||(sumAlpha==NULL)||
      ((FREQSCORR) && (SumEpsilon == NULL)) ||
      (LocPriorLen>0 && (LocPrior==NULL || sumLocPrior==NULL)) ||
      sumIndLikes==NULL || indLikesNorm==NULL) {

    printf ("Error in assigning memory (not enough space?)\n");
    FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno, PreGeno, Recessive,
            Individual, Translation, NumAlleles, Z, Z1, Q, P, LogP, R, sumR, varR, Epsilon, SumEpsilon,
            Fst, FstSum, NumLociPop, PSum, QSum, SiteBySiteSum, AncestDist, UsePopProbs, LocPrior,
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm);
    Kill ();
  }
  /*=========done setting aside memory space=====================*/

  /*initialize variables and arrays */
  Initialization (Geno, PreGeno, Individual, Translation, NumAlleles, Z, Z1, Epsilon, SumEpsilon,
                  Fst, PSum, Q, QSum, SiteBySiteSum, FstSum, AncestDist, UsePopProbs, Alpha,
                  sumAlpha, sumR, varR, &sumlikes, &sumsqlikes, &savefreq, R, lambda,
                  sumlambda,Phase,Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes, indLikesNorm);
  printf ("\n\n--------------------------------------\n\n");
  printf ("Finished initialization; starting MCMC \n");
  printf ("%d iterations + %d burnin\n\n", NUMREPS, BURNIN);
  /*=====Main MCMC loop=======================================*/

  for (rep = 0; rep < (NUMREPS + BURNIN); rep++) {

    UpdateP (P,LogP, Epsilon, Fst, NumAlleles, Geno, Z, lambda, Individual);

    if (LINKAGE && rep >= ADMBURNIN) {
      UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                             Individual, Mapdistance, R, Phase,Phasemodel);
    } else {
      UpdateQ (Geno, PreGeno, Q, P, Z, Alpha, rep, Individual, UsePopProbs, Recessive, LocPrior);
    }

    if (LOCPRIOR && UPDATELOCPRIOR) {
      UpdateLocPrior(Q, LocPrior, Alpha, Individual);
    }
    
    if (RECESSIVEALLELES) {
      UpdateGeno (PreGeno, Geno, P, Z, Recessive, NumAlleles, Q);
    /*The Zs are not correct after UpdateGeno, until UpdateZ is run */
    }
    
    if (LINKAGE && rep > ADMBURNIN) {
      if (!INDIVIDUALR) {
        recomblikelihood = UpdateZandSingleR(Z, SiteBySiteSum, Q, P, Geno,
                                             R, Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN? sumIndLikes : NULL, indLikesNorm);
      } else {
        recomblikelihood = UpdateZandR(Z, SiteBySiteSum, Q, P, Geno, R,
                                       Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN ? sumIndLikes:NULL, indLikesNorm);
      }
    } else {
      UpdateZ (Z, SiteBySiteSum, Q, P, Geno,rep);
      /*      printf("done updatez alpha[2]=%e\n", Alpha[2]); */
    }

    if (LOCPRIOR && NOADMIX==0) {
      UpdateAlphaLocPrior(Q, Alpha, LocPrior, Individual);
    } else if (INFERALPHA) {
      UpdateAlpha (Q, Alpha, Individual, rep);
    }
    
    if (INFERLAMBDA) {
      if  (POPSPECIFICLAMBDA) {
        UpdatePopLambda(LogP,lambda,NumAlleles);
      } else {
        UpdateLambda (LogP,Epsilon,lambda, NumAlleles);
      }
    }


    if (FREQSCORR) {
      UpdateEpsilon(P,LogP,Epsilon,Fst,NumAlleles,lambda[0]);
      UpdateFst (Epsilon, Fst, LogP, NumAlleles);
    }

    /*====book-keeping stuff======================*/
    if (rep + 1 > BURNIN) {
      DataCollection (Geno, PreGeno, Q, QSum, Z, Z1, SiteBySiteSum, P, PSum,
                      Fst, FstSum, NumAlleles,
                      AncestDist, Alpha, sumAlpha, sumR, varR, &like,
                      &sumlikes, &sumsqlikes, R, Epsilon,SumEpsilon,recomblikelihood,
                      lambda, sumlambda, Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes, indLikesNorm, rep);
    }
    
    if ((savefreq) && ((rep + 1) > BURNIN) && (((rep + 1 - BURNIN) % savefreq) == 0)
        && ((rep + 1) != NUMREPS + BURNIN)) {
      OutPutResults (Geno, rep + 1, savefreq, Individual, PSum, QSum,
                     SiteBySiteSum, FstSum, AncestDist, UsePopProbs, sumlikes,
                     sumsqlikes, sumAlpha, sumR, varR,
                     NumAlleles, Translation, 0, Markername, R,
                     SumEpsilon,
                     lambda,sumlambda,sumLocPrior, LocPriorLen,
                     sumIndLikes, indLikesNorm, argc,argv);
    }


    if (PRINTLIKES) {
      PrintLike (like, rep, Geno, PreGeno, Q, P,recomblikelihood);
    }
    
    if (((rep + 1) % UPDATEFREQ) == 0) {
      PrintUpdate (rep + 1, Geno, PreGeno, Alpha, Fst, P, Q, like,
                   sumlikes, sumsqlikes, NumAlleles, R, lambda,Individual,
                   recomblikelihood, Recessive, LocPrior, LocPriorLen);
    }
  }

  /*====final book-keeping====================================*/
  if ((rep % UPDATEFREQ) != 0) {
    PrintUpdate (rep, Geno, PreGeno, Alpha, Fst, P, Q, like, sumlikes,
                 sumsqlikes, NumAlleles,R, lambda, Individual,recomblikelihood,
                 Recessive, LocPrior, LocPriorLen);
  }

  OutPutResults (Geno, rep, savefreq, Individual, PSum, QSum,
                 SiteBySiteSum, FstSum, AncestDist, UsePopProbs,
                 sumlikes, sumsqlikes,
                 sumAlpha, sumR, varR, NumAlleles, Translation, 1,
                 Markername, R, SumEpsilon,
                 lambda,sumlambda,sumLocPrior, LocPriorLen,
                 sumIndLikes, indLikesNorm,
                 argc,argv);

  /*=====Closing everything down==============================*/
  FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno, PreGeno, Recessive,
            Individual, Translation, NumAlleles, Z, Z1, Q, P, LogP, R, sumR, varR, Epsilon, SumEpsilon,
            Fst, FstSum, NumLociPop, PSum, QSum, SiteBySiteSum, AncestDist, UsePopProbs, LocPrior,
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm);
  return (0);
}

/*==========================================

  Notes:


  ============================================*/
