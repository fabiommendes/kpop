/*Part of structure.c.

  This bit of the program is involved in collecting data (DataCollection)
  and printing results (OutPutResults). */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "params.h"
#include "mymath.h"

void UpdateSums (double *Q, double *QSum, int *Z, double *P, double *PSum,
                 double *Fst, double *FstSum, int *NumAlleles,
                 int *AncestDist, double *Epsilon, double *SumEpsilon,
                 double *lambda, double *sumlambda,
                 double *LocPrior, double *sumLocPrior, int LocPriorLen);
double CalcLike (int *Geno, int *PreGeno, double *Q, double *P, int *Recessive,
                 double *sumindlike, double *indlike_norm);
double EstLogProb (double sumlikes, double sumsqlikes, int reps);
double KLDiv (int pop1, int pop2, double *P, double *LogP, int *NumAlleles, int reps);
void PrintKLD (FILE * file, double *P, double *LogP, int *NumAlleles, int reps, int format);
void PrintNET (FILE * file, double *P, int *NumAlleles, int reps, int format);

void PrintBanner (int rep, double *Alpha, double *Fst, double like, double *lambda);
void PrintMainParams (FILE * file, int rep, int argc, char *argv[]);
void PrintQ (FILE * file, int *Geno, int rep, double *QSum, struct IND *Individual,
             int *AncestDist, double *UsePopProbs,double *sumR);
void PrintP (FILE * file, int rep, int *Geno, double *PSum, int *Translation,
             int *NumAlleles, double *SumEpsilon, char *Markername);
void PrintMembership (FILE * file, double *QSum, struct IND *Individual);
void PrintSiteBySite (FILE * file, double *SiteBySiteSum, int rep,
                      char *Markername, double *PSum,
                      int *NumAlleles, int *Geno, int *Translation,double *SumEpsilon, struct IND *Individual);
void PrintSequences (FILE * file, int *Geno, char *Markername, double *SiteBySiteSum, int rep, int *Translation);
void PrintSums (FILE * file, int rep, double sumlikes,
                double sumsqlikes, double *FstSum, double *sumAlpha, double *sumlambda,
                double *sumR,double *varR, struct IND *Individual,
                double *sumLocPriors, int LocPriorLen, double DIC);
void PrintGeneName(FILE * file, int loc, char *Markername);
int EqualGeneNames(int loc1,int loc2,char *Markername);


/*=================================================*/
void
DataCollection (int *Geno, int *PreGeno,
                double *Q, double *QSum, int *Z, int *Z1,
                double *SiteBySiteSum, double *P, double *PSum,
                double *Fst, double *FstSum, int *NumAlleles,
                int *AncestDist, double *Alpha, double *sumAlpha,
                double *sumR, double *varR, double *like,
                double *sumlikes, double *sumsqlikes, double *R,
                double *Epsilon, double *SumEpsilon, double recomblikelihood,
                double *lambda, double *sumlambda, int *Recessive,
                double *LocPrior, double *sumLocPrior, int LocPriorLen,
                double *sumindlikes, double *indlikes_norm, int rep)
{
  int ind, pop, loc, pos;
  UpdateSums (Q, QSum, Z, P, PSum, Fst, FstSum, NumAlleles, AncestDist,
              Epsilon, SumEpsilon, lambda, sumlambda,
              LocPrior, sumLocPrior, LocPriorLen);
  if (LINKAGE) {
    for (ind = 0; ind < NUMINDS; ind++) {
      sumR[ind] += R[ind];
      varR[ind] += R[ind] * R[ind];
    }
  }

  if (LOCPRIOR && NOADMIX==0) {
    for (pop=0; pop<MAXPOPS; pop++)
      for (loc=0; loc<=NUMLOCATIONS; loc++) {
        pos = AlphaPos(loc, pop);
        sumAlpha[pos] += Alpha[pos];
      }
  } else if (!(NOADMIX) && (!(NOALPHA))) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      sumAlpha[pop] += Alpha[pop];
    }
  }


  if (COMPUTEPROB) {
    if (LINKAGE) {
      *like = recomblikelihood;
    }

    if (rep < BURNIN) {
      *like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
    } else {
      *like = CalcLike (Geno, PreGeno, Q, P, Recessive,
                        sumindlikes, indlikes_norm);
    }
    *sumlikes += *like;
    *sumsqlikes += (*like) * (*like);
  }
  /*printf("%f %f %f\n", *like, *sumlikes, *sumsqlikes); */
}

/*---------------------------------------------------*/
void
PrintLike (double like, int rep, int *Geno, int *PreGeno, double *Q,
           double *P,double recomblikelihood,
           int *Recessive)
{
  if (rep + 1 > BURNIN) {        /*already calculated */
    printf ("%6.0f\n", like);
  } else {
    if (LINKAGE) {
      printf("%6.0f\n",recomblikelihood);
    } else {
      printf ("%6.0f\n", CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL));
    }
  }
}


/*---------------------------------------------------*/
double
CalculateRAverage (double *R)
{
  int ind;
  double temp;
  temp = 0.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    temp += R[ind];
  }
  return temp / (double) NUMINDS;
}

/*----------------------------------------------------*/
void
PrintUpdate (int rep, int *Geno, int *PreGeno, double *Alpha, double *Fst, double *P, double *Q,
             double like, double sumlikes, double sumsqlikes, int *NumAlleles,
             double *R, double *lambda, struct IND *Individual,
             double recomblikelihood, int *Recessive,
             double *LocPrior, int LocPriorLen)
{
  /*print a bunch of stuff to screen during run: rep, alpha, f, KLD, likelihood...
    Also occasionally print a header banner to define the variables. */

  double logprob=0.0;
  /*  int i;
   *  int printalign; */
  int pop;

  if ((rep < BURNIN + UPDATEFREQ) && (rep > BURNIN)) {
    printf ("\nBURNIN completed");
  }

  if (((NOADMIX) && (ADMBURNIN > 0)) && !LOCPRIOR ) {
    if ((rep < ADMBURNIN + UPDATEFREQ) && (rep >= ADMBURNIN)) {
      printf ("\nAdmixture Burnin complete.  Current alpha");
      if (POPALPHAS) {
        printf ("s = ");
        for (pop = 0; pop < MAXPOPS; pop++) {
          printf ("%1.3f ", Alpha[pop]);
        }
      } else {
        printf (" = %1.3f ", Alpha[0]);
      }
      printf ("\n\n");
    }
  }
  
  if ((LINKAGE) && (ADMBURNIN> 0)) {
    if ((rep < ADMBURNIN+UPDATEFREQ) && (rep >= ADMBURNIN)) {
      printf ("\nNo recombination Burnin complete.  Current rec = %1.3f",
              CalculateRAverage (R));
      PrintBanner (rep, Alpha, Fst, like, lambda);
    }
  }

  /*calculate some stuff */

  if (LINKAGE) {
    if (rep <= ADMBURNIN) {
      like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL); 
    } else {
      like = recomblikelihood;
    }

    if (rep >= BURNIN + 2) { /* +2 because need 2 observations for variance*/
      logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
    } else if (COMPUTEPROB) { /*not linkage model*/
      if (rep <= BURNIN + 2) {
        like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
      } else {
        logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
      }
    }
  }

  /*possibly print banner to indicate what the numbers are */
  if ((rep == UPDATEFREQ) ||
      ((rep < BURNIN) && ((rep % (BANNERFREQ * UPDATEFREQ)) == 0)) ||
      ((rep > BURNIN) && (((rep - BURNIN) % (BANNERFREQ * UPDATEFREQ)) == 0))) {
    if (rep != NUMREPS + BURNIN) {        /*don't bother for last line of output */
      if ((rep > 1) && (PRINTQSUM) && (MAXPOPS > 1)) {
        PrintMembership (stdout, Q, Individual);
      }
      PrintBanner (rep, Alpha, Fst, like, lambda);
    }
  }
  
  /*print current values to screen */
  printf ("%5d:    ", rep);

  if (LINKAGE && rep >= ADMBURNIN) {
    if (!INDIVIDUALR) {
      printf ("%1.09f  ", R[0]);
    } else {
      printf ("%1.09f   ", CalculateRAverage (R));
    }
  }

  if (PRINTLAMBDA) {
    if (POPSPECIFICLAMBDA) {
      for (pop=0;pop<MAXPOPS;pop++) {
        printf ("%1.2f    ", lambda[pop]);
      }
    } else {
      printf ("%1.2f    ", lambda[0]);
    }
  }

  if ((!(NOADMIX)) && (!(NOALPHA))) {
    if (POPALPHAS) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("%1.3f  ", Alpha[pop]);
        if (pop > 8) {
          printf (" "); /*extra space for number */
        }
      }
    } else {
      printf ("%1.3f  ", Alpha[0]);
    }
  }
  
  if (FREQSCORR) {
    printf ("  ");
    if (ONEFST) {
      printf ("%1.3f ", Fst[0]);
    } else {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("%1.3f ", Fst[pop]);
      }
    }
    printf ("   ");
  } else {
    printf ("  ");
  }

  if (LOCPRIOR) {
    printf ("%1.3f ", LocPrior[0]);
    printf ("   ");
  }

  /*currently it only net distances not KLD ones */
  if (PRINTKLD || PRINTNET) {
    PrintNET (stdout, P, NumAlleles, 1, 0);
  }
  
  if (COMPUTEPROB) {


    /*put correct # of spaces in */
    /*if (((int) log10(like)) < 1) printalign = 4;
      else printalign = 4 + ((int) log10(like));
      for (i=printalign; i<8; i++)
      printf(" "); */
    if (rep > BURNIN + 2) {
      printf ("  %.0f  ", like);
      printf ("  %.0f ", logprob);
    } else {
      printf ("  --  ");
    }
  }

  printf ("\n");

  if (rep == BURNIN) {
    printf ("\nBURNIN completed");
    PrintBanner (rep, Alpha, Fst, like, lambda);
  }
  fflush(stdout);

}
/*----------------------------------------------------*/
void
PrintBanner (int rep, double *Alpha, double *Fst, double like, double *lambda)
    /*print banner to screen during run giving variable names */
{
  int i, j, k;
  int pop;

  printf ("\n");
  for (i = 4; i < ((int) log10 (rep)); i++) {
    printf (" ");
  }
  printf (" Rep#:   ");

  if ((LINKAGE) && (rep >= ADMBURNIN)) {
    printf (" r           ");
  }

  if (PRINTLAMBDA) {
    if (POPSPECIFICLAMBDA) {
      for (pop=0;pop<MAXPOPS;pop++) {
        printf("Lambda%d ",pop+1);
      }
    } else {
      printf ("Lambda  ");
    }
  }

  if (((!(NOADMIX)) && (!(NOALPHA)))) {
    printf (" ");
    if (POPALPHAS) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        printf ("Alpha%d ", pop + 1);
      }
    } else {
      printf ("Alpha  ");
    }
    /*for (i = 7; i < ((int) log10 (Alpha[0])); i++) {
      printf (" ");
    }*/
  }

  if (FREQSCORR) {
    printf ("   ");
    if (ONEFST)
      printf ("Fst   ");
    else
      for (pop = 0; pop < MAXPOPS; pop++)
        printf (" F%d   ", pop + 1);
    printf (" ");
  }
  else
    printf (" ");

  if (LOCPRIOR)
  {
    printf ("  r     ");
  }

  if (PRINTKLD || PRINTNET)
  {
    for (j = 0; j < MAXPOPS - 1; j++)
      for (k = j + 1; k < MAXPOPS; k++)
        printf (" D%d,%d ", j + 1, k + 1);
  }

  if (COMPUTEPROB)
  {
    printf ("   Ln Like ");
    for (i = 8; i < ((int) log10 (rep)); i++)
      printf (" ");

    if (rep >= BURNIN)
      printf (" Est Ln P(D)");
  }
  printf ("\n");


}

/*----------------------------------------------------*/
double
EstLogProb (double sumlikes, double sumsqlikes, int reps)
    /*returns the current estimated Prob of Data.  Reps is the
      number of reps, not including burnin */
{
  double mean = sumlikes / reps;
  double var = SampleVar (sumsqlikes, sumlikes, reps);

  return (mean - var / 2.0);

}

/*-------------------------------------------------------*/
double
NETiv (int pop1, int pop2, double *P, int *NumAlleles, int reps)
{
  /* This function returns the current estimated average net nucleotide
     distance between the allele frequencies in populations 1 and 2.
     Here reps is the number of reps over which the P is an average, rather
     than the number of reps since the start of the program, as elsewhere */
  double sum, d1, d2,d12, norm;
  int loc,allele;
  norm = (double) reps;
  norm *= norm;
  sum=0.0;
  for (loc=0; loc<NUMLOCI;loc++) {
    d1=0.0;
    d2=0.0;
    d12=0.0;
    for (allele=0;allele<NumAlleles[loc];allele++) {
      d1+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop1,allele)]/norm;
      d2+=P[PPos(loc,pop2,allele)]*P[PPos(loc,pop2,allele)]/norm;
      d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
    }
    sum+= 0.5*(d1+d2)-d12;
  }
  return sum/NUMLOCI;
}


/*-------------------------------------------------------*/
double
GROSSiv (int pop1, int pop2, double *P, int *NumAlleles, int reps)
{
  /* This function returns the current estimated average gross nucleotide
     distance between the allele frequencies in populations 1 and 2.
     Here reps is the number of reps over which the P is an average, rather
     than the number of reps since the start of the program, as elsewhere */
  double sum, d12, norm;
  int loc,allele;
  norm = (double) reps;
  norm *= norm;
  sum=0.0;
  for (loc=0; loc<NUMLOCI;loc++)
  {
    d12=0.0;
    for (allele=0;allele<NumAlleles[loc];allele++)
    {
      d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
    }
    sum+= 1.0-d12;
  }
  return sum/NUMLOCI;
}

/*----------------------------------------------------*/
double
KLDiv (int pop1, int pop2, double *P, double *LogP, int *NumAlleles, int reps)
    /*This function returns the current (average) estimated
      Kullback-Leibler divergence between the allele frequencies in
      pops 1 and 2.  Here reps is the number of reps over which the P
      is an average, rather than the number of reps since the start of
      the program, as elsewhere. */
{
  double sum = 0.0;
  int allele, loc;

  for (loc = 0; loc < NUMLOCI; loc++)
    for (allele = 0; allele < NumAlleles[loc]; allele++)
      sum += ((double) P[PPos (loc, pop1, allele)] / reps)
          * log (P[PPos (loc, pop1, allele)] / P[PPos (loc, pop2, allele)]);

  return sum / NUMLOCI;

}

/*----------------------------------------------------*/
void
PrintNET (FILE * file, double *P, int *NumAlleles, int reps, int format)
    /*This function prints the current (average) estimated
      Net-nucleotide divergence between the allele frequencies in the
      different populations, in two formats.  Format 0 prints these in a
      row, averaging D_ij with D_ji, while Format 1 prints the full table */
{
  int i, j;
  /* int k;
   *  double div; */

  if (format == 0) {             /*print in line */
    for (i = 0; i < MAXPOPS - 1; i++) {
      for (j = i + 1; j < MAXPOPS; j++) {
        fprintf (file, "%1.3f ", NETiv(i,j,P,NumAlleles,reps));
      }
    }
  } else {
    /*print 2-D table */
    fprintf (file,
             "\nAllele-freq. divergence among pops (Net nucleotide distance)");
    if (reps > 1) {
      fprintf (file, ",\ncomputed using point estimates of P");
    }
    fprintf (file, ".\n\n");

    fprintf (file, "     ");
    for (j = 0; j < MAXPOPS; j++) {
      fprintf (file, "%2d      ", j + 1);
    }
    fprintf (file, "\n");

    for (i = 0; i < MAXPOPS; i++) {
      fprintf (file, "%2d   ", i + 1);
      for (j = 0; j < MAXPOPS; j++) {
        if (i == j) {
          fprintf (file, "   -    ");
        } else {
          fprintf (file, "%1.4f  ", NETiv (i, j, P, NumAlleles, reps));
        }
      }
      fprintf (file, "\n");
    }

    /* word change population -> cluster, William 03/27/07 */
    fprintf(file,"\nAverage distances (expected heterozygosity) between individuals in same cluster:\n");
    for (i=0;i<MAXPOPS;i++) {
      fprintf(file,"cluster %2d  : %1.4f \n",i+1,GROSSiv(i,i,P,NumAlleles,reps));
    }
    fprintf(file,"\n");
  }
}


/*----------------------------------------------------*/
void
PrintKLD (FILE * file, double *P, double *LogP, int *NumAlleles, int reps, int format)
    /*This function prints the current (average) estimated
      Kullback-Leibler divergence between the allele frequencies in the
      different populations, in two formats.  Format 0 prints these in a
      row, averaging D_ij with D_ji, while Format 1 prints the full table */
{
  int i, j;
  int k;
  double div;

  if (format == 0)              /*print in line */
  {
    for (i = 0; i < MAXPOPS - 1; i++)
      for (j = i + 1; j < MAXPOPS; j++)
      {
        /*printf("D%d%d = %1.2f, ",i,j,
          0.5*KLDiv(i,j,P,NumAlleles)+0.5*KLDiv(j,i,P,NumAlleles) ); */
        div = 0.5 * KLDiv (i, j, P,LogP, NumAlleles, reps) + 0.5 * KLDiv (j, i, P,LogP, NumAlleles, reps);
        fprintf (file, "%1.3f ", div);
        /*add extra spaces when i and j are large, to agree with spacing of banner */
        for (k = 2; k < ((int) log10 (i)) + ((int) log10 (j)); k++)
          fprintf (file, " ");
      }
  }
  else
    /*print 2-D table */
  {
    fprintf (file,
             "\nAllele-freq. divergence among pops (Kullback-Leibler distance)");
    if (reps > 1)
      fprintf (file, ",\ncomputed using point estimates of P");
    fprintf (file, ".\n\n");

    fprintf (file, "     ");
    for (j = 0; j < MAXPOPS; j++)
      fprintf (file, "%2d    ", j + 1);
    fprintf (file, "\n");

    for (i = 0; i < MAXPOPS; i++)
    {
      fprintf (file, "%2d  ", i + 1);
      for (j = 0; j < MAXPOPS; j++)
      {
        if (i == j)
          fprintf (file, "  -   ");
        else
          fprintf (file, "%2.2f  ", KLDiv (i, j, P,LogP, NumAlleles, reps));
      }
      fprintf (file, "\n");
    }
  }
}
/*---------------------------------------------------*/
void
UpdateSums (double *Q, double *QSum, int *Z, double *P, double *PSum,
            double *Fst, double *FstSum, int *NumAlleles,
            int *AncestDist, double *Epsilon, double *SumEpsilon,
            double *lambda, double *sumlambda, double *LocPrior,
            double *sumLocPrior, int LocPriorLen)
{
  int loc, ind, pop, allele, box, i;
  /*  int line; */

  for (pop=0;pop<MAXPOPS;pop++)
    sumlambda[pop] += lambda[pop];

  for (ind = 0; ind < NUMINDS; ind++)
    for (pop = 0; pop < MAXPOPS; pop++)
      QSum[QPos (ind, pop)] += Q[QPos (ind, pop)];


  for (loc = 0; loc < NUMLOCI; loc++)
    for (pop = 0; pop < MAXPOPS; pop++)
      for (allele = 0; allele < NumAlleles[loc]; allele++)
        PSum[PPos (loc, pop, allele)] += P[PPos (loc, pop, allele)];

  if (FREQSCORR)
  {
    for (pop = 0; pop < MAXPOPS; pop++)
      FstSum[pop] += Fst[pop];

    for (loc = 0; loc < NUMLOCI; loc++)
      for (allele = 0; allele < NumAlleles[loc]; allele++)
        SumEpsilon[EpsPos (loc, allele)] += Epsilon[EpsPos (loc, allele)];

  }

  if (ANCESTDIST)               /*store histogram of Q values for each individual */
  {
    for (ind = 0; ind < NUMINDS; ind++)
      for (pop = 0; pop < MAXPOPS; pop++)
      {
        box = ((int) (Q[QPos (ind, pop)] * ((double) NUMBOXES)));
        /*printf("%1.3f__%d  ",Q[QPos(ind,pop)],box); */
        if (box == NUMBOXES)
          box = NUMBOXES - 1;   /*ie, Q = 1.000 */
        AncestDist[AncestDistPos (ind, pop, box)]++;
      }
  }
  if (LOCPRIOR)
    for (i=0; i<LocPriorLen; i++)
      sumLocPrior[i] += LocPrior[i];

}
/*-----------------------------------------*/

/*-----------------------------------------*/
double
CalcLikeIndRecessive(int *Geno, int *PreGeno, double *AncVector, double *P, int ind, int *Recessive)
{
  /*returns log(likelihood) of Data for one individual:  log[ P(Data|p,q) ].
    This version is used only for diploid individuals when there is genotypic
    uncertainty (recessive model or inbreeding [in future])

    note code overlap with CalcLikeRecessive; make any updates to both functions

    Notice use of AncVector (Q for current individual only), not full Q)
  */

  double runningtotal = 1;
  double loglike = 0;
  double term, sum1, sum2;
  int allele1, allele2;
  int loc, pop;
  /*  int line; */

  for (loc = 0; loc < NUMLOCI; loc++) {
    allele1 = PreGeno[GenPos (ind, 0, loc)];
    allele2 = PreGeno[GenPos (ind, 1, loc)];

    if (allele1 != MISSING && allele2 != MISSING) { /* no missing data */
      sum1 = 0.0;
      sum2 = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        /* summing over possible genotypes and possible Z */
        sum1 += AncVector[pop] * P[PPos (loc, pop, allele1)];
        if (Recessive[loc] != MISSING && Recessive[loc] != allele1 && allele1 == allele2) { /* bug fixed 05072007 */
          sum2 += AncVector[pop] * (P[PPos (loc, pop, allele1)] + P[PPos(loc, pop, Recessive[loc])]);
        } else {
          sum2 += AncVector[pop] * P[PPos (loc, pop, allele2)];
        }
      }
      term = sum1 * sum2;
    } else if (allele1!=MISSING)  { /* one allele missing */
      term=0.0;
      for (pop=0;pop<MAXPOPS;pop++) {
        term += AncVector[pop] * P[PPos (loc, pop, allele1)];
      }
    } else if (allele2!=MISSING) {
      term=0.0;
      for (pop=0;pop<MAXPOPS;pop++) {
        term += AncVector[pop] * P[PPos (loc, pop, allele2)];
      }
    } else {
      term=1.0;  /* no data */
    }

    runningtotal *= term;

    if (runningtotal < UNDERFLO) {        /*this is to avoid having to take logs all the time */
      loglike += log (runningtotal);
      runningtotal = 1;
    }
  }

  loglike += log (runningtotal);
  return loglike;
}

/*-----------------------------------------*/
double CalcLikeInd (int *Geno, int *PreGeno, double *AncVector, double *P, int ind, int *Recessive)
{
  /*returns log(likelihood) of Data for one individual:  log[ P(Data|p,q) ]
    See notes 19 March 99 */

  /*This likelihood is used for the metropolis update of Q*/

  /* when there is genotypic ambiguity (eg in the recessive model) the likelihood
     is computed in one of two ways: (1) for diploids it sums over possible
     genotypes, and (2) for other ploidy it is computed based on the current
     imputed genotypes

     note code overlap with CalcLike, CalcLikeIndRecessive;
     make any updates to all functions

     Notice use of AncVector (Q for current individual only), not full Q)
  */

  double runningtotal = 1;
  double loglike = 0;
  double term;
  int allele;
  int line, loc, pop;
  double sqrtunder = sqrt (UNDERFLO);

  if (LINES==2 && RECESSIVEALLELES) {
    loglike = CalcLikeIndRecessive(Geno, PreGeno, AncVector, P, ind, Recessive);
  } else {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)];
        if (allele != MISSING) {
          term = 0.0;
          for (pop = 0; pop < MAXPOPS; pop++) {
            term += AncVector[pop] * P[PPos (loc, pop, allele)];
          }
          
          if (term > sqrtunder) {
            runningtotal *= term;
          } else {
            runningtotal *= sqrtunder;
          }

          if (runningtotal < sqrtunder) { /*this is to avoid having to take logs all the time */
            if (runningtotal == 0.0) {
              printf ("Error in CalcLikeInd\n");
            }
            loglike += log (runningtotal);
            runningtotal = 1;
          }
        }
      }
    }
    loglike += log (runningtotal);
  }
  return loglike;
}

/*-----------------------------------------*/

double CalcLikeRecessive(int *Geno, int *PreGeno, double *Q, double *P, int *Recessive, double *sumindlike, double *indlike_norm)
    /* calculate log likelihood for recessive model.  Only coded for diploids!
     *
     * note code overlap with CalcLikeInd, CalcLikeIndRecessive, CalcLike;
     * make any updates to all functions */

    /*7/12/07: Melissa modified so that it calls CalcLikeIndRecessive
      and there is no more code overlap.  This is more efficient and better
      coding practice, and also helps with calculation of DIC (which
      requires storage of individual likelihoods)*/

{


  double loglike = 0;
  int ind, pop;
  double *AncVector = malloc(MAXPOPS*sizeof(double));
  double indlike;
  if (AncVector==NULL) {
    printf("Error allocating Ancvetor in CalcLikeRecessive: out of memory?\n");
    exit(-1);
  }

  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop=0; pop<MAXPOPS; pop++)
      AncVector[pop] = Q[QPos(ind, pop)];
    indlike = CalcLikeIndRecessive(Geno, PreGeno, AncVector, P, ind,
                                   Recessive);
    if (sumindlike!=NULL) sumindlike[ind] += exp(indlike-indlike_norm[ind]);
    loglike += indlike;
  }
  free(AncVector);
  return loglike;
}
/*-----------------------------------------*/
double CalcLike (int *Geno, int *PreGeno, double *Q, double *P, int *Recessive,
                 double *sumindlike, double *indlike_norm)
{
  /*returns log(likelihood) of Data:  log[ P(Data|p,q) ]
    See notes 19 March 99 */

  /* note code overlap with CalcLikeInd, CalcLikeIndRecessive, CalcLikeRecessive;
   * make any updates to all functions */
  /*Melissa modified 7/12/07 so it calls CalcLikeInd rather than repeats
    the same code.  This helps with DIC calculation...*/

  double loglike = 0, indlike, *AncVector;
  int ind, pop;

  if (LINES == 2 && RECESSIVEALLELES) {
    loglike = CalcLikeRecessive(Geno, PreGeno, Q, P, Recessive, sumindlike, indlike_norm);
  } else {
    AncVector = malloc(MAXPOPS*sizeof(double));
    if (AncVector==NULL) {
      printf("Error allocating Ancvetor in CalcLikeRecessive: out of memory?\n");
      exit(-1);
    }
    for (ind = 0; ind < NUMINDS; ind++) {
      for (pop=0; pop<MAXPOPS; pop++) {
        AncVector[pop] = Q[QPos(ind, pop)];
      }
      indlike = CalcLikeInd(Geno, PreGeno, AncVector, P, ind, Recessive);
      if (sumindlike!=NULL) {
        if (indlike_norm[ind]==0.0) {
          indlike_norm[ind] = indlike;
        }
        sumindlike[ind] += exp(indlike-indlike_norm[ind]);
      }
      loglike += indlike;
    }
    free(AncVector);
  }
  return loglike;
}
/*====================================================*/
/*---------------------------------------------------*/
void
PrintMainParams (FILE * file, int rep, int argc, char *argv[])
{
  int i;
  /*print values of most important parameters */

  fprintf (file, "\n");
  if (argc > 1)
  {
    fprintf (file, "Command line arguments:   ");
    for (i = 0; i < argc; i++)
      fprintf (file, "%s ", argv[i]);
    fprintf (file, "\n");
  }
  fprintf (file, "Input File:    %s\n", DATAFILE);

  if (file == stdout)
    fprintf (file, "Output File:   %s_f\n", OUTFILE);
  fprintf (file, "\n");
  fprintf (file, "Run parameters:\n");
  fprintf (file, "   %d individuals\n", NUMINDS);
  fprintf (file, "   %d loci\n", NUMLOCI);
  fprintf (file, "   %d populations assumed\n", MAXPOPS);
  fprintf (file, "   %d Burn-in period\n", BURNIN);
  fprintf (file, "   %d Reps\n", rep - BURNIN);

  if (USEPOPINFO)
  {
    fprintf (file, "USEPOPINFO turned on\n");
    fprintf (file, "MIGRPRIOR = %1.4f\n", MIGRPRIOR);
  }
  if (RECESSIVEALLELES)
    fprintf(file,"RECESSIVE ALLELES model used\n");
  if (NOADMIX)
    fprintf (file, "NO ADMIXTURE model assumed\n");
  if (STARTATPOPINFO)
    fprintf (file, "STARTATPOPINFO turned on\n");
  if (LOCPRIOR)
    fprintf (file, "LOCPRIOR model used\n");
  if (!(RANDOMIZE))
    fprintf (file, "RANDOMIZE turned off\n");
  fprintf (file, "\n");

}
/*----------------------------------------------------*/
void
PrintAncestDist (FILE * file, int *AncestDist, int ind, int rep)
    /*print credible region for each q in each population */
{
  int pop, box, low;
  double sum;
  /*print the whole histogram */
  /*fprintf(file,"\n");
    for (pop=0; pop<MAXPOPS; pop++)
    {
    fprintf(file,"%d: ",pop);
    for (box=0; box<NUMBOXES; box++)
    fprintf(file,"%d ",AncestDist[AncestDistPos(ind,pop,box)]);
    fprintf(file,"\n");
    } */

  fprintf (file, "    ");
  for (pop = 0; pop < MAXPOPS; pop++)
  {
    sum = 0;
    low = 0;
    for (box = 0; box < NUMBOXES; box++)
    {
      sum += AncestDist[AncestDistPos (ind, pop, box)];
      if ((low == 0) && (sum > (int) ((rep - BURNIN) * (1.0 - ANCESTPINT) / 2.0)))
      {
        /*printf("lo: sum = %d, value = %d\n",
          (int) (rep-BURNIN)*(1.0-ANCESTPINT)/2.0); */
        fprintf (file, "(%1.3f,", (double) (box + 0.0) / NUMBOXES);
        low = 1;
      }
      if (sum > ((rep - BURNIN) * (1 + ANCESTPINT) / 2.0))
      {
        /*{printf("hi: sum = %d, value = %d\n",
          (int) (rep-BURNIN)*(1.0-ANCESTPINT)/2.0); */
        fprintf (file, "%1.3f) ", (double) (box + 1.0) / NUMBOXES);
        break;
      }
    }
  }
}
/*----------------------------------------------------*/
void
PrintGeogQ (FILE * file, int ind, int pop, double *UsePopProbs, int rep)
{
  int homepop, gen;

  homepop = pop - 1;
  fprintf (file, "%1.3f | ",
           (double) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN));
  for (pop = 0; pop < MAXPOPS; pop++)
    if (pop != homepop)
    {
      fprintf (file, "Pop %d: ", pop + 1);
      for (gen = 0; gen < GENSBACK + 1; gen++)
        fprintf (file, "%1.3f ",
                 (double) UsePopProbs[UsePPrPos (ind, pop, gen)] / (rep - BURNIN));
      fprintf (file, " | ");
    }

  fprintf (file, " ");
  if ((double) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.1)
    fprintf (file, "*");
  if ((double) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.3)
    fprintf (file, "*");
  if ((double) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.5)
    fprintf (file, "*");

}
/*----------------------------------------------------*/
void
PrintQ (FILE * file, int *Geno, int rep, double *QSum, struct IND *Individual,
        int *AncestDist, double *UsePopProbs,double *sumR)
    /*print summary of Q to file */
{
  int ind, pop;
  int missing;                  /*(% missing data for this individual) */
  /*  int sum, low, box; */
  int MissingInd (int *Geno, int ind);


  fprintf (file, "\n");
  fprintf (file, "Inferred ancestry of individuals:\n");
  if (USEPOPINFO)
    fprintf (file, "Probability of being from assumed population | prob of other pops\n");
  fprintf (file, "    ");
  if (LABEL)
    fprintf (file, " %8s", "Label");
  fprintf (file, " (%%Miss) ");
  if (POPDATA)
    fprintf (file, "Pop");
  if (!USEPOPINFO)
    fprintf(file, ":  ");
  if (LINKAGE && INDIVIDUALR)
    fprintf(file, " Ind's r ");
  if (!(USEPOPINFO))
  {
    fprintf (file, "Inferred clusters");
    if (ANCESTDIST)
      fprintf (file, " (and %d%c probability intervals)",
               (int) (100.0 * ANCESTPINT), '%');
  }
  fprintf (file, "\n");

  /*
    if (LABEL) fprintf(file,"    %8s ","");
    fprintf(file,"         ");
    if (POPDATA) fprintf(file,"   ");
    fprintf(file,":   ");
    for (pop=0; pop<MAXPOPS; pop++)
    fprintf(file,"%2d    ",pop+1);
    fprintf(file,"\n"); */

  for (ind = 0; ind < NUMINDS; ind++)
  {
    missing = (int) (100 * MissingInd (Geno, ind)) / (LINES * NUMLOCI);

    fprintf (file, "%3d ", ind + 1);
    if (LABEL)
      fprintf (file, "%8s ", Individual[ind].Label);
    if (missing < 10)
      fprintf (file, " ");
    fprintf (file, "  (%d) ", missing);
    fprintf (file, "  ");
    if (POPDATA)
      fprintf (file, "%2d ", Individual[ind].Population);
    fprintf (file, ":  ");
    if (LINKAGE && INDIVIDUALR)
      fprintf(file, " %1.4f  ",(double)sumR[ind]/(rep-BURNIN));
    if ((USEPOPINFO) && (Individual[ind].PopFlag))
      PrintGeogQ (file, ind, Individual[ind].Population, UsePopProbs, rep);

    else
    {
      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, "%1.3f ", (double) QSum[QPos (ind, pop)] / (rep - BURNIN));

      if (ANCESTDIST)   /*Print the credible intervals for ancestry coeffs */
        PrintAncestDist (file, AncestDist, ind, rep);
    }
    fprintf (file, "\n");
  }
}
/*-----------------------------------------------------*/
void
PrintQFile (int rep, double *QSum, struct IND *Individual, double *UsePopProbs)
    /*Produce a file that contains only the label and popinfo (if any) plus Q-hat */
{
  char outname[STRLEN + 20];
  FILE *QHat;
  int ind, pop;
  /*double *QProbs; *//*[MAXPOPS] */

  sprintf (outname, "%s_q", OUTFILE);
  QHat = fopen (outname, "w");
  if (QHat == NULL)
    printf ("WARNING: Unable to open output file %s.\n", outname);

  /*QProbs = calloc(MAXPOPS,sizeof(double));
    if (QProbs==NULL) printf("Warning: unable to assign memory in PrintQFile\n"); */

  for (ind = 0; ind < NUMINDS; ind++)
  {
    if (LABEL)
      fprintf (QHat, "%12s ", Individual[ind].Label);
    else
      fprintf (QHat, "%4d ", ind + 1);
    if (POPDATA)
      fprintf (QHat, "%2d ", Individual[ind].Population);
    /*if ((USEPOPINFO)&&(Individual[ind].PopFlag))
      QFromUsePop(ind,QProbs,UsePopProbs,Individual[ind].Population-1);
      else
      for (pop=0; pop<MAXPOPS; pop++)
      QProbs = (double) QSum[QPos(ind,pop)]/(rep-BURNIN); */
    for (pop = 0; pop < MAXPOPS; pop++)
      fprintf (QHat, "%1.4f ", (double) QSum[QPos (ind, pop)] / (rep - BURNIN));
    fprintf (QHat, "\n");
  }
  fclose (QHat);
  /*free(QProbs); */
}
/*-----------------------------------------------------*/
void
PrintP (FILE * file, int rep, int *Geno, double *PSum, int *Translation, int *NumAlleles,
        double *SumEpsilon, char *Markername)
    /*print summary of P to file */
{
  int loc, pop, allele;
  int MissingLoc (int *Geno, int loc);


  fprintf (file, "\n\nEstimated Allele Frequencies in each cluster\n");
  if (FREQSCORR)
    fprintf (file, "First column gives estimated ancestral frequencies\n");
  fprintf (file, "\n\n");

  for (loc = 0; loc < NUMLOCI; loc++)
  {
    fprintf (file, "Locus %d : ", loc + 1);
    if (MARKERNAMES) PrintGeneName(file, loc, Markername);
    fprintf (file, "\n");
    fprintf (file, "%d alleles\n", NumAlleles[loc]);
    fprintf (file, "%2.1f%c missing data\n",
             (double) 100.0 * MissingLoc (Geno, loc) / (LINES * NUMINDS), '%');  /* JKP changed 2 to LINES (4/06/09) */
    for (allele = 0; allele < NumAlleles[loc]; allele++)
    {
      /*This fine piece of programming uses the allele 29697 to represent the
        recessive allele in the event that it is not observed in the data; used
        only when the RECESSIVEALLELES model is turned on. This is coded in datain.c*/
      if (RECESSIVEALLELES
          && Translation[TransPos (loc, allele)] == 29697)
        fprintf (file, "Null   ");

      else fprintf (file, "%4d   ", Translation[TransPos (loc, allele)]);

      if (FREQSCORR)
        fprintf (file, "(%1.3f) ",
                 (double) SumEpsilon[EpsPos (loc, allele)] / (rep - BURNIN));

      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, "%1.3f ", (double) PSum[PPos (loc, pop, allele)] / (rep - BURNIN));
      fprintf (file, "\n");
    }
    fprintf (file, "\n");
  }
}

int get_location_num(int loc, struct IND *individual) {
  int ind;
  for (ind=0; ind<NUMINDS; ind++)
    if (individual[ind].myloc==loc) return individual[ind].Location;
  printf("error in get_location_num: can't find location %i\n", loc);
  exit(-1);
  return -1;
}

/*-----------------------------------------------------*/
void
PrintSums (FILE * file, int rep, double sumlikes,
           double sumsqlikes, double *FstSum, double *sumAlpha,
           double *sumlambda, double *sumR,double *varR,
           struct IND *Individual, double *sumLocPrior, int LocPriorLen,
           double DIC)
    /*print current value of some averages to file */
{
  int ind, pop, locnum, loc;
  double sumrecs = 0.0;

  fprintf (file, "--------------------------------------------\n");
  if (COMPUTEPROB)
  {
    if (rep - BURNIN > 2)
    {

      fprintf (file, "Estimated Ln Prob of Data   = %1.1f\n",
               EstLogProb (sumlikes, sumsqlikes, rep - BURNIN));
      fprintf (file, "Mean value of ln likelihood = %1.1f\n",
               (double) sumlikes / (rep - BURNIN));
      fprintf (file, "Variance of ln likelihood   = %1.1f\n",
               SampleVar (sumsqlikes, sumlikes, (rep - BURNIN)));
      /*          fprintf (file, "DIC = %.3f\n", DIC); */
    }

    else
      fprintf (file, "Mean value of ln likelihood = %1.1f\n",
               (double) sumlikes / (rep - BURNIN));     /*print this line in either case */
  }

  if (((!(NOADMIX)) && (!(NOALPHA)) && (MAXPOPS > 1)))
  {
    if (POPALPHAS)
    {
      fprintf (file, "\n");
      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, "Mean value of alpha_%d       = %1.4f\n", pop + 1,
                 (double) sumAlpha[pop] / (rep - BURNIN));
    }
    else
      fprintf (file, "Mean value of alpha         = %1.4f\n",
               (double) sumAlpha[0] / (rep - BURNIN));

    if (LOCPRIOR) {
      fprintf(file, "\nMean value of alpha_local for each location:\n");
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        locnum = get_location_num(loc, Individual);
        fprintf(file, "\tlocation %2i:", locnum);
        for (pop=0; pop<MAXPOPS; pop++)
          fprintf(file, "  %1.4f", sumAlpha[AlphaPos(loc, pop)]/(rep - BURNIN));
        fprintf(file, "\n");
      }
    }
  }

  if (INFERLAMBDA)
  {
    if (POPSPECIFICLAMBDA)
      for (pop=0;pop<MAXPOPS;pop++)
        fprintf (file, "\nMean value of lambda%d       = %1.4f\n",pop+1,
                 (double) sumlambda[pop] / (rep - BURNIN));
    else

      fprintf (file, "\nMean value of lambda        = %1.4f\n",
               (double) sumlambda[0] / (rep - BURNIN));
  }
  if (LINKAGE)
  {
    if (!INDIVIDUALR)
    {
      fprintf (file, "Mean value of r              = %1.4f\n",
               (double) sumR[0] / (rep - BURNIN));
      fprintf (file, "Standard deviation of r    = %1.4f\n",
               sqrt((double) varR[0] / (double) (rep - BURNIN) - sumR[0] * sumR[0] / (double) ((rep - BURNIN) * (rep - BURNIN))));
    }
    else
    {
      for (ind = 0; ind < NUMINDS; ind++)
        sumrecs += sumR[ind];
      fprintf (file, "Mean value of r           = %1.6f\n", (double) sumrecs / NUMINDS / (rep - BURNIN));

    }
  }

  if (FREQSCORR)
  {
    if (ONEFST)
      fprintf (file, "Mean value of Fst           = %1.4f\n",
               (double) FstSum[0] / (rep - BURNIN));
    else
    {
      fprintf (file, "\n");
      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, "Mean value of Fst_%d         = %1.4f\n", pop + 1,
                 (double) FstSum[pop] / (rep - BURNIN));
    }
  }
  else
    fprintf (file, "Allele frequencies uncorrelated\n");

  if (PFROMPOPFLAGONLY)
    fprintf (file, "Allele frequencies updated using individuals with POPFLAG=1 ONLY.\n");

  fprintf (file, "\n");

  if (LOCPRIOR) {
    fprintf(file, "Mean value of r = %1.4f\n", sumLocPrior[0]/(rep-BURNIN));
    if (NOADMIX) {
      fprintf(file, "Mean value of nu = ");
      for (pop=0; pop<MAXPOPS; pop++)
        fprintf(file, " %1.4f", sumLocPrior[LocPriorPos(NUMLOCATIONS, pop)]/(rep-BURNIN));
      fprintf(file, "\n");
      fprintf(file, "Mean value of gamma for each location:\n");
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        locnum = get_location_num(loc, Individual);
        fprintf(file, "\tlocation %2i:", locnum);
        for (pop=0; pop<MAXPOPS; pop++)
          fprintf(file, "  %1.4f", sumLocPrior[LocPriorPos(loc, pop)]/(rep-BURNIN));
        fprintf(file, "\n");
      }
    }
  }
}


/*-----------------------------------------------------*/
int
MissingLoc (int *Geno, int loc)
    /*return the number of missing alleles at a locus */
{
  int ind, line;
  int sofar = 0;

  for (ind = 0; ind < NUMINDS; ind++)
    for (line = 0; line < LINES; line++)
    {
      if (Geno[GenPos (ind, line, loc)] == MISSING)
        sofar++;
    }

  return sofar;
}
/*-----------------------------------------------------*/
int
MissingInd (int *Geno, int ind)
    /*return the number of missing alleles in an individual */
{
  int loc, line;
  int sofar = 0;

  for (loc = 0; loc < NUMLOCI; loc++)
    for (line = 0; line < LINES; line++)
    {
      if (Geno[GenPos (ind, line, loc)] == MISSING)
        sofar++;
    }

  return sofar;
}
/*------------------------------------------------------------------------*/
void
PrintSiteBySite (FILE * file, double *SiteBySiteSum, int rep, char *Markername, double *PSum, int *NumAlleles, int *Geno, int *Translation, double *SumEpsilon,struct IND *Individual)
{
  int ind,pop,pop2,line,loc;
  fprintf (file, "\n\n");
  /*
    fprintf(file, "Here are the site by site outputs for each individual.\n");
    fprintf(file, "See user guide for details \n"); */
  for (ind = 0; ind < NUMINDS; ind++)
  {
    for (loc = 0; loc < NUMLOCI; loc++)
    {
      fprintf (file, "%d %d ", ind+1,loc+1);
      /*
        if (MARKERNAMES) PrintGeneName(file, loc, Markername);
        if (LABEL)
        fprintf (file, "%8s ", Individual[ind].Label);
      */
      if (PHASED || !LINKAGE)
        for (line = 0; line < LINES; line++)
          for (pop = 0; pop < MAXPOPS; pop++)
            fprintf (file, "%1.3f ", SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] / (double) (rep - BURNIN));
      else
        for (pop = 0; pop < MAXPOPS; pop++)
          for (pop2=0;pop2<MAXPOPS;pop2++)
            fprintf(file,"%1.3f ", SiteBySiteSum[DiploidSiteBySiteSumPos (ind,pop2,loc,pop)]/(double)(rep-BURNIN));
      fprintf (file, "\n");
    }
    fprintf (file, "\n");
  }
}
/*----------------------------------------------------*/
void PrintGeneName(FILE * file, int loc, char *Markername)
{
  int i;
  for (i=0; i<GENELEN; i++)
  {
    if (Markername[MarkernamePos(loc,i)] != '\0')
      fprintf(file,"%c",Markername[MarkernamePos(loc,i)]);
    else
    {
      if (i==0) fprintf(file,"XXX");
      fprintf(file," ");
      break;
    }
  }
}
/*----------------------------------------------------*/
int EqualGeneNames(int loc1,int loc2,char *Markername)
    /*returns 1 if the gene names are the same, otherwise 0*/
{
  int i;

  for (i=0; i<GENELEN; i++)
  {
    if (Markername[MarkernamePos(loc1,i)] != Markername[MarkernamePos(loc2,i)])
      return 0;
    if (Markername[MarkernamePos(loc1,i)] == '\0')
      return 1;
  }
  return 1;
}
/*----------------------------------------------------*/
void
PrintMembership (FILE * file, double *QSum, struct IND *Individual)
    /*Print summary of relationship between given populations,
      and cluster populations.  Requires POPDATA. An earlier, more
      complicated version of this is stored in backup.c */
{
  double *sumvals;              /*this array stores the sum of QSum for each of the
                                  the cluster populations 0..MAXPOPS-1, for individuals
                                  who are designated as coming from particular population. */
  double rowsum;                /*sum of the values in the array sumvals */
  int ind;
  int minpop;                   /*value of smallest (largest) population number */
  int maxpop;
  int pop, givenpop;
  int numfrompop;               /*number of individuals from each given pop */

  sumvals = calloc (MAXPOPS, sizeof (double));
  if (sumvals == NULL)
    fprintf (file, "Error assigning memory in function PrintMembership\n");
  else
  {

    if (POPDATA)
    {
      minpop = Individual[0].Population;        /*figure out min and max population names */
      maxpop = Individual[0].Population;
      for (ind = 1; ind < NUMINDS; ind++)
      {
        if (Individual[ind].Population < minpop)
          minpop = Individual[ind].Population;
        if (Individual[ind].Population > maxpop)
          maxpop = Individual[ind].Population;
      }

      if (sumvals == NULL)
        fprintf (file, "Error assigning memory in function PrintMembership\n");
      else
      {
        fprintf (file, "\n--------------------------------------------\n");
        fprintf (file, "Proportion of membership of each pre-defined\n");
        fprintf (file, " population in each of the %d clusters\n\n", MAXPOPS);

        fprintf (file, "Given    Inferred Clusters");
        for (pop = 3; pop < MAXPOPS; pop++)
          fprintf (file, "       ");
        fprintf (file, "       Number of\n");

        fprintf (file, " Pop    ");
        for (pop = 0; pop < MAXPOPS; pop++)
          fprintf (file, "  %2d   ", pop + 1);
        for (pop = MAXPOPS; pop < 3; pop++)
          fprintf (file, "       ");
        fprintf (file, "   Individuals\n\n");

        for (givenpop = minpop; givenpop <= maxpop; givenpop++)
        {
          for (pop = 0; pop < MAXPOPS; pop++)
            sumvals[pop] = 0.0;

          numfrompop = 0;
          for (ind = 0; ind < NUMINDS; ind++)
          {
            if (givenpop == Individual[ind].Population)
            {
              numfrompop++;
              for (pop = 0; pop < MAXPOPS; pop++)
                sumvals[pop] += QSum[QPos (ind, pop)];
            }
          }
          rowsum = 0.0;
          for (pop = 0; pop < MAXPOPS; pop++)
            rowsum += sumvals[pop];

          if (rowsum > 0.0)
          {
            fprintf (file, "%3d:     ", givenpop);
            for (pop = 0; pop < MAXPOPS; pop++)
              fprintf (file, "%1.3f  ", sumvals[pop] / rowsum);
            for (pop = MAXPOPS; pop < 3; pop++) /*number of individuals */
              fprintf (file, "       ");
            fprintf (file, "    %3d\n", numfrompop);
          }

        }
      }
    }
    else
      /* no popdata */
    {
      for (pop = 0; pop < MAXPOPS; pop++)
        sumvals[pop] = 0.0;

      for (ind = 0; ind < NUMINDS; ind++)
        for (pop = 0; pop < MAXPOPS; pop++)
          sumvals[pop] += QSum[QPos (ind, pop)];

      rowsum = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++)
        rowsum += sumvals[pop];

      fprintf (file, "\n--------------------------------------------\n");
      fprintf (file, "Overall proportion of membership of the\n");
      fprintf (file, "sample in each of the %d clusters\n\n", MAXPOPS);

      fprintf (file, "Inferred Clusters\n");
      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, " %2d    ", pop + 1);
      fprintf (file, "\n");
      for (pop = 0; pop < MAXPOPS; pop++)
        fprintf (file, "%1.3f  ", sumvals[pop] / rowsum);
      fprintf (file, "\n\n");


    }

    free (sumvals);
    fprintf (file, "--------------------------------------------\n");
  }
}

double CalcDIC(int rep, double sumlikes, double* sumindlikes,
               double* indlikes_norm) {
  double sumind=0.0, dic;
  int ind;
  for (ind=0; ind<NUMINDS; ind++)
    sumind += log(sumindlikes[ind]/(rep-BURNIN))+indlikes_norm[ind];
  dic = -4.0*sumlikes/(rep-BURNIN)+2.0*sumind;
  /*  return -4.0*sumlikes/(rep-BURNIN) + 2.0*log(sumindlikes)/(rep-BURNIN); */
  return dic;
}


/*====================================================*/
/*Melissa modified 7/12/07 to incorporate poppriors and DIC*/
void
OutPutResults (int *Geno, int rep, int savefreq,
               struct IND *Individual,
               double *PSum, double *QSum, double *SiteBySiteSum,
               double *FstSum, int *AncestDist, double *UsePopProbs,
               double sumlikes, double sumsqlikes, double *sumAlpha,
               double *sumR, double *varR,
               int *NumAlleles, int *Translation, int final,
               char *Markername, double *R, double *SumEpsilon,
               double *lambda, double *sumlambda,
               double *sumLocPrior, int LocPriorLen,
               double *sumindlikes, double *indlikes_norm,
               int argc, char *argv[])

    /*final indicates that the program is terminating.  Stuff gets
      printed to the screen at this stage. */
{
  /*print a bunch of stuff to file, and a subset of this to the screen:
    P, Q, Fsts, Net distances, likelihood results, all parameters, amount of
    missing data for each individual... */

  char outname[STRLEN + 20];
  char outname2[STRLEN + 20];

  FILE *RESULTS;
  FILE *RESULTS2;
  /*  int outputoption; */
  double DIC = 0.0;
  /*  DIC = CalcDIC(rep, sumlikes, sumindlikes, indlikes_norm); */
  if (final) {
    sprintf (outname, "%s_f", OUTFILE);
  } else {
    sprintf (outname, "%s_%d", OUTFILE, (rep - BURNIN) / savefreq);
  }
  
  RESULTS = fopen (outname, "w");
  if (RESULTS == NULL) {
    printf ("WARNING: Unable to open output file %s.\n", outname);
  } else {
    Welcome (RESULTS);
    if (final) {
      printf ("\nMCMC completed\n");
    }
    PrintMainParams (RESULTS, rep, argc, argv);
    PrintMembership (RESULTS, QSum, Individual);
    PrintNET (RESULTS, PSum, NumAlleles, rep - BURNIN, 1);
    /*if (final) PrintNET(stdout,PSum,NumAlleles,rep-BURNIN,1); */
    PrintSums (RESULTS, rep, sumlikes, sumsqlikes, FstSum, sumAlpha, sumlambda, sumR,varR, Individual, sumLocPrior, LocPriorLen, DIC);

    PrintQ (RESULTS, Geno, rep, QSum, Individual, AncestDist, UsePopProbs,sumR);
    PrintP (RESULTS, rep, Geno, PSum, Translation, NumAlleles, SumEpsilon, Markername);
    if (final) {
      PrintQ (stdout, Geno, rep, QSum, Individual, AncestDist, UsePopProbs,sumR);
      if (PRINTQHAT) {
        PrintQFile (rep, QSum, Individual, UsePopProbs);
      }
      PrintMainParams (stdout, rep, argc, argv);
      PrintSums (stdout, rep, sumlikes, sumsqlikes, FstSum, sumAlpha, sumlambda, sumR,varR, Individual, sumLocPrior, LocPriorLen, DIC);
      PrintMembership (stdout, QSum, Individual);
    }

    PrintAllParams (RESULTS);
    if (final) {
      printf ("Final results printed to file %s\n\n", outname);
    }
  }
  fclose (RESULTS);
    
  if ((final) && (SITEBYSITE)) {
    sprintf (outname2, "%s_ss", OUTFILE);
    RESULTS2 = fopen (outname2, "w");
    if (RESULTS2 == NULL) {
      printf ("WARNING: Unable to open sitebysite file %s.\n", outname2);
    } else {
      PrintSiteBySite (RESULTS2, SiteBySiteSum, rep, Markername, PSum, NumAlleles, Geno, Translation,SumEpsilon,Individual);
      printf("sitebysite results printed to file %s\n\n",outname2);
    }
  }
  
}
