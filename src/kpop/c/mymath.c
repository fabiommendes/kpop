
/*

VERSION 4-4-00

 */

#include "mymath.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "ran.h"

 
/*                                 RETURNS:
Square(x)                          square
SampleVar(sumsq,sum,numreps)       sample variance
PickAnOption                       A random number between 1 and n, according
                                       to a list of probabilities.
LDirichletProb                     Probability of a vector given Dirichlet prior
LGammaDistProb                     Probability of a draw from a Gamma distribution
FindAveLogs                        Estimated average of a series of log numbers.
RandomOrder                        Return a list of integers from 0 to L-1 in 
                                        random order.
double Factorial(int n)            Return n!, as a double

double ChiSq(int *list1,int len1,int *list2,int len2,int mincount,int missing,int *df)
       This function takes two lists of integers (allele scores),
and computes a chi^sq for the equality of frequencies.  list1 and
list2 are the two lists, of length len1, and len2, respectively.
Missing data are marked by the value "missing".  Classes with fewer
than "mincount" observations are pooled.  The function returns the
chisq value, and a pointer to the resulting degrees of freedom (df).

extern double LDirichletProb(double prior[],double post[],int length)
extern double Square(double x);
extern double SampleVar(double sum,double sumsq,double numreps);
extern int PickAnOption(int total,double sum,double Probs[]);
extern double FindAveLogs(double *logmax,double *sum, double lognext,int rep);
extern void RandomOrder(int list[],int length);
extern double Factorial(int n);
*/

#define NUMUNDERFLOW 1e-200           /*maximum resolution of doubles*/

/*---------------------------------*/
double Square(double x)
{
  return x*x;
}
/*---------------------------------*/
double SampleVar(double sumsq, double sum, long num)
/*returns the value of the sample variance, computed using the sum of
x_i^2 (sumsq), the sum of x_i (sum), and the sample size (num)*/
{
    if (num==1)
	{
	    printf("Division by zero in variance calculation.\n");
	    return 0;
	}
    else return (sumsq - sum*sum/num)/(num-1);
}

/*---------------------------------*/
double SD(double sumsq, double sum, long num)
{
/*returns the value of the square root of the sample variance, 
computed using the sum of x_i^2 (sumsq), the sum of x_i (sum), 
and the sample size (num)*/

  double var; 

  var = SampleVar(sumsq,sum,num);
  /*printf("var: %1.3f\n",var);*/

return  exp(log(var)/2) ;
  
}
/*---------------------------------*/
int PickAnOption(int total,double sum,double Probs[])
{
  /*Returns a random number between 0 and n-1, according to a list of
    probabilities.  The function takes a (possibly) unnormalized
    vector of probabilities of the form (p1, p2, p3,....,p_n).  There
    are "total" possible options, and sum is the sum of the
    probabilities.  This comes up in the Gibbs sampler context.*/

  int option;
  double random;
  double sumsofar = 0.0;

  random = RandomReal(0,sum);     /*Get uniform random real in this range*/
  for (option=0; option<total; option++) /*Figure out which popn this is*/
    {
      sumsofar += Probs[option];
	if (random <= sumsofar) break;
    }

  return option;

}
/*---------------------------------*/
double LDirichletProb(double prior[],double post[],int length)
/*returns the log probability of a vector "post" of length "length", 
  given a Dirichlet process with prior "prior". */
{
  double sumprior = 0.0;
  double logsum;
  int i;

  for (i=0; i<length; i++)
    sumprior += prior[i];

  logsum = mylgamma(sumprior);
  for (i=0; i<length; i++)
    logsum += (prior[i]-1.0)*log(post[i]) - mylgamma(prior[i]);
  
  return logsum;
}
/*---------------------------------*/
double LGammaDistProb(double alpha,double beta, double y)
/*returns the log probability of a gamma-distributed random
variable "y", with parameters alpha and beta, where the mean
is alpha*beta, and the variance is alpha*beta*beta*/
{

  double logsum;

  logsum = -1*mylgamma(alpha) - alpha*log(beta) + (alpha-1)*log(y) - y/beta;

  return logsum;

}
/*----------------------------------*/
double FindAveLogs(double *logmax,double *sum, double lognext,int rep)
/*This function is for use in estimating the mean of a set of really
large or small numbers, where it is necessary to use logs to avoid
numerical problems.  This function returns the log of the current
estimate of the mean.

logmax is the maximum value observed so far; 
lognext is the current value;
max*sum/reps is the current estimate of the mean;
rep is the number of reps so far.  

Mean = Pmax*Sum/rep, where Sum = sum( P_i/Pmax ).
Log Mean = logmax -log(rep) + log(sum).  */
{

  if (lognext <= (*logmax))
    {
      if (exp(lognext-(*logmax)) > NUMUNDERFLOW)
	*sum += exp(lognext-(*logmax));
    }
  
  else
    {
      if (*sum*exp((*logmax)-lognext) < NUMUNDERFLOW)
	*sum = 1.0;
      else
	*sum = 1.0 + (*sum)*exp((*logmax)-lognext);
      *logmax = lognext;
      /*printf("New logmax = %5.1f\n",*logmax); */

    }

  /*printf("max = %5.0f, next = %5.0f, sum = %f\n ",
   *logmax,lognext,*sum);*/

  return *logmax - log(rep+1) + log(*sum);
}
/*----------------------------------*/
void RandomOrder(int list[],int length)
     /*Returns a list of integers from 0 to Length-1 in random order. This is
      like having a bag of numbers.  Each time we pull one out and write it down.*/

{
  int i,pick,dummy;

  for (i=0;i<length;i++)
    list[i] = i;

  for (i=0;i<length;i++)
    {
      pick = RandomInteger(i,length-1);
      dummy = list[i];
      list[i] = list[pick];
      list[pick] = dummy;
    }

}
/*----------------------------------*/
double Factorial(int n)
     /*return n!, as a double*/
{
  double product = 1;
  int i;
  if (n < 0) printf("WARNING: trying to compute %d! \n", n);
  
  for (i=2; i<=n; i++)
    product *= i;

  return product;
}
/*------------------------------------------------------------*/
double ChiSq(int *list1,int len1,int *list2,int len2,int mincount,int missing,int *df)
/*This function takes two lists of integers (allele scores), and
computes a chi^sq for the equality of frequencies.  list1 and list2
are the two lists, of length len1, and len2, respectively.  Missing
data are marked by the value "missing".  Classes with fewer than
"mincount" observations are pooled.  The function returns the chisq
value, and a pointer to the resulting degrees of freedom (df).*/
{
  int *alleles;
  int *count1;
  int *count2;
  int *totcount;
  int numalleles=0;
  int i,allele;
  int data;
  int min1,min2,min3,all1,all2;
  double chisq = 0.0;
  int x1,x2;
  /*  double xbar; */
  int t1, t2; /*total allele counts*/
  double e1, e2; /*exp allele counts*/

  alleles = calloc(len1+len2,sizeof(int));
  count1 = calloc(len1+len2,sizeof(int));
  count2 = calloc(len1+len2,sizeof(int));
  totcount= calloc(len1+len2,sizeof(int));
  if ((alleles==NULL)||(count1==NULL)||(count2==NULL)||(totcount==NULL))
    {printf("Warning error assigning memory in ChiSq\n"); return 0;}
  
  for (i=0; i<len1+len2; i++)     /*compute allele counts*/
    {
      if (i<len1) data = list1[i];
      else data = list2[i-len1];
      if (data != missing)
	{
	  /*printf("%d ",data);*/
	  for (allele=0; allele<numalleles; allele++)
	    {if (data==alleles[allele]) break;}
	  if (allele==numalleles)
	    {
	      alleles[allele] = data;
	      numalleles++;	
	      totcount[allele] = 0;
	      count1[allele] = 0;
	      count2[allele] = 0;
	    }
	  totcount[allele]++;
	  if (i<len1) count1[allele]++;
	  else count2[allele]++;
	}
    }

  /*printf("\n\n\nNEXTLOCUS\n");
  printf("%d: ",numalleles);   
  for (allele=0; allele<numalleles; allele++)
    printf("%d: (%d,%d)  ",totcount[allele],count1[allele],count2[allele]);
    printf("\n"); */

  do        /*combine undersized classes*/
    { 
      min1 = min2 = min3 = len1+len2;     /*figure out 3 lowest counts*/
      all1 = all2 = 0;
      for (allele=0; allele<numalleles; allele++)  
	{
	  if (totcount[allele] < min1) 
	    {
	      min3 = min2;
	      min2 = min1; all2 = all1;
	      min1 = totcount[allele]; all1 = allele;
	    }
	  else if (totcount[allele] < min2) 
	    {
	      min3 = min2;
	      min2 = totcount[allele]; all2 = allele;
	    }
	  else if (totcount[allele] < min3) 
	    min3 = totcount[allele];
	}
      if (min1<=mincount)   /*combine classes 1 and 2*/
	{
	  if (all1 > all2) {allele = all2; all2 = all1; all1 = allele;}
	  totcount[all1] += totcount[all2];
	  count1[all1] += count1[all2];
	  count2[all1] += count2[all2];
	  totcount[all2] = totcount[numalleles-1];
	  count1[all2] = count1[numalleles-1];
	  count2[all2] = count2[numalleles-1];
	  numalleles--;
	}
      /*j++;
	if (j>10) printf("%d %d %d %d %d %d\n",min1,min2,min3,all1,all2,numalleles);*/
      /* melissa added parentheses 7/12/07 to fix ambiguous and/or statement */
    } while ((min1+min2 <= mincount || min3 <= mincount) && (numalleles>=1));

  /*printf("\nIn chisq: ");
  for (allele=0; allele<numalleles; allele++)
    printf("%d, ",totcount[allele]);
    printf("\n"); 
 
  printf("%d: ",numalleles);
  for (allele=0; allele<numalleles; allele++)
    printf("%d: (%d,%d)  ",totcount[allele],count1[allele],count2[allele]);
    printf("\n"); */

  *df = numalleles-1;  /*compute chisq*/
  t1 = t2 = 0;
  for (allele=0; allele<numalleles; allele++)
    {
      t1+= count1[allele]; 
      t2+= count2[allele]; 
    }
  /*printf("Phenotype totals = %d %d\n",t1,t2);*/

  for (allele=0; allele<numalleles; allele++)
    {
      x1 = count1[allele]; x2 = count2[allele];
      e1 = ((double) (x1+x2)*t1/(t1+t2));
      e2 = ((double) (x1+x2)*t2/(t1+t2));
      /*printf("O-E=%d--%1.3f   O-E=%d--%1.3f\n",x1,e1,x2,e2);  */
      
      chisq += ((double) (x1 - e1)*(x1 - e1)/e1);
      chisq += ((double) (x2 - e2)*(x2 - e2)/e2);
      
    }

  /*printf("Final chisq= %1.3f [mymath.c]\n",chisq);*/
  return chisq;
}
/*----------------------------------------------------*/

double mylgamma(double z)
{
/* LGAMMA function

   double_value = lgamma(<double_value > 0.>)

   returns the natural log of the gamma function

Uses Lanczos-type approximation to ln(gamma) for z > 0.   
Reference:                                          
 Lanczos, C. 'A precision approximation of the gamma   
    function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.   

Original was in FORTRAN
Accuracy: About 14 significant digits except for small regions   
          in the vicinity of 1 and 2.  
Programmer: Alan Miller   
          CSIRO Division of Mathematics & Statistics   
Latest revision - 17 April 1988   
  
Translated and modified into C by Peter Beerli 1997 
Tested against Mathematica's Log[Gamma[x]]
*/

  double a[9] = { 0.9999999999995183, 676.5203681218835,
    -1259.139216722289, 771.3234287757674, -176.6150291498386,
    12.50734324009056, -0.1385710331296526, 9.934937113930748e-6,
    1.659470187408462e-7 };
  double lnsqrt2pi = 0.9189385332046727;
  double result;
  long j;
  double tmp;
  if (z <= 0.)
    {
      fprintf(stderr,"lgamma function failed with wrong input (%f)\n",z);
      assert(0);
      exit(-1);
    }
  result = 0.;
  tmp = z + 7.;
  for (j = 9; j >= 2; --j)
    {
      result += a[j - 1] / tmp;
      tmp -= 1.;
    }
  result += a[0];
  result = log (result) + lnsqrt2pi - (z + 6.5) + (z - 0.5) * log (z + 6.5);
  return result;
}  /* lgamma */
