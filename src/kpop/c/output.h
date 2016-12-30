/* #include "structure.h" */
void
DataCollection (int *Geno,int *PreGeno, double *Q, double *QSum, int *Z, int *Z1, double *SiteBySiteSum, double *P, double *PSum,
		double *Fst, double *FstSum, int *NumAlleles,
		int *AncestDist, double *Alpha, double *sumAlpha,
		double *sumR, double *varR, double *like,
		double *sumlikes, double *sumsqlikes, double *R,
		double *Epsilon, double *SumEpsilon, double recomblikelihood,
		double *lambda, double *sumlambda, int *Recessive,
		double *PopPrior, double *PopPriorSum, int PopPriorLen,
		double *sumindlikes, double *indlikes_norm, int rep);

void PrintLike (double like, int rep, int *Geno, int *PreGeno, double *Q, double *P, 
		double recomblikelihood);
void PrintUpdate (int rep, int *Geno, int *PreGeno, double *Alpha, double *Correls, 
		  double *P, double *Q, double like, double sumlikes, 
		  double sumsqlikes, int *NumAlleles, double *R, double *lambda, 
		  struct IND *Individual,  double recomblikelihood, int *Recessive, double *PopPrior, int PopPriorLen);
void
OutPutResults (int *Geno, int rep, int savefreq, struct IND *Individual,
	  double *PSum, double *QSum, double *SiteBySiteSum, double *FstSum,
	       int *AncestDist, double *UsePopProbs,
	       double sumlikes, double sumsqlikes, double *sumAlpha,
	       double *sumR, double *varR,
	    int *NumAlleles, int *Translation, int final,
	       char *Markername, double *R, double *SumEpsilon,
	       double *lambda, double *sumlambda, double *PopPriorSum,
	       int PopPriorLen, double *sumindlikes, double *indlikes_norm, 
	       int argc, char *argv[]);
double CalcLikeInd(int *Geno, int *PreGeno, double *Q, double *P,
		   int ind, int *Recessive);
