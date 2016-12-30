
/*

   Part of structure.c.  

   This bit is in charge of reading in the information from the datafile, and 
   preparing it. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structure.h"
#include "params.h"

FILE *INPUT;			/*Data file */

void OpenData ();
void ReadData (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase,
                int *Recessive);
void WarnEOF (int ind);
int CheckIfValidInt (char intstr[ALLELLEN], int ind, int loc);
void Mismatch (int ind, char description[30]);

void ExtraData ();
void PrintSomeData (int *Geno, struct IND *Individual, FILE * file);
void CountLineLens ();
void RecessiveDataCheck(int *Geno, int *Recessive);


/*================================================*/
void 
ReadInputFile (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase,
                int *Recessive)
{
  OpenData ();			/*open input file */
  printf ("Reading file \"%s\".\n", DATAFILE);
  ReadData (Geno, Mapdistance, Markername, Individual, Phase,Recessive);
  fclose (INPUT);
  if (ECHODATA)
    PrintSomeData (Geno, Individual, stdout);
  if (RECESSIVEALLELES) RecessiveDataCheck(Geno,Recessive);

}
/*-------------------------------------------*/
void 
OpenData ()
{
  /*  int trouble = 0; */
  INPUT = fopen (DATAFILE, "r");

  if (INPUT == NULL)
    {
      printf ("Unable to open the file %s.\n", DATAFILE);
      Kill ();
    }

}


/*-------------------------------------------*/
void 
ReadData (int *Geno, double *Mapdistance, char *Markername, struct IND *Individual, double *Phase, int *Recessive)
{
  int ind;
  int loc;
  int line;
  int col;
  char label[LABELLEN];
  char intstr[ALLELLEN];
  char gene[GENELEN];
  char next;
  int pop=-1;
  int popflag=-1;
  int phen=-1;
  int strlength;
  int valid = 1;
  /*  int polymorphismcounter; */
  /* int toggle; */
  int rowsperind;  /*the next 2 variables allow single-row input for data*/
  int colsperloc;  /*number of columns per locus = 1 or PLOIDY depending on format*/ 
  int i,j;
  /*  int loc1,loc2,temp; */
  int *inputorder = malloc(NUMLOCI*sizeof(int));
  int location=-1;
  int phenotypecol=-1;
  /*  int numcols; */
  int CheckIfValidDouble (char intstr[ALLELLEN], int ind,int loc);


  ind = -1;

 /*inputorder is initialised with the loci in the right order */
  for (loc=0;loc<NUMLOCI;loc++)
  inputorder[loc]=loc;

/* The following part (which should normally be commented out) randomizes the order of the loci*/
/* while leaving mapdistances, etc, intact. */
/* this can be useful when testing if there is significant map distance info in the sample. */
/*
  for (ind=0;ind<NUMLOCI*100;ind++)
{
loc1=RandomInteger(0,NUMLOCI-1);
loc2=RandomInteger(0,NUMLOCI-1);
temp=inputorder[loc1];
inputorder[loc1]=inputorder[loc2];
inputorder[loc2]=temp;
}
*/
  /*Read in locus information*/
  
    if (MARKERNAMES) {
      for (loc = 0; loc < NUMLOCI; loc++)
      {
	strlength = ReadString (gene, GENELEN, INPUT);
	if (strlength == 0) {
	  printf("%i MARKERNAMES\n", loc);WarnEOF (ind);}
	for (j=0; j<strlength; j++)
	  Markername[MarkernamePos(loc,j)] = gene[j];
      }
    }
  if (RECESSIVEALLELES)
    for (loc = 0; loc < NUMLOCI; loc++)
      {
	strlength = ReadString (intstr, ALLELLEN, INPUT);
	if (strlength == 0)
	  WarnEOF (ind);
	Recessive[loc] = (double) atoi (intstr);
      }
  if (MAPDISTANCES)
    for (loc = 0; loc < NUMLOCI; loc++)
      {
	strlength = ReadString (intstr, ALLELLEN, INPUT);
	if (strlength == 0)
	  WarnEOF (ind);
	
	/* daniel made the change: change from atoi to atof on Dec 24, 2002 */
	Mapdistance[loc] = (double) atof (intstr);
      }
  /*End reading in locus information-------------------------------
    Read in individual/genotype data-------------------------------*/
  
  if (ONEROWPERIND) {rowsperind = 1; colsperloc = LINES;}
  else {rowsperind = LINES; colsperloc = 1;}
  
  for (ind = 0; ind < NUMINDS; ind++)
    {
      for (line = 0; line < rowsperind; line++)
	{
	  if (LABEL)		/*read in extra data for each individual */
	    {
	      strlength = ReadString (label, LABELLEN, INPUT);
	      if (strlength == 0)
		WarnEOF (ind);
	    }
	  if (POPDATA)
	    {
	      strlength = ReadString (intstr, ALLELLEN, INPUT);
	      if (strlength == 0)
		WarnEOF (ind);
	      valid = valid * CheckIfValidInt (intstr, ind, loc);
	      pop = atoi (intstr);
	    }
	  if (POPFLAG)
	    {
	      strlength = ReadString (intstr, ALLELLEN, INPUT);
	      if (strlength == 0)
		WarnEOF (ind);
	      valid = valid * CheckIfValidInt (intstr, ind, loc);
	      popflag = atoi (intstr);
	    }
	  /* melissa added 7/12/07 */
	  if (LOCDATA==1) {
	    strlength = ReadString(intstr, ALLELLEN, INPUT);
	    if (strlength ==0)
	      WarnEOF(ind);
	    valid = valid + CheckIfValidInt(intstr, ind, loc);
	    location = atoi(intstr);
	  }
	  if (PHENOTYPE) /*-------------------------*/
	    {
	      if (PHENOTYPECOL == -9) phenotypecol = 0;   /*if STRATPARAMS is not read, then
							    this is given a value of -9*/

	      else if (PHENOTYPECOL <= LABEL + POPDATA + POPFLAG)  /*too small*/
		{
		  printf("Error: PHENOTYPECOL (now set in STRATPARAMS as %d) must be at least %d\n",
			PHENOTYPECOL, LABEL + POPDATA + POPFLAG + 1);
		  printf("given your current values of LABEL, POPDATA and POPFLAG\n");
		  Kill();
		}
	      else if (PHENOTYPECOL > LABEL + POPDATA + POPFLAG + 1 + EXTRACOLS)  /*too large*/
		{
		  printf("Error: PHENOTYPECOL (set in STRATPARAMS as %d) must be at most %d,\n",
			 PHENOTYPECOL, LABEL + POPDATA + POPFLAG + 1 + EXTRACOLS);
		  printf("  given your current values of LABEL, POPDATA, POPFLAG and EXTRACOLS\n");
		  Kill();
		}
	      else phenotypecol = PHENOTYPECOL - LABEL - POPDATA - POPFLAG -1;   /*figure out which
									           column to read*/

	      /*printf("phenotypecol = %d; EXTRACOLS = %d\n",phenotypecol, EXTRACOLS);*/
	      for (col = 0; col < EXTRACOLS + 1; col++) 
		{
		  if (col == phenotypecol)
		    {
		      strlength = ReadString (intstr, ALLELLEN, INPUT);
		      if (strlength == 0)
			WarnEOF (ind);
		      valid = valid * CheckIfValidInt (intstr, ind, loc);
		      phen = atoi (intstr);
		      /*printf("phenotypecol = %d, col = %d, phen = %d\n",phenotypecol,col,phen);*/
		    }
		  else ReadString (intstr, ALLELLEN, INPUT);   /*skip these data*/
		}
	    
	    }          /*---End if (Phenotype)----------------------*/

	  else for (col = 0; col < EXTRACOLS; col++)  /*no phenotypes*/
	    ReadString (intstr, ALLELLEN, INPUT);


	  
	  if (line == 0)		/*save the preliminary data */
	    {
	      if (LABEL)
		strcpy (Individual[ind].Label, label);
	      if (POPDATA)
		Individual[ind].Population = pop;
	      if (POPFLAG)
		Individual[ind].PopFlag = popflag;
	      if (PHENOTYPE)
		Individual[ind].Phenotype = phen;

	      /* melissa added 7/12/07 */
	      if (LOCDATA)
		Individual[ind].Location = location;
	      if (LOCISPOP) 
		Individual[ind].Location = pop;
	    }
	  if (line > 1)		/*check for consistency across lines */
	    {
	      if (LABEL)
		if ((strcmp (Individual[ind].Label, label)))
		  {
		    Mismatch (ind, "label");
		    valid = 0;
		  }
	      /*printf("The labels are %s and %s\n",Individual[ind].Label,label);} */
	      if (POPDATA)
		if (Individual[ind].Population != pop)
		  {
		    Mismatch (ind, "pop");
		    valid = 0;
		  }
	      if (POPFLAG)
		if (Individual[ind].PopFlag != popflag)
		  {
		    Mismatch (ind, "popflag");
		    valid = 0;
		  }
	      if (PHENOTYPE)
		if (Individual[ind].Phenotype != phen)
		  {
		    Mismatch (ind, "phenotype");
		    valid = 0;
		  }
	    }
	  
	  for (loc = 0; loc < NUMLOCI; loc++)    /*read in genotype data here*/
	    for (i=0; i<colsperloc; i++)
	      {
		strlength = ReadString (intstr, ALLELLEN, INPUT);
		if (strlength == 0) {printf("readlociEOF\n");
		WarnEOF (ind);}
		valid = valid * CheckIfValidInt (intstr, ind, loc);
		if (ONEROWPERIND) Geno[GenPos (ind, i, inputorder[loc])] = atoi (intstr);
		else Geno[GenPos (ind, line, inputorder[loc])] = atoi (intstr);
	      }
	/* printf(" % 4ld % 4ld % 4ld     .....    % 4ld % 4ld \n",Geno[GenPos(ind,line,0)],Geno[GenPos(ind,line,1)],Geno[GenPos(ind,line,2)],Geno[GenPos(ind,line,NUMLOCI-2)],Geno[GenPos(ind,line,NUMLOCI-1)]); */
	}

      if (PHASEINFO) /*possibly read in row of phase information*/
	for (loc=0;loc<NUMLOCI;loc++)
	  {
	    strlength = ReadString (intstr, ALLELLEN, INPUT);
	    if (strlength == 0) WarnEOF (ind);
	    valid=valid * CheckIfValidDouble(intstr,ind,loc);
	    Phase[PhasePos(ind,loc)]=atof(intstr);
	    /*check that probabilities are in [0,1]*/
	    if (Phase[PhasePos(ind,loc)]>1.0 || Phase[PhasePos(ind,loc)]<0.0)
	      {
		printf("Phase information for individual %d locus %d (%1.3f) is not a real in [0.0,1.0] \n",ind,loc,Phase[PhasePos(ind,loc)]);
		valid=0;
	      }
	  }
    } /*end of loop over individuals*/


  /*check if anything else left in the input file*/
  do
    {
      next = fgetc (INPUT);
      if ((!Whitespace (next)) && (next != EOF))
	{
	  ExtraData ();
	  valid = 0;
	  break;
	}
    }
  while (next != EOF);

  if (!(valid))
    {
      CountLineLens ();
      Kill ();
    }
  free(inputorder);
}
/*------------------------------------*/
void RecessiveDataCheck(int *Geno, int *Recessive)
    /* this function checks whether any genotypes have both a recessive and dominant
     * allele.  If any do, it terminates the program.
     * in the polyploid case, it also checks whether the NOTAMBIGUOUS code is anywhere 
     * in the datset. If it is it also terminates the program. */
{
  int ind, line, loc;
  int rec, dom;
  int error=0;

  int recessive_shown = 0;
  int recessive_allele_shown = 0;
  int recessive_allele_not_shown = 0;
  int unambiguous_loci = 0;
  int ambiguous_norecessive = 0;
  for (loc=0; loc<NUMLOCI; loc++)
    {
      
      if (Recessive[loc] != MISSING)
	{
	  recessive_shown=0;
	  for (ind=0; ind<NUMINDS; ind++)
	    {
	      rec=0; 
	      dom=0;
	      for (line=0; line<LINES; line++)
		{
		  if (Geno[GenPos (ind,line,loc)] == Recessive[loc]) {
		    rec=1;
		    recessive_shown=1;
		  }
		  else if (Geno[GenPos (ind,line,loc)] != MISSING) dom=1;
		}
	      if (rec*dom==1 && error<100) 
		{
		  printf("WARNING: both recessive and dominant alleles in genotype: ind=%d, loc=%d\n",ind+1,loc+1);
		  error++;
		}

	    }
	  
	  if(recessive_shown) 
	    recessive_allele_shown++;
	  else
	    if(Recessive[loc]!= NOTAMBIGUOUS)
	      recessive_allele_not_shown++;
	}
    }
  if (error>100) printf("Total of %d such errors\n",error);
  if (error>0) {printf("Terminating program.\n"); Kill();}
  if (LINES>2)
    for (loc=0;loc<NUMLOCI;loc++){
      
      if(Recessive[loc] == NOTAMBIGUOUS)
	unambiguous_loci++;
      if(Recessive[loc] == MISSING)
	ambiguous_norecessive++;
      for (ind=0;ind<NUMINDS;ind++)
	for (line=0;line<LINES;line++)
	  if (Geno[GenPos(ind,line,loc)]==NOTAMBIGUOUS)
	    {
	      printf("WARNING: the code for NOTAMBIGUOUS alleles, %d, appears  at least once in the dataset at : ind=%d, loc=%d\n", NOTAMBIGUOUS, ind+1,loc+1);
	      Kill(); /* modify by William - kill (small case) is not a standard ANSI C function */
	    }
  
    }
  printf("\n");
  if(LINES>2){
    printf("Number of loci without genotypic ambiguity: %d \n",unambiguous_loci);
    printf("Number of loci with ambiguous copy numbers but no recessive alleles: %d \n",ambiguous_norecessive);
  }
  printf("Number of loci with recessive allele present in data: %d\n", recessive_allele_shown);
  printf("Number of loci with recessive allele absent in data: %d\n", recessive_allele_not_shown);
  fflush(stdout);

  
}
/*------------------------------------*/
void 
WarnEOF (int ind)
{
  /*This function is called from ReadData if there is an unexpected end of file */

  printf ("\n\nWARNING:  Unexpected end of input file.  The details of the\n");
  printf ("input file are set in mainparams.  I ran out of data while reading\n");
  printf ("the data for individual %d.\n\n", ind + 1);

  CountLineLens ();
  Kill ();
}
/*------------------------------------*/
int 
CheckIfValidInt (char intstr[ALLELLEN], int ind, int loc)
{
  /*This function checks for non-numeric data in the input file (not
     used when reading the labels).  Returns 1 if valid, otherwise 0. */

  int i;
  int ok = 1;

  for (i = 0; i < ALLELLEN; i++)
    {
      if (intstr[i] == '\0')
	break;
      if (((intstr[i] < '0') || (intstr[i] > '9')) && (intstr[i] != '-'))
	ok = 0;
    }

  if (ok == 0)
    {
      printf ("\nWARNING! Probable error in the input file.  \n");
      printf ("Individual %d, locus %d;  ", ind + 1, loc + 1);
      printf ("encountered the following data \n");
      printf ("\"%s\" when expecting an integer\n\n",intstr);
    }

  return ok;
}
/*------------------------------------*/
int 
CheckIfValidDouble (char intstr[ALLELLEN], int ind,int loc)
{
  /*This function checks for non-numeric data in the input file (not
     used when reading the labels).  Returns 1 if valid, otherwise 0. */

  int i;
  int ok = 1;
  int count=0;

  for (i = 0; i < ALLELLEN; i++)
    {
      if (intstr[i] == '\0')
	break;
      if (((intstr[i] < '0') || (intstr[i] > '9')) && (intstr[i] != '-') && (intstr[i] != '.'))
	ok = 0;
      if (intstr[i]=='.')
      count=count+1;
    }
  if (count>1)
    ok=0;
  if (ok == 0)
    {
      printf ("\nWARNING! Possible error in the input file.  Non-Real \n");
      printf ("number in an unexpected place: encountered %s while reading\n", intstr);
      printf ("the data for individual %d\n", ind + 1);
      printf ("and locus %d\n", loc + 1);
    }

  return ok;
}
/*------------------------------------*/
void 
Mismatch (int ind, char description[30])
{
  printf ("\nWARNING! Possible error in the input file. The value of %s\n", description);
  printf ("does not agree across lines for individual %d\n", ind + 1);
}
/*------------------------------------*/
void 
ExtraData ()
{
  printf ("\n\nWARNING:  There may be more data in the input file\n");
  printf ("than indicated by the program constants.  Check the values\n");
  printf ("entered for NUMLOCI and NUMINDS, etc, in the program constants.\n\n\n");
}
/*------------------------------------*/
void 
PrintSomeData (int *Geno, struct IND *Individual, FILE * file)
{
  int ind, line;
  int col = 0;
  int precols;
  /*  int loc; */

  fprintf (file, "\n\nData file \"%s\" (truncated) --\n\n", DATAFILE);
  fprintf (file, "Ind:   ");
  if (LABEL)
    {
      fprintf (file, "%5s ", "Label");
      col++;
    }
  if (POPDATA)
    {
      fprintf (file, "Pop ");
      col++;
    }
  if (POPFLAG)
    {
      fprintf (file, "Flag ");
      col++;
    }
  if (PHENOTYPE)
    {
      fprintf (file, "Phen ");
      col++;
    }
  precols = col;
  if (precols > 1)
    fprintf (file, ": ");
  fprintf (file, "Genotype_data . . . .\n");

  for (ind = 0; ind < NUMINDS; ind++)
    for (line = 0; line < LINES; line++)
      {
	fprintf (file, "%3d: ", ind + 1);
	if (LABEL)
	  fprintf (file, "%5s ", Individual[ind].Label);
	if (POPDATA)
	  fprintf (file, "%4d ", Individual[ind].Population);
	if (POPFLAG)
	  fprintf (file, "%4d ", Individual[ind].PopFlag);
	if (PHENOTYPE)
	  fprintf (file, "%4d ", Individual[ind].Phenotype);
	if (precols > 1)
	  fprintf (file, " :");
	col = precols;
	while ((col < 9) && (col < precols + NUMLOCI))
	  {
	    fprintf (file, "%3d ", Geno[GenPos (ind, line, col - precols)]);
	    col++;
	  }
	if (NUMLOCI > col - precols)
	  fprintf (file, " . . . . %3d\n", Geno[GenPos (ind, line, NUMLOCI - 1)]);
	else
	  fprintf (file, "\n");

	if (NUMINDS > LINES * MAXINDS)
	  if ((ind == MAXINDS) && (line == LINES-1) && (NUMINDS-MAXINDS > ind))
	    {
	      ind = NUMINDS - MAXINDS;
	      fprintf (file, "\n      *******   \n\n");
	    }
      }
  fprintf (file, "\n");
}
/*------------------------------------*/
void CountLineLens ()
/*This function is called from function ReadData if there
   are problems with the input file.  It counts the number of
   words in each row of the input file, and then prints this
   to the screen */
{
  int label, popdata, popflag, phenotype;
  char last, next;
  int row, word;
  int numrows;
  int *numwords;
  int datalines;
  int totalwords;
  int printed;
  int consec;
  int min, max;
  if (LABEL)
    label = 1;
  else
    label = 0;
  if (POPDATA)
    popdata = 1;
  else
    popdata = 0;
  if (POPFLAG)
    popflag = 1;
  else
    popflag = 0;
  if (PHENOTYPE)
    phenotype = 1;
  else
    phenotype = 0;

  printf ("----------------------------------\n");
  printf ("There were errors in the input file (listed above). According to \n");
  printf ("\"mainparams\" the input file should contain ");
  if (MARKERNAMES)
    printf("one row of markernames with %d entries,\n",NUMLOCI);
  if (RECESSIVEALLELES)
    printf ("one row indicating recessive alleles with %d entries,\n", NUMLOCI);
  if(MAPDISTANCES)
    printf ("one row of map distances with %d entries,\n", NUMLOCI);

  if (ONEROWPERIND)
    printf(" %d rows with %d entries. ",NUMINDS, label + popdata + popflag + phenotype + EXTRACOLS + LINES*NUMLOCI);
  else printf(" %d rows with %d entries ",LINES * NUMINDS, label + popdata + popflag + phenotype + EXTRACOLS + NUMLOCI);
  if (PHASEINFO)
    printf("\nplus an additional row for each individual (after the first genotypes) \ncontaining %d entries of phase information",NUMLOCI);
printf (".\n\n");  

  fclose (INPUT);
  OpenData ();			/*start reading input file from the beginning */
  numrows = 0;			/*count number of rows */
  do
    {
      next = fgetc (INPUT);
      if (next == '\n')
	numrows++;
    }
  while (next != EOF);
  numrows++;

  numwords = calloc (numrows, sizeof (int));
  if (numwords == NULL)
    printf ("Warning: error assigning memory in CountLineLens\n");
  else
    {
      fclose (INPUT);
      OpenData ();		/*start reading input file from the beginning */
      last = ' ';
      word = 0;
      row = 0;
      do
	{
	  next = fgetc (INPUT);
	  if (Whitespace (last) && (!(Whitespace (next))) && (next != EOF))
	    word++;
	  if (next == '\n')
	    {
	      numwords[row] = word;
	      row++;
	      word = 0;
	    }
	  last = next;
	}
      while (next != EOF);
      numwords[row] = word;
      for (row = numrows - 1; row >= 0; row--)	/*get rid of extra zeroes at end of file */
	if (numwords[row] > 0)
	  break;
      numrows = row + 1;

      datalines = 0;
      totalwords = 0;
      min = max = numwords[0];
      for (row = 0; row < numrows; row++)
	{
	  if (numwords[row] > 0)
	    datalines++;
	  totalwords += numwords[row];
	  if (((min == 0) || (min > numwords[row])) && (numwords[row] > 0))
	    min = numwords[row];
	  if (max < numwords[row])
	    max = numwords[row];
	}

      printf ("There are %d rows of data in the input file, with an average of %1.2f\n",
	      datalines, (double) totalwords / datalines);
      printf ("entries per line.  The following shows the number of entries in each\n");
      printf ("line of the input file:\n\n");
      printf ("# Entries:   Line numbers\n");

      /*for (row=0; row<numrows; row++)
         printf("%d %d\n",row,numwords[row]);
         printf("%d\n",numrows); */

      for (word = 0; word <= max; word++)
	{
	  printed = 0;
	  consec = 0;
	  for (row = 0; row < numrows; row++)
	    {
	      if (numwords[row] == word)
		{
		  if (!(printed))	/*print number of words on line */
		    if (!((row == numrows - 1) && (word == 0)))
		      printf ("     %4d:   ", word);

		  if (row == 0)	/*print value of first row */
		    printf ("%d", row + 1);
		  else if (numwords[row - 1] != word)	/*new value */
		    {
		      if (!((row == numrows - 1) && (word == 0)))	/*except last line empty */
			{
			  if (printed)
			    printf (", ");	/*not first occurrence */
			  printf ("%d", row + 1);
			}
		    }
		  else if (row == numrows - 1)
		    {		/*print if last row in file */
		      if (consec == 1)
			printf (", %d", row + 1);
		      if (consec > 1)
			printf ("--%d", row + 1);
		      consec = 0;
		    }

		  consec++;
		  printed = 1;

		}
	      else if ((row > 0) && (numwords[row - 1] == word))
		{
		  if (consec == 2)
		    printf (", %d", row);
		  if (consec > 2)
		    printf ("--%d", row);
		  consec = 0;
		}

	    }
	  if (printed)
	    printf ("\n");
	  if (word == 0)
	    word = min - 1;
	}
      free (numwords);
    }
  printf ("----------------------------------\n");
}
/*---------------------------------------*/
int FindMaxAlleles (int *Geno, int *Recessive)
     /*Count the number of alleles at each locus, and return the largest
        number */
{
  int pos, k;
  int ind, line, loc;
  int value, recessivecoded;
  int mink = LINES * NUMINDS + 1;	/*start at max possible */
  int maxk = 0;			/*start at min possible */
  int sumk = 0;
  int *Alleles;

  Alleles = calloc (LINES * NUMINDS, sizeof (int));

  for (loc = 0; loc < NUMLOCI; loc++)
    {
      k = 0;
      recessivecoded = 0;
      for (ind = 0; ind < NUMINDS; ind++)
	for (line = 0; line < LINES; line++)
	  {
	    value = Geno[GenPos (ind, line, loc)];
	    if (value != MISSING)
	      {
		for (pos = 0; pos < k; pos++)
		  if (Alleles[pos] == value)
		    break;

		if (pos == k)
		  {
		    Alleles[k] = value;
		    k++;
		    if (RECESSIVEALLELES && value == Recessive[loc])
		      recessivecoded = 1;
		  }
	      }
	  }
      /*add one if an unobserved recessive needs to be encoded */
      if (RECESSIVEALLELES && Recessive[loc] != MISSING
	  && recessivecoded == 0 && Recessive[loc]!=NOTAMBIGUOUS)
	k += 1;
      if (maxk < k)
	maxk = k;
      if (mink > k)
	mink = k;
      sumk += k;
    }

  printf ("Number of alleles per locus: min=%2d; ave=%2.1f; max=%2d\n",
	  mink, (double) sumk / NUMLOCI, maxk);
  free (Alleles);
  return maxk;
}

/*---------------------------------------*/
void CountAlleles (int *Geno, int *NumAlleles, int *Translation, int *Recessive)
     /*This function records the number of alleles at each
        locus, and recodes the alleles in Geno to be in {0..k-1}.
        Translation stores the coding information. */
{
  int pos, k;			/*pos= position in array; k= alleles so far */
  int ind, line, loc;
  int value;			/*value of current allele */
  int maxk = 0;
  int newmissing,newnotambiguous, recessivecoded;
  /*worry about whether missing data value interferes with recoding */
  ORIGMISSING = MISSING;	/*store value of MISSING in original data */
  if ((MISSING >= 0) && (MISSING < MAXALLELES+1))
    newmissing = -9;
  else
    newmissing = MISSING;
    /*set new value of notambiguous 1 less than missing */
    newnotambiguous=newmissing-1;

  /*recode all data */
  for (loc = 0; loc < NUMLOCI; loc++)
    {
      if (RECESSIVEALLELES && Recessive[loc] == MISSING)
	Recessive[loc] = newmissing;
      if (RECESSIVEALLELES && Recessive[loc]==NOTAMBIGUOUS) 
{
  Recessive[loc]=newnotambiguous;
  /* printf("locus %d is NOTAMBIGUOUS\n" ,loc); */
}
      recessivecoded = 0;
      k = 0;
      for (ind = 0; ind < NUMINDS; ind++)
	for (line = 0; line < LINES; line++)
	  {
	    value = Geno[GenPos (ind, line, loc)];
	    if (value == MISSING)
	      Geno[GenPos (ind, line, loc)] = newmissing;
	    else
	      {
		for (pos = 0; pos < k; pos++)
		  if (Translation[TransPos (loc, pos)] == value)
		    break;
		if (pos == k)
		  {
		    Translation[TransPos (loc, pos)] = value;
		    k++;
		    if (RECESSIVEALLELES && value == Recessive[loc])
		      {
			Recessive[loc] = pos;
			recessivecoded = 1;
		      }
		  }
		Geno[GenPos (ind, line, loc)] = pos;	/*recoding */
	      }
	  }
/*add an extra allele null allele if there are null alleles at the locus but 
the homozygous null is present nowhere in the sample*/
      if (RECESSIVEALLELES && Recessive[loc] != newmissing
	  && recessivecoded == 0  && k>0 && Recessive[loc]!=newnotambiguous)
	{
	  Recessive[loc] = k;
	  Translation[TransPos (loc, k)] = 29697;
	  k++;
/*this value chosen as unlikely to be a number that anyone is going to choose */
	}
      NumAlleles[loc] = k;
      if (k==0) 
	printf("WARNING: locus %d has all missing data.  This might cause unexpected behaviour.\n",loc);

      if (maxk < k)
	maxk = k;
    }

  if (maxk > MAXALLELES)
    {
      printf ("maxk = %d, MAXALLELES = %d\n", maxk, MAXALLELES);
      printf ("Program error in function CountAlleles\n");
      Kill ();
    }
  MISSING = newmissing;
  NOTAMBIGUOUS=newnotambiguous;
}
