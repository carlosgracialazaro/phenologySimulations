// Description: Dynamics of a mutualist evologyc system

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "eigenvalues.h"

////////////////////////////////
// PARAMETERS:
////////////////////////////////

#define heterogeneousBeta 0 // heterogeneous (or not) interspecific competitive term
#define maximumN_A 57 // maximum number of animal species
#define maximumN_P 27 // maximum number of plant species
#define transient 2000 // transient Time (transient/delta_t = number of time steps in differential equations integration)
//#define intervalStepsForEvolutionFiles 100 // outpout of evolution into files every intervalStepsForEvolutionFiles steps
//#define rewiringsPerSwitchOfNetwork 10
#define nullModelNetwork 0 // 1 for randomization of the network according to a probabilistic null model algorithm
#define rewiringsTowardLowerNestedness 0 // if 1: {we accept a rewiring if the nestedness is lower than the previous one, otherwise it will reject it}

#define networks 100 // number of different randomizations of phenology coefficients


float delta_t=0.01; // time increment in differential equations integration

// initialized through parameters to main ():
char asymmetric; // additional competitive term for plants, not driven by pollinators (light, soil, etc...). // asymmetric = 1 -> both mechanisms in plants // asymmetric = 2, only MF in plants
char phenology; // parameter #15 to main, if (phenology=1) we consider time overlaps (phenology coefficients)
//float min_betaP_0; // plants competition constant (cross-species), minimum value to sweep the interval
//float max_betaP_0; //  plants competition constant (cross-species), maximum value to sweep the interval
//float delta_betaP_0; // plants competition constant (cross-species), increment to sweep the interval
float betaP_ii_0; // plants competition constant (intra-specie)
//float min_gammaP_0; // plants mutualism constant, minimum value to sweep the interval
//float max_gammaP_0; // plants mutualism constant, maximum value to sweep the interval
//float delta_gammaP_0; // plants mutualism constant, increment to sweep the interval
float min_alphaP; // plants growing constant, minimum value not to sweep the interval but to assign at random
float max_alphaP; // plants growing constant, maximum value not to sweep the interval but to assign at random
float hP; // holling term for plants (constant) ; ~(N_A*gammaP_0)?
int rewiringsBeforeFirstNetwork; //number of rewirings before first realization
char competitionMeanField; // 0 for our formula, 1 for meanfield considering projection matrices, 2 for Bastolla et al. formula (asymmetric must be 0 for competitionMeanField 1 or 2)
float betaP_i[maximumN_P]; //plants competition constant (cross-species) per specie
float betaA_i[maximumN_P]; // plants mutualism constant (cross-species) per specie


// swept values:
float betaP_0; //plants competition constant (cross-species)
float gammaP_0; // plants mutualism constant

// inicialited by symmetry (plants constant = animal constant):
float betaA_0; // animals competition constant (cross-species), to sweep the interval
float betaA_ii_0; // animals competition constant (intra-specie)
float gammaA_0; // animals mutualism constant, to sweep the interval
float min_alphaA; // animals growing constant, minimum value not to sweep the interval but to assign at random
float max_alphaA; // animals growing constant, maximum value not to sweep the interval but to assign at random
float hA; // holling term for animals (constant) ; ~(N*gammaA_0)?

////////////////////////////////
// GLOBAL VARIABLES:
////////////////////////////////

int N_A; // number of animal species
int N_P; // number of plant species
char K [maximumN_P][maximumN_A]; // bipartite adjacency network
int kP [maximumN_P]; // connectivity of a plant in the bipartite graph, i.e., number of animal species that pollinate this plant
int kA [maximumN_A]; // connectivity of an animal in the bipartite graph, i.e., number of plant species that this animal pollinates
int WP [maximumN_P][maximumN_P]; // weighted projection matrix (plants)
int WA [maximumN_A][maximumN_A]; // weighted projection matrix (animals)
float SP [maximumN_P]; // abundances vector (plants)
float SA [maximumN_A]; // abundances vector (animals)
float alphaP [maximumN_P]; // growing factor (plants)
float alphaA [maximumN_A]; // growing factor (animals)

// Phenology:
float OmegaPijk [maximumN_P][maximumN_P][maximumN_A];
float PhiPik [maximumN_P][maximumN_A];
float OmegaAkli [maximumN_A][maximumN_A][maximumN_P];
float PhiAki [maximumN_A][maximumN_P];

////////////////////////////////
// DECLARATIONS:
////////////////////////////////

struct rewiringKey 
{
   char plant; // 1 for plants, 0 for animals
   int who;
   int withWho;
   int unplug;
   int plug; 
};

float randomFloat (); // random [0,1)
void readBipartiteMatrix (char* networkName, int nullNumber); // read bipartite matrix file
void initialConditions (); // to be called after reading bipartite matrix
void equationDynamics (); // differential equations integration
float plantStep (int i);
float animalStep (int i);
//void printAbundancesVectorEvolution (int step);
void printAbundancesVector (float nestedness, float nestednessNODF, int diversity, int networkNumber); // it writes into file
void printAverage (float bd);

void printProjectedMatrices (char* matrixFile, int re); // parameter: argv [1], # rewiring
void globalVariablesInitialization (int argc, char **argv);
int biodiversity (); // it returns the total number of species
void projectionMatrices (); //it computes weighted projection matrices Wp, WA
struct rewiringKey rewiring ();
void undoRewiring (struct rewiringKey r);
float largestEigenvalue (); // it returns the maximum eigenvalue of the adjacency matrix
void nullModel (); // randomize the network according to a probabilistic NULL MODEL
float NODF (); // it computes the nestesned according to NODF measure
float swap (); // it makes an effective rewiring


///////////////////////////////
// ALGORITHMS:
//////////////////////////////

float randomFloat () // random [0,1)
{
  float r;
  int i;
  i=(rand()% RAND_MAX);
  r=1.*i/RAND_MAX;
  return r;
}


void initialConditions () // to be called after reading bipartite matrix
{
	int j;

	for	(j=0; j<N_P; j++)
	{
		alphaP [j] = randomFloat () * (max_alphaP - min_alphaP) + min_alphaP;
		SP [j] = randomFloat () * 0.9 + 0.05; // abundances vector (plants), from 0.05 to 0.95
		if (heterogeneousBeta)
		    betaP_i[j]= betaP_0 * (1.8 * randomFloat() + 0.1);
		else
		    betaP_i[j]= betaP_0;
	}
	for	(j=0; j<N_A; j++)
	{
		alphaA [j] = randomFloat () * (max_alphaA - min_alphaA) + min_alphaA;
		SA [j] = randomFloat () * 0.9 + 0.05; // abundances vector (animals), from 0.05 to 0.95
		if (heterogeneousBeta)
		    betaA_i[j]= betaA_0 * (1.8 * randomFloat() + 0.1);
		else
		    betaA_i[j]= betaA_0;
	}
}


float plantStep (int i)
{
    int j, k, l;
    float denominator = 0.;
    float suma, delta;

    // growing term
	float increment = alphaP [i];

    // cross-species competition term:
    if (asymmetric == 0)
    {
     delta = 0.;
	 if (competitionMeanField>0)
	 {
	  for (j=0; j<N_P; j++)
	   if (i!=j)
	   {
	    if (competitionMeanField == 2)
			delta -= SP [j];
		else // competitionMeanField == 1
		  if (WP [i][j]>0)
			delta -= SP [j];
	   }
	  increment += delta * betaP_i[i];
	 }
	 else	
	 {
	  for (j=0; j<N_P; j++)
	  {
	   suma = 0.;
	   if (i!=j)
		for	(k=0; k<N_A; k++)
			if (K [i][k] && K [j][k])
				suma += SA[k] * OmegaPijk [i][j][k]; 
	   delta -= suma * SP[j];
	  }
	  for	(k=0; k<N_A; k++)
		if (K[i][k])
			denominator += SA[k];
	  increment += delta * betaP_i[i] / denominator;
	 }
    } // if (asymmetric == 0)
    else
     if (asymmetric == 1)
     {
      delta = 0.;
      for (j=0; j<N_P; j++)
	  {
	    suma = 0.;
	    if (i!=j)
	     for	(k=0; k<N_A; k++)
			if (K [i][k] && K [j][k])
				suma += SA[k] * OmegaPijk [i][j][k];
	    delta -= suma * SP[j];
	  }
	  for	(k=0; k<N_A; k++)
	  if (K[i][k])
			denominator += SA[k];
	  increment += delta * betaP_i[i] * 0.5 / denominator;
      delta = 0.;
	  for (j=0; j<N_P; j++)
	   if (i!=j)
			delta -= SP [j];
	  increment += delta * betaP_i[i] * 0.5;
     } // if (asymmetric == 1)
     else // i.e., asymmetric == 2
     {
      delta = 0.;
	  for (j=0; j<N_P; j++)
	   if (i!=j)
			delta -= SP [j];
	  increment += delta * betaP_i[i];
     } // if (asymmetric == 0)


    // intra-specie competition term:
	increment -= SP [i] * betaP_ii_0;

    // mutualism term:
	delta = 0.;
	for	(k=0; k<N_A; k++)
		delta += K [i][k] * SA[k] * PhiPik [i][k];
	increment += delta *  gammaP_0 / (1. + hP  * gammaP_0 *  delta);


	// common for all the contributions:
	increment = increment * SP[i];
	increment = increment * delta_t;
	if (SP[i] + increment > 0.)
		return SP[i] + increment;
	return 0.;
}


float animalStep (int i)
{
    int j, k, l;
    float denominator = 0.; 
	float suma;
    // growing term
	float increment = alphaA [i];

    // cross-species competition term:
	float delta = 0.;
	if (competitionMeanField)
	{
	 for (j=0; j<N_A; j++)
	  if (i!=j)
	  {
	    if (competitionMeanField == 2)
			delta -= SA [j];
		else // competitionMeanField == 1
		  if (WA [i][j]>0)
			delta -= SA [j];
	  }
	 increment += delta * betaA_i[i];
	}
	else	
	{
	 for (j=0; j<N_A; j++)
	 {
	  suma = 0.;
	  if (i!=j)
		for	(k=0; k<N_P; k++)
			if (K [k][i] && K [k][j])
				suma += SP[k] * OmegaAkli [i][j][k];; 
	  delta -= suma * SA[j];
	 }
	 for	(k=0; k<N_P; k++)
		if (K[k][i])
			denominator += SP[k];
	 increment += delta * betaA_i[i] / denominator;
	}
    // intra-specie competition term:
	increment -= SA [i] * betaA_ii_0;

    // mutualism term:
	delta = 0.;
	for	(k=0; k<N_P; k++)
		delta += K [k][i] * SP[k] * PhiAki [i][k];
	increment += delta *  gammaA_0 / (1. + hA  * gammaA_0 *  delta);

	// common for all the contributions:
	increment = increment * SA[i];
	increment = increment * delta_t;
	if (SA[i] + increment > 0.)
		return SA[i] + increment;
	return 0.;
}


void readBipartiteMatrix (char* networkName, int nullNumber) // read bipartite matrix file
{
	FILE *file1;
	char filename [128], nullString [16];
    sprintf (filename, "matrix%s.dat", networkName);
    if (nullNumber > 0)
        sprintf (nullString, "null%d_", nullNumber);
	
    char c;
    int i, j, k, l;
	int animalIndex=0;
	if ( ( file1 = fopen( filename, "r" ) ) == NULL )
	{
		printf ("Error: can't read bipartite_matrix file\n");
		exit (-3);
	}

   	N_A=N_P=0;
	for	(l=0; l<maximumN_A; l++)
		kA [l] = 0;
	for	(j=0; j<maximumN_P; j++)
		kP [j] = 0;

   	do
   	{
   		if (fscanf(file1, "%c", &c) == EOF)
   			break;
   		if (c==10)
   		{
   			N_P++;
   			animalIndex = 0;
			if (N_P  >= maximumN_P)
			{
				printf ("Error: maximimN_P is too small\n");
				exit (-5);
			}
   		}
   		else
   			if (c==48||c==49)
   			{
   				K [N_P][animalIndex++] = c-48;
   				if (N_P == 0)
   					N_A++;
   				if (N_A  >= maximumN_A)
				{
					printf ("Error: maximumN_A is too small\n");
					exit (-4);
				}
   			}
   	} while (1);
    fclose (file1);

	for	(j=0; j<N_P; j++)
		for	(l=0; l<N_A; l++)
			if (K [j][l])
			{
   					kP [j]++;
   					kA [l]++;
   			}

	fscanf(file1, "%c",&c);
	projectionMatrices ();

    if (phenology)
    {
     // PhiPik:
     if (nullNumber == 0)
      sprintf (filename, "ov%srows.res", networkName);
     else
      sprintf (filename, "%sov%srows.res", nullString, networkName);
     
	 if ( ( file1 = fopen( filename, "r" ) ) == NULL )
	 {
		printf ("Error: can't read ov_rows.res file \n");
		exit (-3);
	 }
     for(i = 0; i < N_P; i++)
        for (k = 0 ; k < N_A; k++)
            fscanf(file1,"%f",&PhiPik[i][k]);
     fclose (file1);



     // PhiAki:

     if (nullNumber == 0)
      sprintf (filename, "ov%scols.res", networkName);
     else
      sprintf (filename, "%sov%scols.res", nullString, networkName);

	 if ( ( file1 = fopen( filename, "r" ) ) == NULL )
	 {
		printf ("Error: can't read ov_col.res file \n");
		exit (-5);
	 }
     for(k = 0; k < N_A; k++)
        for (i = 0 ; i < N_P; i++)
            fscanf(file1,"%f",&PhiAki[k][i]);
     fclose (file1);



     // OmegaPijk:
     
     if (nullNumber == 0)
      sprintf (filename, "ov%srorows.res", networkName);
     else
      sprintf (filename, "%sov%srorows.res", nullString, networkName);
     
     
	 if ( ( file1 = fopen( filename, "r" ) ) == NULL )
	 {
		printf ("Error: can't read ov_rorows.res file \n");
		exit (-6);
	 }
     for(i = 0; i < N_P; i++)
        for (j = 0 ; j < N_P; j++)
            for(k = 0; k < N_A; k++)
                fscanf(file1,"%f",&OmegaPijk[i][j][k]);
     fclose (file1);


     // OmegaAkli:	

     if (nullNumber == 0)
      sprintf (filename, "ov%scocols.res", networkName);
     else
      sprintf (filename, "%sov%scocols.res", nullString, networkName);

	 if ( ( file1 = fopen( filename, "r" ) ) == NULL )
	 {
		printf ("Error: can't read ov_cocols.res file \n");
		exit (-4);
	 }
     for(k = 0; k < N_A; k++)
        for(l = 0; l < N_A; l++)
            for(i = 0; i < N_P; i++)
                fscanf(file1,"%f",&OmegaAkli[k][l][i]);
     fclose (file1);
    } // if (phenology)
    else
    {
     // PhiPik:
     for(i = 0; i < N_P; i++)
        for (k = 0 ; k < N_A; k++)
            PhiPik[i][k] = 1.;
     // PhiAki:
     for(k = 0; k < N_A; k++)
        for (i = 0 ; i < N_P; i++)
            PhiAki[k][i] = 1.;
     // OmegaPijk:
     for(i = 0; i < N_P; i++)
        for (j = 0 ; j < N_P; j++)
            for(k = 0; k < N_A; k++)
                OmegaPijk[i][j][k] = 1.;
     // OmegaAkli:	
     for(k = 0; k < N_A; k++)
        for(l = 0; l < N_A; l++)
            for(i = 0; i < N_P; i++)
                OmegaAkli[k][l][i] = 1.;
    } // else // if (phenology)

}


void nullModel () // randomize the network according to a probabilistic NULL MODEL
{
    int j, l;

	for	(j=0; j<N_P; j++)
		for	(l=0; l<N_A; l++)
			K[j][l] = (randomFloat() < 0.5*(kP [j]/N_P+kA [l]/N_A)) ? 1: 0;

	for	(j=0; j<N_P; j++)
		kP [j] = 0;
	for	(l=0; l<N_A; l++)
		kA [l] = 0;

	for	(j=0; j<N_P; j++)
		for	(l=0; l<N_A; l++)
			if (K [j][l])
			{
 				kP [j]++;
   				kA [l]++;
   			}
	projectionMatrices ();
}




void printAbundancesVector (float nestedness, float nestednessNODF, int diversity, int networkNumber) // it writes into file
{
	int i;
	FILE *file1;
	char fileName [256];

    /*	
	// with betaP_ii_0 and betaA_ii_0 as fields:
	sprintf (fileName, "N_P%dN_A%d_alphaP_%f_%f_alphaA_%f_%f_hP_%f_hA_%f_asymmetric%dphenology%d.dat",
	         N_P, N_A, min_alphaP, max_alphaP, min_alphaA, max_alphaA, hP, hA, asymmetric, phenology);
	if ( ( file1 = fopen( fileName, "a" ) ) == NULL )
		exit (-2);
	fprintf (file1, "%f %f %f %f %f %f %f %d %d ", betaP_0, betaA_0, gammaP_0, gammaA_0, betaP_ii_0, betaA_ii_0, nestedness, diversity, networkNumber); 
	fprintf (file1, "\n");
	fclose (file1);
    */

	// with betaP_ii_0 and betaA_ii_0 into the file name:
	sprintf (fileName, "N_P%dN_A%d_betaP_ii_0%f_betaA_ii_0%f_alphaP_%f_%f_alphaA_%f_%f_hP_%f_hA_%f_asymmetric%dphenology%d.dat",
	         N_P, N_A, betaP_ii_0, betaA_ii_0, min_alphaP, max_alphaP, min_alphaA, max_alphaA, hP, hA, asymmetric, phenology);
	if ( ( file1 = fopen( fileName, "a" ) ) == NULL )
		exit (-2);
	fprintf (file1, "%f %f %f %f %f %d %f %d ", betaP_0, betaA_0, gammaP_0, gammaA_0, nestedness, diversity, nestednessNODF, networkNumber);
	for (i=0; i<N_P; i++)
		fprintf (file1, "%4.10f ", SP [i]); // abundance vector (plants)
	for (i=0; i<N_A; i++)
		fprintf (file1, "%4.10f ", SA [i]); // abundance vector (animals)
	fprintf (file1, "\n");
	fclose (file1);
    /*
	// with betaP_0, betaA_0, gammaP_0, gammaA_0, betaP_ii_0 and betaA_ii_0 into the file name:
	sprintf (fileName, "N_P%dN_A%d_betaP_0%f_betaA_0%f_gammaP_0%f_gammaA_0%f_betaP_ii_0%f_betaA_ii_0%f_alphaP_%f_%f_alphaA_%f_%f_hP_%f_hA_%f-asymmetric%dphenology%d.dat",
	         N_P, N_A, betaP_0, betaA_0, gammaP_0, gammaA_0, betaP_ii_0, betaA_ii_0, min_alphaP, max_alphaP, min_alphaA, max_alphaA, hP, hA,asymmetric, phenology);
	if ( ( file1 = fopen( fileName, "a" ) ) == NULL )
		exit (-2);
	fprintf (file1, "%f %d %d %f ", nestedness, diversity, networkNumber, nestednessNODF);
	for (i=0; i<N_P; i++)
		fprintf (file1, "%4.10f ", SP [i]); // abundance vector (plants)
	for (i=0; i<N_A; i++)
		fprintf (file1, "%4.10f ", SA [i]); // abundance vector (animals)
	fprintf (file1, "\n");
	fclose (file1);
	*/
}

void printAverage (float bd)
{
	FILE *file1;
	char fileName [256];
	sprintf (fileName, "NullModelPhenology_N_P%dN_A%d_betaP_ii_0%f_betaA_ii_0%f_alphaP_%f_%f_alphaA_%f_%f_hP_%f_hA_%f_asymmetric%dphenology%d.dat",
	         N_P, N_A, betaP_ii_0, betaA_ii_0, min_alphaP, max_alphaP, min_alphaA, max_alphaA, hP, hA, asymmetric, phenology);
	if ( ( file1 = fopen( fileName, "a" ) ) == NULL )
		exit (-2);
	fprintf (file1, "%f %f %f\n", betaP_0, gammaP_0, bd);
	fclose (file1);
}


void printProjectedMatrices (char* matrixFile, int re) // parameter: argv [1]
{
	int i, j;
	FILE *file1;
	char fileName [256];

	// write projected matrix of plants:
	sprintf (fileName, "projected_matrix_plants_%s_rewiring_%d.dat", matrixFile, re);
	if ( ( file1 = fopen( fileName, "w" ) ) == NULL )
		exit (-10);
	for	(i=0; i<N_P; i++)
	{
		for	(j=0; j<N_P; j++)
			fprintf (file1, "%2d ",WP [i][j]);
		fprintf (file1, "\n");
	}
	fclose (file1);

	// write projected matrix of animals:
	sprintf (fileName, "projected_matrix_animals_%s_rewiring_%d.dat", matrixFile, re);
	if ( ( file1 = fopen( fileName, "w" ) ) == NULL )
		exit (-11);
	
	for	(i=0; i<N_A; i++)
	{
		for	(j=0; j<N_A; j++)
			fprintf (file1, "%2d ",WA [i][j]);
		fprintf (file1, "\n");
	}
	fclose (file1);
}

void equationDynamics ()  // differential equations integration
{
	int i, step;
	float localSA[N_A], localSP[N_P];
	for (step=0; 1.*step*delta_t < 1.*transient; step++)
	{
		for (i=0; i<N_A; i++)
			localSA[i] = animalStep(i);
		for (i=0; i<N_P; i++)
			localSP[i]=plantStep (i);
		for (i=0; i<N_A; i++)
			SA [i] = localSA[i];
		for (i=0; i<N_P; i++)
			SP[i] = localSP[i];
	}
}

int biodiversity () // it returns the total number of species
{
	int i;
	int b = 0;
	for (i=0; i<N_P; i++)
		if (SP [i] > 0.0000001)
			b++;
	for (i=0; i<N_A; i++)
		if (SA [i] > 0.0000001)
			b++;
	return b;
}


struct rewiringKey rewiring ()
{
	struct rewiringKey r;
	int i;
	char ok;
	r.plant = rand() % 2; // 1 for plants, 0 for animals
	if (r.plant)
	{
		do
		{
			ok = 0;
			do { r.who = rand () % N_P;} while (kP[r.who] == N_A || kP[r.who] == 0); // error in fully connected networks!
			do { r.withWho = rand () % N_P;} while (r.who == r.withWho || kP[r.withWho] == N_A || kP[r.withWho] == 0);  // error in fully connected networks!
			for (i=0; i<N_A; i++)
				if (K[r.who][i] && !K[r.withWho][i]) 
				{
					ok = 1;
					break;
				}
			if (ok==1)
			  for (i=0; i<N_A; i++)
				if (!K[r.who][i] && K[r.withWho][i]) 
				{
					ok = 2;
					break;
				}
		} while(ok<2);
		do { r.unplug = rand() % N_A; } while (!K[r.who][r.unplug] || K[r.withWho][r.unplug]); //
		do { r.plug = rand() % N_A; } while  (K[r.who][r.plug] || !K[r.withWho][r.plug]); //
		K[r.who][r.unplug] = 0;
		K[r.who][r.plug] = 1;
		K[r.withWho][r.unplug] = 1;
		K[r.withWho][r.plug] = 0;
	}
	else
	{
		do
		{
			ok = 0;
			do { r.who = rand () % N_A;} while (kA[r.who] == N_P || kA[r.who] == 0); // error in fully connected networks!
			do { r.withWho = rand () % N_A;} while (r.who == r.withWho || kA[r.withWho] == N_P || kA[r.withWho] == 0);  // error in fully connected networks!
			for (i=0; i<N_P; i++)
				if (K[i][r.who] && !K[i][r.withWho]) 
				{
					ok = 1;
					break;
				}
			if (ok==1)
			  for (i=0; i<N_P; i++)
				if (!K[i][r.who] && K[i][r.withWho]) 
				{
					ok = 2;
					break;
				}
		} while(ok<2);
		do { r.unplug = rand() % N_P; } while (!K[r.unplug][r.who] || K[r.unplug][r.withWho]); //
		do { r.plug = rand() % N_P; } while  (K[r.plug][r.who] || !K[r.plug][r.withWho]); //
		K[r.unplug][r.who] = 0;
		K[r.plug][r.who] = 1;
		K[r.unplug][r.withWho] = 1;
		K[r.withWho][r.plug] = 0;
	}

	/* // old algorithm, it doesn't preserve k_i
	if (r.plant)
	{
		do { r.who = rand () % N_P;} while (kP[r.who] == N_A || kP[r.who] == 0); // error in fully connected networks!
		do { r.unplug = rand() % N_A; } while (!K[r.who][r.unplug]); //
		do { r.plug = rand() % N_A; } while (K[r.who][r.plug]); //
		K[r.who][r.unplug] = 0;
		K[r.who][r.plug] = 1;
        kA[r.unplug]--;
		kA[r.plug]++;        		
	}
	else
	{
		do { r.who = rand () % N_A;} while (kA[r.who] == N_P || kA[r.who] == 0); // error in fully connected networks!
		do { r.unplug = rand() % N_P; } while (!K[r.unplug][r.who]); //
		do { r.plug = rand() % N_P; } while (K[r.plug][r.who]); //
		K[r.unplug][r.who] = 0;
		K[r.plug][r.who] = 1;
        kP[r.unplug]--;
		kP[r.plug]++;        		
	}
	*/
	projectionMatrices ();
	return r;
}

void undoRewiring (struct rewiringKey r)
{
	if (r.plant)
	{
		K[r.who][r.unplug] = 1;
		K[r.who][r.plug] = 0;
		K[r.withWho][r.unplug] = 0;
		K[r.withWho][r.plug] = 1;
	}
	else
	{
		K[r.unplug][r.who] = 1;
		K[r.plug][r.who] = 0;
		K[r.unplug][r.withWho] = 0;
		K[r.withWho][r.plug] = 1;
	}

    /* old algorithm, it doesn't preserve k_i
	if (r.plant)
	{
		K[r.who][r.unplug] = 1;
		K[r.who][r.plug] = 0;
        kA[r.unplug]++;
		kA[r.plug]--;        		
	}
	else
	{
		K[r.unplug][r.who] = 1;
		K[r.plug][r.who] = 0;
        kP[r.unplug]++;
		kP[r.plug]--;        		
	}
	projectionMatrices ();
	*/
}


void projectionMatrices ()
{
	int i, j, k, l;
	for	(i=0; i<N_A; i++)
		for	(j=0; j<N_A; j++)
			WA [i][j] = 0;
	for	(i=0; i<N_P; i++)
		for	(j=0; j<N_P; j++)
			WP [i][j] = 0;

	//weighted projection matrix (plants)
	for	(i=0; i<N_P; i++)
		for	(j=i+1; j<N_P; j++)
			for	(k=0; k<N_A; k++)
			{
				if (K[i][k] && K[j][k])
					WP [i][j]++;
				WP [j][i] = WP [i][j];
			}
	// weighted projection matrix (animals)
	for	(k=0; k<N_A; k++)
		for	(l=k+1; l<N_A; l++)
			for	(i=0; i<N_P; i++)
			{
				if (K[i][k] && K[i][l])
					WA [k][l]++;
				WA [l][k] = WA [k][l];
			}

}

float largestEigenvalue (int n) // it returns the maximum eigenvalue of the adjacency matrix
{
  double a[n*n];
  double d[n];
  double error_frobenius;
  int it_max, it_num;
  int rot_num;
  double v[n*n];
  int i, j, k;
  int l = 0;
  for (i=0; i<N_P; i++)
  {
  	for (j=0; j<N_P; j++)
  		a [l++] = 0.;
  	for (j=0; j<N_A; j++)
  		a [l++] = K [i][j];
  }
  for (i=0; i<N_A; i++)
  {
  	for (j=0; j<N_P; j++)
  		a [l++] = K [j][i];
  	for (j=0; j<N_A; j++)
  		a [l++] = 0.;
  }

  it_max = 1000;

  jacobi_eigenvalue ( n, a, it_max, v, d, &it_num, &rot_num );

  printf ( "\n" );
  printf ( "  Number of iterations = %d\n", it_num );
  printf ( "  Number of rotations  = %d\n", rot_num );

  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  printf ( "\n" );
  printf ( "  Frobenius norm error in eigensystem A*V-D*V = %g\n",
    error_frobenius );
    
  for (i=0; i<n; i++)
    printf ("eigenvalue %d = %f \n", i,d[i]);
  return  d[n-1];
}

float NODF () // it computes the nestesned according to NODF measure
{
  float nestednessNODF=0.;
  int i, j;
  for (i=0; i<N_P; i++)
  {
  	for (j=i+1; j<N_P; j++)
  	  if (kP[i]!=kP[j]  && kP[i] > 0 && kP[j]>0)
  	    nestednessNODF += 1.* WP[i][j]/(kP[i]<kP[j]?kP[i]:kP[j]);
  }
  for (i=0; i<N_A; i++)
  {
  	for (j=i+1; j<N_A; j++)
  	  if (kA[i]!=kA[j] && kA[i] > 0 && kA[j]>0)
  	    nestednessNODF += 1.* WA[i][j]/(kA[i]<kA[j]?kA[i]:kA[j]);
  }
  return (nestednessNODF/((N_P*(N_P-1))+(N_A*(N_A-1)))/2);
}

void globalVariablesInitialization (int argc, char **argv)
{
	// the first parameter argv[1] is the bipartite matrix's file name 
    if (argc != 12)
    {
      printf ("\nError: there are not 11 parameters.\n");
      exit (-1);
    }
    
	//min_betaP_0=1.*atoi(argv[2])/100; // plants competition constant (cross-species), minimum value to sweep the interval
	//max_betaP_0=1.*atoi(argv[3])/100; // plants competition constant (cross-species), maximum value to sweep the interval
	//delta_betaP_0=1.*atoi(argv[4])/100; // plants competition constant (cross-species), increment value to sweep the interval
	betaP_0 = 1. * atoi(argv[2])/100; 
	betaP_ii_0=1.*atoi(argv[3])/100; // plants competition constant (intra-specie)
	//min_gammaP_0=1.*atoi(argv[6])/100; // plants mutualism constant, minimu//m value to sweep the interval
	//max_gammaP_0=1.*atoi(argv[7])/100; // plants mutualism constant, maximum value to sweep the interval
	//delta_gammaP_0=1.*atoi(argv[8])/100; // plants mutualism constant, increment value to sweep the interval
	gammaP_0 = 1. * atoi(argv[4])/100; 
	min_alphaP=1.*atoi(argv[5])/100; // plants growing constant
	max_alphaP=1.*atoi(argv[6])/100; // animals growing constant
	hP=1.*atoi(argv[7])/100; // holling term for plants (constant) ; ~(N_A*gammaP_0)?
    rewiringsBeforeFirstNetwork = atoi(argv[8]);
 	competitionMeanField = atoi(argv[9]);
	asymmetric = atoi(argv[10]);
	phenology = atoi(argv[11]);


	// inicialited by symmetry (plants constant = animal constant):
	betaA_ii_0 = betaP_ii_0; // animals competition constant (intra-specie)
	min_alphaA = min_alphaP; // plants growing constant
	max_alphaA = max_alphaP; // animals growing constant
	hA = hP; // holling term for animals (constant) ; ~(N*gammaA_0)?

	betaA_0 = betaP_0;// inicialited by symmetry (plants constant = animal constant)
	gammaA_0 = gammaP_0; // inicialited by symmetry (plants constant = animal constant)


	if (max_alphaP < min_alphaP || max_alphaA < min_alphaA)
    {
      printf ("\nError: in parameters, it can't be:  max_alphaP < min_alphaP || max_alphaA < min_alphaA.\n");
      exit (-1);
    }
	if (competitionMeanField < 0 || competitionMeanField > 2)
    {
      printf ("\nError: in parameters, competitionMeanField must be equal to 0, 1 or 2\n");
      exit (-1);
    }

    
}

float swap () // it makes an effective rewiring
{
	struct rewiringKey r;
	int attempt=0;
	float tmpNODF;
	float oldNODF = NODF ();
	printf ("before swap NODF= %f\n", oldNODF);
	
	if (rewiringsTowardLowerNestedness)
	  do
	  {
		r=rewiring ();
		if (attempt++%1000==0)
			printf ("oldNODF = %f  ;    Swaping, attempt number = %d\n", oldNODF, attempt+1);
		tmpNODF = NODF ();
		if (tmpNODF > oldNODF)
			undoRewiring (r);
		else
			return tmpNODF;
	  } while (1);
	rewiring ();
	printf ("after swap NODF= %f\n", NODF ());
	return NODF ();
}

int main (int argc, char **argv)
{
	int biod;
	int networkNumber, rewiringNumber;
	float nestedness, nestednessNODF;
	float AvgBiodiversity = 0.;

	globalVariablesInitialization (argc, argv);
	srand (getpid()+time(NULL)); // random seed
	readBipartiteMatrix (argv[1], 0);
	nestedness = largestEigenvalue (N_P+N_A);
    nestednessNODF = NODF ();
	for (networkNumber=0; networkNumber<networks ; networkNumber++)
	{
        if (networkNumber > 0)
            readBipartiteMatrix (argv[1], networkNumber);
		initialConditions ();
		equationDynamics ();
		biod = biodiversity ();
        if (networkNumber > 0)
            AvgBiodiversity += biod;
		printAbundancesVector (nestedness, nestednessNODF, biod,  networkNumber);
	}
    if (networks > 1)
     printAverage (AvgBiodiversity / (networks-1));
	return (1);
}
