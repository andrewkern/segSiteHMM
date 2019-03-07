/* segSiteHMM.h
/
/
*/

void usage();
double emAlg(HMM *hmm, gsl_vector *startGammas, int snpNumber,\
	     struct snp data[MAXSNPS], int transitionPowerFlag);
void getParameters(int argc, char *argv[]);
void calculateEmissions(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *gammas,\
			double theta,int snpNumber, struct snp theData[MAXSNPS]);
int **elementArray(gsl_vector *vits, double state, int *pCount);
void scoreElements(HMM *h, gsl_matrix *emisLogs, gsl_matrix *vits, char *outfileName, int transitionPowerFlags);


void initializeEmissions(gsl_matrix *emisLogs, gsl_matrix *emisMatrix, int snpNumber, struct snp theData[MAXSNPS]);
