/* popGenTools.h
/
/
/ Andrew Kern
/
/
*/


#define MAXSNPS 1000000
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct snp{ /*SNP struct definition- meant to hold position, derived freq, sample size */
  int pos, i, n, ascSize;
};
	
struct my_F_params{
  int i, n, ascSize;
};

struct my_snpProb_params{
  int i, n, ascSize;
  double denom;
};

struct my_f_params{
  int i,n, ascSize;
  double beta;
};

struct my_lik_params{
  int snpNumber, maxSampleSize, ascSize;
  gsl_vector *weights, *sampleSizeVector, *betas;
  struct snp *data;
  gsl_matrix *sfsBools;

};
struct my_likCI_params{
  int snpNumber, maxSampleSize;
  gsl_vector *weights, *sampleSizeVector;
  struct snp *data;
  gsl_matrix *sfsBools;
  double lMax, logUnits, beta_hat;
};


struct my_jointLik_params{
  int snpNumber, maxSampleSize, ascSize, nstates;
  gsl_vector *weights, *sampleSizeVector, *betas, *startingPoint;
  struct snp *data;
  gsl_matrix *sfsBools, *posts;

};

double my_f(double q, void * params);
double my_F(double beta, void * params);
double my_F2(double beta, void * p);
double my_lik(double beta, void * params);
double sfsLikBetaVector(gsl_vector *beta, void * p);
double weightedLik(double beta, void * params);
double weightedLikLook(double beta, void * params);
double weightedLikLookCI(double beta, void * params);
double likWrap(double beta);
double snpProb(double beta, void * params);
double snpProbDenomLookup(double beta, void * params);
double snpProbDenom(double beta, int sampleSize);
gsl_vector *makeSnpProbDenomVector(double beta, int maxSampleSize, gsl_vector *sampleSizeVec);
double ml_est(double * lik_beta_hat);
double weighted_ml_est(double * lik_beta_hat, void * p);
double weighted_ml_est_lookup(double * lik_beta_hat, void * p);
double weighted_ml_est_lookup_errReport(double * lik_beta_hat, void * p);
double weighted_ml_CILower_lookup(double * lik_beta_hat, void * p);
double weighted_ml_CIUpper_lookup(double * lik_beta_hat, void * p);
double wlikWrap(double beta, void * p);
double  wlikWrapLook(double beta, void * p);
double  wlikCIWrapLook(double beta, void * p);
gsl_vector *sampleSizeVector(struct snp data[], int snpNumber, int maxSampleSize);
gsl_matrix *summarizeSFS(int maxSampleSize, struct snp data[], int snpNumber);
gsl_matrix *summarizeSFSBool(int maxSampleSize, struct snp data[], int snpNumber);
gsl_matrix *snpProbMatrix(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools);
gsl_matrix *snpProbMatrixNotLog(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools);
gsl_matrix *snpProbMatrixNotLogFull(double beta,  int maxSampleSize,gsl_vector *sampleSizeVector);
gsl_vector *snpProbVectorNotLog(double beta,  int sampleSize);
int simulateSiteFreq(int sampleSize, double alpha, void *r);
double simpleAscertain(double sampleFreq, void *p);
double probAscertainmentGivenModel(double beta, void *p);
gsl_matrix *snpAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools, int ascSize);
gsl_matrix *estimateAscSFS(int maxSampleSize, gsl_vector *sampleSizeVector,  gsl_matrix *sfsBools, int ascSize, gsl_matrix *sfsSummary);
double weightedLikLookAsc(double beta, void * p);
double  wlikWrapLookAsc(double beta, void * p);
double weighted_ml_est_lookup_asc(double * lik_beta_hat, void * p);
double probAscertainmentGivenModelLookup(double beta, int sampSize, gsl_matrix *snpProbs, gsl_matrix *ascProbs);
double probAscertainmentGivenModelHemiLookup(double beta, int sampSize, gsl_vector *snpProbs, gsl_matrix *ascProbs);
double outgroupAscertain(double sampleFreq, void * p);
gsl_matrix *snpOutgroupAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, int ascSize);
double weightedLikLookOutgroupAsc(double beta, void * p);
double  wlikWrapLookOutgroupAsc(double beta, void * p);
double weighted_ml_est_lookup_outgroup_asc(double * lik_beta_hat, void * p);
double probOutgroupAscertainmentGivenModel(double beta, void *p);
double sfsLikBetaVectorOutgroupAsc(gsl_vector *betas, void * p);

