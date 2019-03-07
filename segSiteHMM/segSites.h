//segSites.h

#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_vector.h> 
double my_h(double q, void * p);
double my_H(gsl_vector *pVector, void * p);
void siteProbMatrix(gsl_matrix *probs, gsl_vector *pVector,  int maxSampleSize, gsl_vector *sampleSizeVector);
double weightedLikLookSites(gsl_vector *pVector, void * p);
gsl_vector *jointMLEst(double *lik, void * p);
double weightedLikLookSitesNState(gsl_vector *pVector, void * p);
void jointMLEstNState(gsl_vector *results, double *lik, void * p);
