// ********************  Functions for segSiteHMM


#include "stdio.h"
#include "math.h"
#include "segSites.h"
#include "../hmm/popGenTools.h"
#include "../hmm/numerical.h"
#include "assert.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf.h>


//core function for probs
double my_h(double q, void * p){
  struct my_f_params * params = (struct my_f_params *) p;
  double beta = params->beta;
  int i = params->i;
  int n = params->n;
	
   return ((1.0 - exp(-2.0 * beta * (1.0 - q))) / (1.0 - exp(-2.0 * beta))) \
     * (2.0 / (q * (1.0 - q))) * gsl_ran_binomial_pdf(i, q, n);
}



/* integral of f over allele freq -- returns Prob(Segsite | mu, beta) */
double my_H(gsl_vector *pVector, void * p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_f_params fParams;
  gsl_function f;
  gsl_integration_workspace * w;
  double result,error;
  
  if (gsl_vector_get(pVector, 0) == 0){
    return gsl_vector_get(pVector, 1) / params->i;
  }
  else{
    w = gsl_integration_workspace_alloc (100000);
    fParams.beta = gsl_vector_get(pVector, 0);
    fParams.i = params->i;
    fParams.n = params->n;
    f.function = &(my_h);
    f.params = &fParams;
    gsl_integration_qags(&f,0.0,1.0,1e-5,1e-5, 100000,w, &result,&error);
    gsl_integration_workspace_free(w);
    return result * gsl_vector_get(pVector,1); 
  }
}

/* siteProbMatrix- fills up a matrix of log siteProbs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
   variable sampleSizes (2 to max), columns represent freqs (0 to max -1), freq 0 is for monomorphic sites */
void siteProbMatrix(gsl_matrix *probs, gsl_vector *pVector,  int maxSampleSize, gsl_vector *sampleSizeVector){
  struct my_snpProb_params FParams;
  int i, j;
  double tot, prob;
  
  //check & init probs
  assert(probs->size1 == maxSampleSize + 1 && probs->size2 == maxSampleSize);
  gsl_matrix_set_zero(probs);
 
  //go through sampleSizes
  for(i = 2; i <= maxSampleSize; i++){
    //have sample size i?
    if (gsl_vector_get(sampleSizeVector, i)){
      //set tot to zero for p{monomorphic}
      tot = 0;
      //go through freqs
      for(j = 1; j < maxSampleSize; j++){
	//calc prob
	prob = 0;
	FParams.i = j;
	FParams.n = i;
	//need to impose constraint that p{SegSite = 0 | theta < 0}
	if(gsl_vector_get(pVector,1) < 0){
	  prob = 0;
	}
	else{
	  prob = my_H(pVector, &FParams);
	}
	gsl_matrix_set(probs, i, j, log(prob));
	tot += prob;
      }
      //set p{monomorphic}
      gsl_matrix_set(probs,i, 0, log(1.0 - tot));
    }
  }
}

double weightedLikLookSites(gsl_vector *pVector, void * p){
  struct my_jointLik_params * params = (struct my_jointLik_params *) p;
  gsl_matrix *probs;
  double lik;
  int i;
	
  //make prob matrix (log space)
  probs = gsl_matrix_alloc(params->maxSampleSize + 1, params->maxSampleSize);
  siteProbMatrix(probs, pVector,  params->maxSampleSize, params->sampleSizeVector);

  //tally likelihood
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)  * gsl_vector_get(params->weights,i);
  }
  gsl_matrix_free(probs);
  return -lik;
}
//  weightedLikLookSites2State -- for use in 2 state HMM
double weightedLikLookSites2State(gsl_vector *pVector, void * p){
  struct my_jointLik_params * params = (struct my_jointLik_params *) p;
  gsl_matrix *probs, *probsN;
  gsl_vector *npVector;
  double lik;
  int i;
  
  //setup neutral params
  npVector = gsl_vector_alloc(2);
  gsl_vector_set(npVector,0,0);
  gsl_vector_set(npVector,1,gsl_vector_get(pVector,1));

  //alloc & make probs and probsN matrix (log space)
  probs = gsl_matrix_alloc(params->maxSampleSize + 1, params->maxSampleSize);
  siteProbMatrix(probs, pVector,  params->maxSampleSize, params->sampleSizeVector);
  probsN = gsl_matrix_alloc(params->maxSampleSize + 1, params->maxSampleSize);
  siteProbMatrix(probsN,npVector,  params->maxSampleSize, params->sampleSizeVector);

  //tally likelihood
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)  * gsl_vector_get(params->weights,i);
    lik += gsl_matrix_get(probsN,params->data[i].n,params->data[i].i)  * (1.0 - gsl_vector_get(params->weights,i));
  }
  gsl_matrix_free(probs);
  gsl_matrix_free(probsN);
  gsl_vector_free(npVector);
  return -lik;
}


gsl_vector *jointMLEst(double *lik, void * p){
  size_t np = 2;
  struct my_jointLik_params * params = (struct my_jointLik_params *) p;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  gsl_vector *ss, *x, *results;

  size_t iter = 0, i;
  int status;
  double size;

  /* Initial vertex size vector */
  ss = gsl_vector_alloc (np);
  //initialize results vector
  results = gsl_vector_alloc (np);
  /* Set all step sizes to 1 */
  gsl_vector_set_all (ss, 0.01);
  
  /* Starting point */
  x = gsl_vector_alloc (np);
  gsl_vector_set (x, 0, gsl_vector_get(params->startingPoint,0));
  gsl_vector_set (x, 1, gsl_vector_get(params->startingPoint,1));
  
  /* Initialize method and iterate */
  minex_func.f = &weightedLikLookSites2State;
  minex_func.n = np;
  minex_func.params = params;
  
  s = gsl_multimin_fminimizer_alloc (T, np);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	break;
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-3);
    }
  while (status == GSL_CONTINUE && iter < 500);
  for (i = 0; i < np; i++){
    gsl_vector_set(results,i,gsl_vector_get (s->x, i));
  }
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  return(results);
}

//  weightedLikLookSitesNState -- for use in n state HMM
// pVector now consists of n parameters, not just two
// pVector = [theta_hat, beta_hat_1, beta_hat_2,
double weightedLikLookSitesNState(gsl_vector *pVector, void * p){
  struct my_jointLik_params * params = (struct my_jointLik_params *) p;
  gsl_matrix *probs;
  gsl_vector *currentVector;
  double lik;
  int i,j;
  
  //setup neutral params
  currentVector = gsl_vector_alloc(params->nstates);
  gsl_vector_set(currentVector,0,0);
  gsl_vector_set(currentVector,1,gsl_vector_get(pVector,0));

  //make prob matrix (log space)
  probs = gsl_matrix_alloc(params->maxSampleSize + 1, params->maxSampleSize);
  siteProbMatrix(probs, currentVector,  params->maxSampleSize, params->sampleSizeVector);

  //tally likelihood due to neutral
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)  * gsl_matrix_get(params->posts,0,i);
  }
  //go through states tallying lik
  for (i = 1; i < params->nstates; i++){
    gsl_vector_set(currentVector,0,gsl_vector_get(pVector,i));
    gsl_vector_set(currentVector,1,gsl_vector_get(pVector,0));
    siteProbMatrix(probs,currentVector,  params->maxSampleSize, params->sampleSizeVector);
    for(j = 0; j < params->snpNumber; j++){
      lik += gsl_matrix_get(probs,params->data[j].n,params->data[j].i)  * gsl_matrix_get(params->posts,i,j);
    }
  }
  gsl_matrix_free(probs);
  gsl_vector_free(currentVector);
  return -lik;
}

//N state version
void jointMLEstNState(gsl_vector *results, double *lik, void * p){
  struct my_jointLik_params * params = (struct my_jointLik_params *) p;
  size_t np = params->nstates;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  gsl_vector *ss, *x;

  size_t iter = 0, i;
  int status,j;
  double size;
  
  //check results vector
  assert(results->size == np);
  
  /* Initial vertex size vector */
  ss = gsl_vector_alloc (np);
  
  /* Set all step sizes to 1 */
  gsl_vector_set_all (ss, 0.01);
  
  /* Starting point */
  x = gsl_vector_alloc (np);
  for(j = 0; j < np; j++){
    gsl_vector_set (x, j, gsl_vector_get(params->startingPoint,j));
  }
  
  /* Initialize method and iterate */
  minex_func.f = &weightedLikLookSitesNState;
  minex_func.n = np;
  minex_func.params = params;
  
  s = gsl_multimin_fminimizer_alloc (T, np);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	break;
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-5);
    }
  while (status == GSL_CONTINUE && iter < 500);
  for (i = 0; i < np; i++){
    gsl_vector_set(results,i,gsl_vector_get (s->x, i));
  }
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

