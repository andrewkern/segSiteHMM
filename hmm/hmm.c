/* hmm.c v0.0 
 
 	Andrew D. Kern	8/24/05
*/


#include "hmm.h"
#include "adkGSL.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>



/* returns a pointer to a new instance of an HMM. 
Currently an HMM is defined by a transition matrix (pMatrix), 
and a vector representing the starting probabilities.
use of initHMM() fills the log score matrices. 	*/

HMM* newHMM(gsl_matrix *pMatrix, gsl_vector *piStart, gsl_vector *piEnd, gsl_vector *otherData ) {
  HMM *hmm = (HMM*)malloc(sizeof(HMM));
  int i;

  hmm->transition_matrix = pMatrix;
  hmm->piStart = piStart;
  hmm->piEnd = piEnd;
  hmm->nstates = pMatrix->size1;
  hmm->transition_score_matrix = NULL;
  hmm->piStart_scores = NULL;
  hmm->piEnd_scores = NULL;
  hmm->otherData = otherData;

  /* if piStart are NULL, make them uniform */
  if (piStart == NULL) {
    hmm->piStart = gsl_vector_alloc(pMatrix->size1);
    for (i = 0; i < pMatrix->size1; i++) 
      gsl_vector_set(hmm->piStart, i, 1.0/pMatrix->size1);
  }
  if (piEnd == NULL) {
    hmm->piEnd = gsl_vector_alloc(pMatrix->size1);
    for (i = 0; i < pMatrix->size1; i++) 
      gsl_vector_set(hmm->piEnd, i, 1.0/pMatrix->size1);
  }

  return hmm;
}

/* initialize an HMM, this makes sure that things are kosher,
it checks that the row sums of the pMatrix == 1, that the sum
of the piStart vector == 1. then it initializes scores */

void initHMM(HMM *hmm){
	gsl_matrix *pMatrixLogs;
	gsl_vector *piStartLogs, *piEndLogs;
	int i, j;
	double piEndSum, piStartSum, prob;
	
	//check piStartSum == 1
	piStartSum = gsl_vector_sum(hmm->piStart, hmm->nstates);
	assert(piStartSum == 1.0);
	
	//make log transformed start probs
	piStartLogs = gsl_vector_alloc(hmm->nstates);
	for (i = 0; i < hmm->nstates; i++){
		prob = gsl_vector_get(hmm->piStart, i);
		gsl_vector_set(piStartLogs, i, log(prob));
		}
	
	//set hmm piEnd_scores
	hmm->piStart_scores = piStartLogs;
	
	//check piEndSum == 1
	piEndSum = gsl_vector_sum(hmm->piEnd, hmm->nstates);
	assert(piEndSum == 1.0);
	
	//make log transformed start probs
	piEndLogs = gsl_vector_alloc(hmm->nstates);
	for (i = 0; i < hmm->nstates; i++){
		prob = gsl_vector_get(hmm->piEnd, i);
		gsl_vector_set(piEndLogs, i, log(prob));
		}
	
	//set hmm piEnd_scores
	hmm->piEnd_scores = piEndLogs;
	
	//check row sums for markov matrix
	for(i = 0; i < hmm->nstates; i++){
		assert(gsl_matrix_row_sum(hmm->transition_matrix, i, hmm->nstates) == 1.0);
		}
	//make log transformed transition matrix
	pMatrixLogs = gsl_matrix_alloc(hmm->nstates, hmm->nstates);
	for (i = 0; i < hmm->nstates; i++){
		for (j = 0; j < hmm->nstates; j++){
			prob = gsl_matrix_get(hmm->transition_matrix, i, j);
			gsl_matrix_set(pMatrixLogs, i, j, log(prob));
			}
		}
		
	//set hmm transition_score_matrix	
	hmm->transition_score_matrix = pMatrixLogs;
}
	

/* frees an entire hmm */
void freeHMM(HMM *hmm){
	
	if (hmm->transition_matrix)
		gsl_matrix_free(hmm->transition_matrix);
	if (hmm->transition_score_matrix)
		gsl_matrix_free(hmm->transition_score_matrix);
	if (hmm->piStart)
		gsl_vector_free(hmm->piStart);
	if (hmm->piStart_scores)
		gsl_vector_free(hmm->piStart_scores);
	if (hmm->piEnd)
		gsl_vector_free(hmm->piEnd);
	if (hmm->piEnd_scores)
		gsl_vector_free(hmm->piEnd_scores);
	if (hmm->otherData)
		gsl_vector_free(hmm->otherData);
	free(hmm);
}

/*output hmm values: nstates, transition matrix, piStarts */
void printHMM(HMM *hmm){
  printf("hmm output\nstates: %d \ntransition matrix:\n", hmm->nstates);
  gsl_matrix_fprintf(stdout, hmm->transition_matrix, "%f");
  printf("piStart vector:\n");
  gsl_vector_fprintf(stdout, hmm->piStart, "%f");
  printf("piEnd vector:\n");
  gsl_vector_fprintf(stdout, hmm->piEnd, "%f");
  if (hmm->transition_score_matrix != NULL){
    printf("transition matrix (logs):\n");
    gsl_matrix_fprintf(stdout, hmm->transition_score_matrix, "%f");
  }
  if (hmm->piStart_scores != NULL){
    printf("piStart vector (logs):\n");
    gsl_vector_fprintf(stdout, hmm->piStart_scores, "%f");
  } 
  if (hmm->piEnd_scores != NULL){
    printf("piEnd vector (logs):\n");
    gsl_vector_fprintf(stdout, hmm->piEnd_scores, "%f");
  } 
  if (hmm->otherData != NULL){
    //    printf("other data:\n");
    //gsl_vector_fprintf(stdout, hmm->otherData, "%f");
  }
  
}

/*output transition matrix of hmm values to file */
void printTransitions(HMM *h, char *outfileName){
  FILE *outfile;
  int i,j;

   outfile = fopen(outfileName, "w");
   if (outfile == NULL){
     fprintf(stderr,"Error opening outfile! ARRRRR!!!!\n");
     exit(1);
   }
   for(i = 0; i < h->nstates; i++){
     for(j = 0; j < h->nstates; j++){
       fprintf(outfile,"%f\t",gsl_matrix_get(h->transition_matrix,i,j));
     }
     fprintf(outfile,"\n");
   }
   fclose(outfile);
}

/* this normalizes the transition matrix so that rows sum to one */
void hmm_normalize_transitions(HMM *hmm){
  int i, j;
  double rowSum, tempVal;

  for(i=0; i < hmm->nstates; i++){
    rowSum = 0;
    for(j=0; j < hmm->nstates; j++){
      rowSum += gsl_matrix_get(hmm->transition_matrix,i,j);
    }
    if(rowSum == 0){
      gsl_matrix_set(hmm->transition_matrix, i, i, 1.0);
    }
    else{
      for(j= 0; j < hmm->nstates;j++){
	tempVal =  gsl_matrix_get(hmm->transition_matrix, i, j) / rowSum;
	gsl_matrix_set(hmm->transition_matrix, i, j, tempVal);
      }
    }
  }
}

void hmm_logify_transitions(HMM *hmm){
  int i, j;
  double prob;
  //set log transition scores
  for (i = 0; i < hmm->nstates; i++){
    for (j = 0; j < hmm->nstates; j++){
      prob = gsl_matrix_get(hmm->transition_matrix, i, j);
      gsl_matrix_set(hmm->transition_score_matrix, i, j, log(prob));
    }
  }
  //set log start and end scores
  for (i = 0; i < hmm->nstates; i++){
    prob = gsl_vector_get(hmm->piStart, i);
    gsl_vector_set(hmm->piStart_scores, i, log(prob));
    prob = gsl_vector_get(hmm->piEnd, i);
    gsl_vector_set(hmm->piEnd_scores, i, log(prob));
  }
}

/* forward algorithm- this fills a matrix (alphas) of Nstates * L observations
with the forward log probabilities. a matrix of emission probs (log!!!) with the same
dimensions as the forward matrix must be allocated and passed to this function.
returns probFor. */

double forwardAlg(HMM *hmm, gsl_matrix *alphas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
	int i, j, power, last;
	double prob, emit, prevSum, pFor;
	gsl_vector_view prevCol;
	gsl_vector *tempProbs;
	gsl_matrix *logTransitions;

	assert(length > 0 && hmm != NULL);
	assert(alphas->size1 == hmm->nstates && alphas->size2 == length);
	assert(emisLogs->size1 == hmm->nstates && emisLogs->size2 == length);

//initialization
	for (i = 0; i < hmm->nstates; i++){
		prob = gsl_vector_get(hmm->piStart_scores, i);
		emit = gsl_matrix_get(emisLogs, i, 0);
		gsl_matrix_set(alphas, i, 0, (prob + emit));
	}

//are we using powers of the transition matrix? if so check to see we have location data
	if (transitionPowerFlag){
		assert(hmm->otherData != NULL);
	//set "last" equal to the first location (stored in otherData) here's a totally non general assumption... FIX ME!!!
		last = gsl_vector_get(hmm->otherData,0);
	}
	else {
		last = 0; /* not used */
	}

//fill em up!
	for (j = 1; j < length; j++){			//iterate over obs
		for (i = 0; i < hmm->nstates; i++){		//iterate over states
			emit = gsl_matrix_get(emisLogs, i, j);
			/* are we taking powers here? */
			if (transitionPowerFlag){
			//get power from data 
				power = gsl_vector_get(hmm->otherData,j) - last;

				//get power of matrix in log scale
				logTransitions = gsl_matrix_power_logs(hmm->transition_matrix, power);
				prevSum = forwardSumPrevious(hmm, alphas, logTransitions, i, j);
				gsl_matrix_free(logTransitions);
			}
			/* no powers */
			else{
				prevSum = forwardSumPrevious(hmm, alphas, hmm->transition_score_matrix, i, j);
			//fprintf (stderr,"forwardAlg state %2d, index %2d, log.prevSum %g, emit %g\n",i,j,prevSum,emit);
			}
			gsl_matrix_set(alphas, i, j, (emit + prevSum));
		}
		if (transitionPowerFlag){
			last = gsl_vector_get(hmm->otherData,j);
		}
	}	
//return P(x)
	tempProbs = gsl_vector_alloc(hmm->nstates);
	prevCol = gsl_matrix_column(alphas,  length - 1);
	for(i = 0; i < hmm->nstates; i++){
		gsl_vector_set(tempProbs, i, gsl_vector_get(&prevCol.vector, i));
	}
	pFor = log_sum(tempProbs);
	gsl_vector_free(tempProbs);	
	return(pFor);
}


/* forward algorithm offset- identical to the normal forward algorithm except starts at an arbitrary
offset with respect to the observations. this is useful for getting probabilities from path subsets */

double forwardAlgOffset(HMM *hmm, gsl_matrix *alphas, gsl_matrix *emisLogs, int length, int offset, int transitionPowerFlag){
  int i, j, power, last;
  double prob, emit, prevSum, pFor;
  gsl_vector_view prevCol;
  gsl_vector *tempProbs;
  gsl_matrix *logTransitions;
	
  assert(length > 0 && hmm != NULL);
  assert(alphas->size1 == hmm->nstates && alphas->size2 == length);
	
  //initialization
  for (i = 0; i < hmm->nstates; i++){
    prob = gsl_vector_get(hmm->piStart_scores, i);
    emit = gsl_matrix_get(emisLogs, i, offset);
    gsl_matrix_set(alphas, i, 0, (prob + emit));
  }
  
  //are we using powers of the transition matrix? if so check to see we have location data
  if (transitionPowerFlag){
    assert(hmm->otherData != NULL);
    //set "last" equal to the first location (stored in otherData) here's a totally non general assumption... FIX ME!!!
    last = gsl_vector_get(hmm->otherData,offset);
  }
  else {
      last = 0; /* not used */
  }

  //fill em up!
  for (j = 1; j < length; j++){			//iterate over obs
    for (i = 0; i < hmm->nstates; i++){		//iterate over states
      emit = gsl_matrix_get(emisLogs, i, j + offset);
      /* are we taking powers here? */
      if (transitionPowerFlag){
	//get power from data 
	power = gsl_vector_get(hmm->otherData,j + offset) - last;

	//get power of matrix in log scale
	logTransitions = gsl_matrix_power_logs(hmm->transition_matrix, power);
	prevSum = forwardSumPrevious(hmm, alphas, logTransitions, i, j);
	gsl_matrix_free(logTransitions);
      }
      /* no powers */
      else{
	prevSum = forwardSumPrevious(hmm, alphas, hmm->transition_score_matrix, i, j);
      }
      gsl_matrix_set(alphas, i, j, (emit + prevSum));
    }
    if (transitionPowerFlag){
	last = gsl_vector_get(hmm->otherData,j+offset);
    }
  }	
  //return P(x)
  tempProbs = gsl_vector_alloc(hmm->nstates);
  prevCol = gsl_matrix_column(alphas,  length - 1);
  for(i = 0; i < hmm->nstates; i++){
    gsl_vector_set(tempProbs, i, gsl_vector_get(&prevCol.vector, i));
  }
  pFor = log_sum(tempProbs);
  gsl_vector_free(tempProbs);	
  return(pFor);
}

/* backward algorithm- this fills a matrix (betas) of Nstates * L observations
with the backward log probabilities. a matrix of emission probs (log!!!) with the same
dimensions as the backward matrix must be allocated and passed to this function.
returns probBack. */

double backwardAlg(HMM *hmm, gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  int i, j, last, power;
  gsl_vector_view nextCol;
  gsl_vector *tempProbs;
  gsl_matrix *logTransitions;
  double pBack;
	
  assert(length > 0 && hmm != NULL);
  assert(betas->size1 == hmm->nstates && betas->size2 == length);
  assert(emisLogs->size1 == hmm->nstates && emisLogs->size2 == length);
	
  //initialization
  for (i = 0; i < hmm->nstates; i++){
    //prob = gsl_vector_get(hmm->piEnd_scores, i);
    //emit = gsl_matrix_get(emisLogs, i, length - 1);
    gsl_matrix_set(betas, i, length - 1, 0);		//initialize to prob = 1.0, log = 0
  }

  //are we using powers of the transition matrix?
  if (transitionPowerFlag){
    assert(hmm->otherData != NULL);
    //set "last" equal to the last location (stored in otherData) here's a totally non general assumption... FIX ME!!!
    last = gsl_vector_get(hmm->otherData,length - 1);
  }
  else {
      last = 0; /* not used */
  }
	
  //fill em up!
  for (j =  length - 2; j >=0; j--){			//iterate over obs
    for (i = 0; i < hmm->nstates; i++){		//iterate over states
      /* are we taking powers here? */
      if (transitionPowerFlag){
	//get power from data
	power = last - gsl_vector_get(hmm->otherData,j);

	//get power of matrix in log scale
	logTransitions = gsl_matrix_power_logs(hmm->transition_matrix, power);
	gsl_matrix_set(betas, i, j, backwardSumNext(hmm, betas,logTransitions,emisLogs, i, j));
	gsl_matrix_free(logTransitions);
      }
      else{
	gsl_matrix_set(betas, i, j, backwardSumNext(hmm, betas,hmm->transition_score_matrix,emisLogs, i, j));
      }
    }
    if (transitionPowerFlag){
	last = gsl_vector_get(hmm->otherData,j);
    }
  }
  //return P(x)
  tempProbs = gsl_vector_alloc(hmm->nstates);
  nextCol = gsl_matrix_column(betas, 0);
  for(i = 0; i < hmm->nstates; i++){
    gsl_vector_set(tempProbs, i, (gsl_vector_get(&nextCol.vector, i) + gsl_matrix_get(emisLogs, i, 0) + gsl_vector_get(hmm->piStart_scores, i)));
  }
  pBack = log_sum(tempProbs); 
  gsl_vector_free(tempProbs);	
  return(pBack);
}


/* posterior probs- this fills a matrix (postProbMat) of nstates * L observations with the 
posterior prob (not log !!) of being in that state at time i given the observations and the model.
calls forwardAlg and backwardAlg, checks that prob of observations determined by
forward and backward are (nearly) equal, and also (nearly) equal the sum of each column
before normalization to postProbMat.
a matrix of emission probs (log!!!) of dimensions nstates * L observations must be supplied. */

void posteriorProbs(HMM *hmm, gsl_matrix *postProbMat, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
	gsl_matrix *alphas, *betas;
	gsl_vector *tempVec;
	int i, j;
	double pFor, pBack, sumProb;
	
	//allocate matrices for forwards and backwards algs.
	alphas = gsl_matrix_alloc(hmm->nstates, length);
	betas = gsl_matrix_alloc(hmm->nstates, length);
	
	//run forwards/backwards
	pFor = forwardAlg(hmm, alphas, emisLogs, length, transitionPowerFlag);
	pBack = backwardAlg(hmm, betas, emisLogs, length, transitionPowerFlag);
	
	//check that Forward and Backward look okay
	//assert(fabs(pFor - pBack) < 0.01);
	
	//do the posterior probability thing- okay sir pull the switch
	tempVec = gsl_vector_alloc(hmm->nstates);
	for(j=0; j < length; j++){
		//first get sum for denominator
		for(i = 0; i < hmm->nstates; i++){
			gsl_vector_set(tempVec, i, (gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j))); 
		}
		sumProb = log_sum(tempVec);
		
		//check that sum looks like pFor
		assert(fabs(pFor - sumProb) < 0.01);
		
		//then calculate marginal and set value (not log)
		for(i = 0; i < hmm->nstates; i++){
			gsl_matrix_set(postProbMat, i, j, exp(gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j) - sumProb));
		}
	}
	//cleanup
	gsl_matrix_free(alphas);
	gsl_matrix_free(betas);
	gsl_vector_free(tempVec);
}

/* posterior probs reduced- this fills a matrix of nstates * L observations with the 
posterior probs as above except this time takes alphas and betas as args */

void posteriorProbsReduced(HMM *hmm, gsl_matrix *postProbMat, double pFor, gsl_matrix *alphas, gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  gsl_vector *tempVec;
  int i, j;
  double  sumProb;
  
  //do the posterior probability thing- okay sir pull the switch
  tempVec = gsl_vector_alloc(hmm->nstates);
  for(j=0; j < length; j++){
    //first get sum for denominator
    for(i = 0; i < hmm->nstates; i++){
      gsl_vector_set(tempVec, i, (gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j))); 
    }
    sumProb = log_sum(tempVec);
		
    //check that sum looks like pFor
      //assert(fabs(pFor - sumProb) < 0.01);
    
    //then calculate marginal and set value (not log)
    for(i = 0; i < hmm->nstates; i++){
      gsl_matrix_set(postProbMat, i, j, exp(gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j) - sumProb));
    }
  }
  //cleanup
   gsl_vector_free(tempVec);
}

/* posterior probs reduced log- this fills a matrix of nstates * L observations with the 
posterior probs as above except this time takes alphas and betas as args */

void posteriorProbsReducedLog(HMM *hmm, gsl_matrix *postProbMat, double pFor, gsl_matrix *alphas, gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  gsl_vector *tempVec;
  int i, j;
  double  sumProb;
  
  //do the posterior probability thing- okay sir pull the switch
  tempVec = gsl_vector_alloc(hmm->nstates);
  for(j=0; j < length; j++){
    //first get sum for denominator
    for(i = 0; i < hmm->nstates; i++){
      gsl_vector_set(tempVec, i, (gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j))); 
    }
    sumProb = log_sum(tempVec);
		
    //check that sum looks like pFor
    //assert(fabs(pFor - sumProb) < 0.01);
    
    //then calculate marginal and set value
    for(i = 0; i < hmm->nstates; i++){
      gsl_matrix_set(postProbMat, i, j, gsl_matrix_get(alphas,i,j) +  gsl_matrix_get(betas,i,j) - sumProb);
    }
  }
  //cleanup
   gsl_vector_free(tempVec);
}


/* meant to computer prob{state_i_t, state_j_t+1 | data, model} -- for restimating transitions, takes as an arg *jpMatArray an
array of gsl_matrix structs representing joint posterior probabilities -- log scale so need log posts */

void jointPosteriorProbsReduced(HMM *hmm, gsl_matrix **jpMatArray, double pFor, gsl_matrix *alphas, gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  
  int i, j, t;
  double  sumProb, val = 0;
  gsl_matrix *powTrans;
  
  //for each obs
  for(t=0; t < length - 1; t++){
    gsl_matrix_set_all(jpMatArray[t],0);
    //iterate i,j states    
    for(i = 0; i < hmm->nstates; i++){
      for(j = 0; j < hmm->nstates; j++){
	//transition stuff
	if (transitionPowerFlag){
	  
	  if (gsl_vector_get(hmm->otherData,t+1) -  gsl_vector_get(hmm->otherData, t) > 1){
	    powTrans = gsl_matrix_power_logs(hmm->transition_matrix, gsl_vector_get(hmm->otherData,t+1) \
					     -  gsl_vector_get(hmm->otherData, t));
	    val = gsl_matrix_get(alphas, i, t) + gsl_matrix_get(powTrans, i, j) \
	      + gsl_matrix_get(emisLogs, j, t+1) + gsl_matrix_get(betas, j, t+1) - pFor;
	    gsl_matrix_free(powTrans);
	  }
	  else{
	    val = gsl_matrix_get(alphas, i, t) + gsl_matrix_get(hmm->transition_score_matrix, i, j) \
	      + gsl_matrix_get(emisLogs, j, t+1) + gsl_matrix_get(betas, j, t+1) - pFor;
	  }
	}
	else{
	  val = gsl_matrix_get(alphas, i, t) + gsl_matrix_get(hmm->transition_score_matrix, i, j) \
		    + gsl_matrix_get(emisLogs, j, t+1) + gsl_matrix_get(betas, j, t+1) - pFor;
	}
	//set i,j value in gsl_matrix at t
	gsl_matrix_set(jpMatArray[t], i, j, val);
      }
    } 
  }
}

/*restimateTransitions-- uses posteriors and joint posteriors to reestimate the transition matrix based on EM */
void reestimateTransitions(HMM *hmm, gsl_matrix *postProbMat, gsl_matrix **jpMatArray, double pFor, gsl_matrix *alphas, \
			   gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  double sum_gamma, sum_cosi, val;
  int i, j, t,k;
  gsl_vector *tmpVec1, *tmpVec2;
  gsl_matrix *tempTrans, *powTrans;
  gsl_vector *rowSums;

  /* this is the joint posterior method- currently not operational */
  /*
  //calculate joints, no reason to do it twice...
  //  jointPosteriorProbsReduced(hmm, jpMatArray, pFor, alphas, betas, emisLogs, length ,transitionPowerFlag);
  //allocate tmpVecs
  tmpVec1 = gsl_vector_alloc(length - 2);
  tmpVec2 = gsl_vector_alloc(length - 2);
  //iterate over states 
  for(i = 0; i < hmm->nstates; i++){
    for(j = 0; j < hmm->nstates; j++){
      gsl_vector_set_all(tmpVec1,0);
      gsl_vector_set_all(tmpVec2,0);
      //sum up numerator and denom
      for(t = 1; t < length - 1; t++){
	//	gsl_vector_set(tmpVec1, t - 1, gsl_matrix_get(postProbMat,i,t));
	//gsl_vector_set(tmpVec2, t - 1, gsl_matrix_get(jpMatArray[t],i,j));
	sum_gamma += gsl_matrix_get(postProbMat,i,t);
	sum_cosi += gsl_matrix_get(jpMatArray[t],i,j);
      }
      
      //gsl_matrix_set(hmm->transition_matrix, i, j, exp(log_sum(tmpVec1) - log_sum(tmpVec2)));
      // hmm_logify_transitions(hmm);
      gsl_matrix_set(hmm->transition_matrix, i, j, sum_cosi / sum_gamma);
       gsl_matrix_set(hmm->transition_score_matrix, i, j, log(sum_cosi / sum_gamma));
    }
  }
  */

  /* older method */
  tempTrans = gsl_matrix_alloc(hmm->nstates,hmm->nstates);
  rowSums = gsl_vector_alloc(hmm->nstates);
  
  gsl_matrix_set_all(tempTrans,0.0);
  gsl_vector_set_all(rowSums, 0.0);
  for(i = 1; i < length - 1; i++){
    for(j = 0; j < hmm->nstates; j++){
      for(k = 0; k < hmm->nstates; k++){
	if (transitionPowerFlag){
	  if (gsl_vector_get(hmm->otherData,i + 1) -  gsl_vector_get(hmm->otherData, i ) > 1){
	    powTrans = gsl_matrix_power_logs(hmm->transition_matrix, \
					     gsl_vector_get(hmm->otherData,i + 1) -  gsl_vector_get(hmm->otherData, i));
	    val = exp(gsl_matrix_get(alphas, j, i) + gsl_matrix_get(powTrans, j, k) \
		      + gsl_matrix_get(emisLogs, k, i+1) + gsl_matrix_get(betas, k, i+1) - pFor);
	    gsl_matrix_free(powTrans);
	  }
	  else{
	    val = exp(gsl_matrix_get(alphas, j, i) + gsl_matrix_get(hmm->transition_score_matrix, j, k) \
		      + gsl_matrix_get(emisLogs, k, i+1) + gsl_matrix_get(betas, k, i+1) - pFor);
	  }
	}
	else{
	    val = exp(gsl_matrix_get(alphas, j, i) + gsl_matrix_get(hmm->transition_score_matrix, j, k) \
		      + gsl_matrix_get(emisLogs, k, i+1) + gsl_matrix_get(betas, k, i+1) - pFor);
	}
      	gsl_matrix_set(tempTrans, j, k, gsl_matrix_get(tempTrans,j,k) + val);
	gsl_vector_set(rowSums, j, gsl_vector_get(rowSums, j) + val);
      }
    }
  }
  
  //update transition scores
  for(j = 0; j < hmm->nstates; j++){
    for(k = 0; k < hmm->nstates; k++){      
      val = gsl_matrix_get(tempTrans, j, k) / gsl_vector_get(rowSums, j);
      gsl_matrix_set(hmm->transition_matrix, j, k, val);
      gsl_matrix_set(hmm->transition_score_matrix, j, k, log(val));
    }
  }
  gsl_matrix_free(tempTrans);
  gsl_vector_free(rowSums);
}

/*restimateTransitions_vanilla-- uses posteriors to reestimate the transition matrix based on EM */
void reestimateTransitions_vanilla(HMM *hmm, double pFor, gsl_matrix *alphas, \
gsl_matrix *betas, gsl_matrix *emisLogs, int length){
	double val;
	int i, j,k;
	gsl_vector *rowSums;
	gsl_matrix *tempTrans;

/* older method */
  	tempTrans = gsl_matrix_alloc(hmm->nstates,hmm->nstates);
	gsl_matrix_set_all(tempTrans,0.0);
	rowSums = gsl_vector_alloc(hmm->nstates);
	gsl_vector_set_all(rowSums, 0.0);
	for(i = 1; i < length - 1; i++){
		for(j = 0; j < hmm->nstates; j++){
			for(k = 0; k < hmm->nstates; k++){
				val = exp(gsl_matrix_get(alphas, j, i) + gsl_matrix_get(hmm->transition_score_matrix, j, k) \
					+ gsl_matrix_get(emisLogs, k, i+1) + gsl_matrix_get(betas, k, i+1) - pFor);

				gsl_matrix_set(tempTrans, j, k, gsl_matrix_get(tempTrans,j,k) + val);
				gsl_vector_set(rowSums, j, gsl_vector_get(rowSums, j) + val);
			}
		}
	}

//update transition scores
	for(j = 0; j < hmm->nstates; j++){
		for(k = 0; k < hmm->nstates; k++){      
			val = gsl_matrix_get(tempTrans, j, k) / gsl_vector_get(rowSums, j);
			gsl_matrix_set(hmm->transition_matrix, j, k, val);
			gsl_matrix_set(hmm->transition_score_matrix, j, k, log(val));
		}
	}
	gsl_matrix_free(tempTrans);
	gsl_vector_free(rowSums);
}

/* posteriorDecode-- fills a vector (posteriorPath) of length L observations
with the most probably state sequence based on posterior probabilities */
void posteriorDecode(HMM *hmm, gsl_vector *posteriorPath, gsl_matrix *posts, int length){
	int i, j, maxState;
	double maxProb;
	
	assert(length > 0 && hmm != NULL);
	assert(posteriorPath->size == length);
	assert(posts->size1 == hmm->nstates && posts->size2 == length);
	
	for(i=0; i < length; i++){
		maxProb = 0;
		maxState = 0;
		for(j=0;j < hmm->nstates;j++){
			if (gsl_matrix_get(posts, j, i) > maxProb){
				maxProb = gsl_matrix_get(posts, j, i);
				maxState = j;
			}
		}
		gsl_vector_set(posteriorPath,i,maxState);
	}
	
}

/* viterbi algorithm- this fills a vector (viterbiPath) of length L observations with
the most probable state sequence. it returns the prob of the most probable path.
a matrix of emission probs (log!!!) of dimensions nstates * L observations must be supplied. */

  double viterbi(HMM *hmm, gsl_vector *viterbiPath, gsl_matrix *emisLogs, int length, int transitionPowerFlag){
  gsl_matrix *deltas, *psi, *logTransitions;
	gsl_vector_view col;
	int i, j, lastState, tmp, last, power;
	double tempProb,prob, emit;
	int point;
	
	assert(length > 0 && hmm != NULL);
	assert(viterbiPath->size == length);
	assert(emisLogs->size1 == hmm->nstates && emisLogs->size2 == length);
	

	//if using transition powers, then check that we have other data
	if (transitionPowerFlag){
	  assert(hmm->otherData != NULL);
	  //initialize last
	  last = gsl_vector_get(hmm->otherData,0);
	} 
	else {
	    last = 0; /* not used */
	}
	
	//allocate matrices for deltas (log prob of most prob path to stateN,obsL) and psis (the backpointers)
	deltas = gsl_matrix_alloc(hmm->nstates, length);
	psi    = gsl_matrix_alloc(hmm->nstates, length);
	
	//initialization
	for (i = 0; i < hmm->nstates; i++){
		prob = gsl_vector_get(hmm->piStart_scores, i);
		emit = gsl_matrix_get(emisLogs, i, 0);
		gsl_matrix_set(deltas, i, 0, (prob + emit));
		gsl_matrix_set(psi, i, 0, 0);
		}

	point = 666;
	//recursion
	for (j = 1; j < length; j++){
	  for(i = 0; i < hmm->nstates; i++){
	    /* transition powers? */
	    if (transitionPowerFlag){
	      power = gsl_vector_get(hmm->otherData,j) - last;
	      logTransitions = gsl_matrix_power_logs(hmm->transition_matrix, power);
	      tempProb = viterbiMaxPrevious(hmm, deltas, logTransitions, i, j, &point);
	      gsl_matrix_free(logTransitions);
	    }
	    else{
	      tempProb =  viterbiMaxPrevious(hmm, deltas, hmm->transition_matrix, i, j, &point);
	    }
	    gsl_matrix_set(deltas, i, j, tempProb + gsl_matrix_get(emisLogs, i, j));
	    gsl_matrix_set(psi, i ,j , point);
	  }
	  if (transitionPowerFlag){
	      last = gsl_vector_get(hmm->otherData,j);
	  }
	}

	//termination
	col = gsl_matrix_column(deltas, length - 1);
	lastState = (int) gsl_vector_max_index(&col.vector);
  	
	//backtrack
	gsl_vector_set(viterbiPath, length - 1, lastState);
	for (j = length - 2; j >= 0; j--){
		tmp = gsl_matrix_get(psi, lastState, j + 1);
		lastState =  tmp;
		gsl_vector_set(viterbiPath, j , (int) lastState);
		}
	
		
	return gsl_vector_max(&col.vector);
	gsl_matrix_free(deltas);
	gsl_matrix_free(psi);
}	
	

/* this returns the sum of the previous states for the forward probabilities */

double forwardSumPrevious(HMM *hmm, gsl_matrix *alphas, gsl_matrix *transition_score_matrix_i, int state, int obsIndex){
	gsl_vector *tempProbs;
	gsl_vector_view prevCol;
	int i;
	double tempProb;
	
	tempProbs = gsl_vector_alloc(hmm->nstates);
	prevCol = gsl_matrix_column(alphas, obsIndex - 1);
	for(i = 0; i < hmm->nstates; i++){
		tempProb = gsl_vector_get(&prevCol.vector, i) + gsl_matrix_get(transition_score_matrix_i, i, state);
		gsl_vector_set(tempProbs, i, tempProb);
		}
	tempProb = log_sum(tempProbs); 
	gsl_vector_free(tempProbs);
	return(tempProb);
}

/* this returns the sum of the "next" states for the backward probabilities */

double backwardSumNext(HMM *hmm, gsl_matrix *betas, gsl_matrix *transition_score_matrix_j,gsl_matrix *emisLogs, int state, int obsIndex){
	gsl_vector *tempProbs;
	gsl_vector_view nextCol;
	int j;
	double tempProb;
	
	tempProbs = gsl_vector_alloc(hmm->nstates);
	nextCol = gsl_matrix_column(betas, obsIndex + 1);
	for(j = 0; j < hmm->nstates; j++){		//sum over begin states...
		tempProb = gsl_vector_get(&nextCol.vector, j) + gsl_matrix_get(transition_score_matrix_j, state, j) + gsl_matrix_get(emisLogs, j,obsIndex + 1) ;
		gsl_vector_set(tempProbs, j, tempProb);
		}
	tempProb = log_sum(tempProbs);
	gsl_vector_free(tempProbs);
	return(tempProb); 
}	

/* this returns the max of the previous states for the viterbi algorithm. 
it also sets the back pointer (point) to that state. */

double viterbiMaxPrevious(HMM *hmm, gsl_matrix *deltas, gsl_matrix *transition_score_matrix_i, int state, int obsIndex, int *point){
	gsl_vector *tempProbs;
	gsl_vector_view prevCol;
	int i;
	double tempProb;
	
	tempProbs = gsl_vector_alloc(hmm->nstates);
	prevCol = gsl_matrix_column(deltas, obsIndex - 1);
	for(i = 0; i < hmm->nstates; i++){
		tempProb = gsl_vector_get(&prevCol.vector, i) + gsl_matrix_get(transition_score_matrix_i, i, state);
		gsl_vector_set(tempProbs, i, tempProb);
		}
	*point = gsl_vector_max_index(tempProbs);
	tempProb = gsl_vector_max(tempProbs);
	gsl_vector_free(tempProbs);
	return(tempProb);
}	
	
/* This computes the log likelihood of an subpath of the data, for a specified set of the states defined
    in the vector stateInclusion (1 for include, 0 for exclude). */
double hmm_subpath_score(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *stateInclusion, int begin, int length, int transitionPowerFlag){
  gsl_matrix *alphas, *real_trans ;
  gsl_vector *real_begin;
  int i, j, stateCount;
  double score;
    
  //initialize some matrices
  alphas = gsl_matrix_alloc(hmm->nstates, length);
  real_trans = gsl_matrix_alloc(hmm->nstates, hmm->nstates);
  real_begin = gsl_vector_alloc(hmm->nstates);

  //count number of states under consideration
  stateCount = 0;
  for(i = 0; i < hmm->nstates; i++){
    stateCount += gsl_vector_get(stateInclusion, i);
  }

  //keep track of original transition matrix and "begin" vector
  gsl_matrix_memcpy(real_trans,hmm->transition_matrix);
  gsl_vector_memcpy(real_begin, hmm->piStart);
  
  //set begin states to uniform distribution over states considered in subset
  for(i = 0; i < hmm->nstates; i++){
    if (gsl_vector_get(stateInclusion,i)){
      gsl_vector_set(hmm->piStart, i, 1.0 / stateCount);
    }
    else{
      gsl_vector_set(hmm->piStart, i, 0.0);
    }
  }

  //renomalize transition matrix to reflect states included.
  //go through transition matrix, shutting off those which aren't in subset
  for(i = 0; i < hmm->nstates; i++){
    for(j = 0; j < hmm->nstates; j++){
      if (gsl_vector_get(stateInclusion,i) == 0 || gsl_vector_get(stateInclusion,j) == 0){
	gsl_matrix_set( hmm->transition_matrix, i, j, 0.0);
      }
    }
  }
 
  //one is the magic number...
  hmm_normalize_transitions(hmm);
  hmm_logify_transitions(hmm);

  //run the forward algorithm to get our answer!
  score = forwardAlgOffset(hmm, alphas, emisLogs, length, begin, transitionPowerFlag);

  //cleanup
  gsl_matrix_memcpy(hmm->transition_matrix, real_trans);
  gsl_vector_memcpy(hmm->piStart, real_begin);
  hmm_logify_transitions(hmm);
  gsl_matrix_free(alphas);
  gsl_matrix_free(real_trans);
  gsl_vector_free(real_begin);
  return(score);
}
  
/* Calculates the log odds score for two competing sets of states, over a subset of the data
   by comparing the ratio of their likelihoods */
double hmm_subpath_log_odds(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *stateInclusion1,
			     gsl_vector *stateInclusion2,int begin, int length,
			    int transitionPowerFlag){
  double score1, score2;
  score1 = hmm_subpath_score(hmm, emisLogs, stateInclusion1, begin,length,transitionPowerFlag);
  score2 =  hmm_subpath_score(hmm, emisLogs, stateInclusion2, begin,length,transitionPowerFlag);
  //printf("scores: %f %f\n",score1,score2);
  return(score1 - score2);
}


/* simulatePath-- returns a vector snpNumber long that represents an instance
of the simulated markov state path through the data. takes as arguments,an hmm,
 and a pointer to void which actually is an gsl random numb. gen, as well as the observation
 number */
gsl_vector *simulatePath(HMM *h, void *r, int obsNumber){
  int i, j, pastState;
  double sum, rand;
  gsl_rng * rn = (gsl_rng *) r;
  gsl_vector *path;

  path = gsl_vector_alloc(obsNumber);

  /*choose initial state, set paths[0] */
  rand = gsl_rng_uniform(rn);
  sum = 0;
  for(i = 0; i < h->nstates; i++){
    sum += gsl_vector_get(h->piStart, i);
    printf("%f %f\n",rand, sum);
    if (rand <= sum){
      printf("true\n");
      gsl_vector_set(path,0,i);
      rand = 2.0;
    }
  }
   
  /* now go through obs, using transition matrix to choose states */
  
  for(i = 1; i < obsNumber; i++){
    pastState = gsl_vector_get(path, i-1);
    rand = gsl_rng_uniform(rn);
    sum = 0;
    for(j = 0; j < h->nstates; j++){
      sum += gsl_matrix_get(h->transition_matrix, pastState, j);
      if (rand <= sum){
	gsl_vector_set(path,i,j);
	rand = 2.0;
      }
    }
  }
  return(path);
}


/* simulatePathSpaced-- just like simulatePath except takes two
added parameters, a pointer to a gsl_vector which will contain the 
simulated snp locations, and a double for the exponential mean */
gsl_vector *simulatePathSpaced(HMM *h, void *r, int obsNumber, gsl_vector *locs, double mu){
  int i, j, pastState;
  double sum, rand;
  gsl_rng * rn = (gsl_rng *) r;
  gsl_vector *path;
  gsl_matrix *matPower;

  path = gsl_vector_alloc(obsNumber);

  /*choose initial state, set paths[0] */
  rand = gsl_rng_uniform(rn);
  sum = 0;
  for(i = 0; i < h->nstates; i++){
    sum += gsl_vector_get(h->piStart, i);
    printf("%f %f\n",rand, sum);
    if (rand <= sum){
      printf("true\n");
      gsl_vector_set(path,0,i);
      rand = 2.0;
    }
  }
 
  /*set initial location to 1 by convention */
  gsl_vector_set(locs,0,1);
  //now go through obs and simulate spacing -- curently drawn from exp
  for(i = 1; i < obsNumber; i++){
    rand = gsl_ran_exponential(rn, mu);
    gsl_vector_set(locs, i, gsl_vector_get(locs, i - 1) + ceil(rand));
  }

  /* now go through obs, using transition matrix and powers to choose states */
  for(i = 1; i < obsNumber; i++){
    pastState = gsl_vector_get(path, i-1);
    matPower = gsl_matrix_power(h->transition_matrix, gsl_vector_get(locs,i) - gsl_vector_get(locs,i - 1));
    rand = gsl_rng_uniform(rn);
    sum = 0;
    for(j = 0; j < h->nstates; j++){
      sum += gsl_matrix_get(matPower, pastState, j);
      if (rand <= sum){
	gsl_vector_set(path,i,j);
	rand = 2.0;
      }
    }
    gsl_matrix_free(matPower);
  }
  return(path);
}
