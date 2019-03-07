/* hmm class header */
// Andrew Kern 8/24/05
//

#include "adkGSL.h"
#include <math.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

typedef struct {
  int nstates;
  gsl_matrix *transition_matrix;
  gsl_matrix *transition_score_matrix; /* entries are logs of entries
                                        * in transition matrix */
  gsl_vector *piStart, *piEnd, *piStart_scores, *piEnd_scores;
  gsl_vector *otherData; //this is meant to hold distances between obs
} HMM;

HMM* newHMM(gsl_matrix *mat, gsl_vector *piStart, gsl_vector *piEnd, gsl_vector *otherData);
void freeHMM(HMM *hmm);
void printHMM(HMM *hmm);
void initHMM(HMM *hmm);
double forwardSumPrevious(HMM *hmm, gsl_matrix *alphas, gsl_matrix *transition_score_matrix_i, int state, int obsIndex);
double viterbiMaxPrevious(HMM *hmm, gsl_matrix *deltas, gsl_matrix *transition_score_matrix_i, int state, int obsIndex, int *point);
double backwardSumNext(HMM *hmm, gsl_matrix *alphas,  gsl_matrix *transition_score_matrix_i, gsl_matrix *emissions, int state, int obsIndex);
double forwardAlg(HMM *hmm, gsl_matrix *alphas, gsl_matrix *emissions, int length, int transitionFlag);
double backwardAlg(HMM *hmm, gsl_matrix *betas, gsl_matrix *emissions, int length,int transitionPowerFlag);
void posteriorProbs(HMM *hmm, gsl_matrix *postProbMat, gsl_matrix *emissions, int length, int transitionPowerFlag);
void posteriorProbsReduced(HMM *hmm, gsl_matrix *postProbMat, double pFor, gsl_matrix *alphas, gsl_matrix *betas, gsl_matrix *emissions, int length, int transitionPowerFlag);
void reestimateTransitions(HMM *hmm, gsl_matrix *postProbMat, gsl_matrix **jpMatArray, double pFor, gsl_matrix *alphas, \
			   gsl_matrix *betas, gsl_matrix *emisLogs, int length, int transitionPowerFlag);
void reestimateTransitions_vanilla(HMM *hmm, double pFor, gsl_matrix *alphas, \
	gsl_matrix *betas, gsl_matrix *emisLogs, int length);
		
double viterbi(HMM *hmm, gsl_vector *viterbiPath, gsl_matrix *emissions,  int length, int transitionPowerFlag);
void hmm_logify_transitions(HMM *hmm);
void hmm_normalize_transitions(HMM *hmm);
double forwardAlgOffset(HMM *hmm, gsl_matrix *alphas, gsl_matrix *emissions, int length, int offset, int transitionPowerFlag);
double hmm_subpath_score(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *stateInclusion, int begin, int length,  int transitionPowerFlag);
double hmm_subpath_log_odds(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *stateInclusion1,gsl_vector *stateInclusion2,int begin, int length, int powerTransitionFlag);
gsl_vector *simulatePath(HMM *h, void *r, int obsNumber);
gsl_vector *simulatePathSpaced(HMM *h, void *r, int obsNumber, gsl_vector *locs, double mean);
void posteriorDecode(HMM *hmm, gsl_vector *posteriorPath, gsl_matrix *posts, int length);
void printTransitions(HMM *h, char *outfileName);