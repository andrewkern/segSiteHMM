
/* segSiteHMM v. 0.1 !!!
/
/
/ Andrew Kern 3/14/07
*/


#include "../hmm/hmm.h"
#include "../hmm/adkGSL.h"
#include "../hmm/popGenTools.h"
#include "segSites.h"
#include "segSiteHMM.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>


#define MAXSNPS 100000000

char modeFlag, *outfileName, smoothFlag;
int snpNumber, maxSampleSize, states, verboseFlag = 0, fixedEmisFlag = 0, sampleSize;
gsl_vector *sampleSizes, *piStart, *piEnd;
struct snp data[MAXSNPS];
gsl_matrix *m, *sfsCounts, *emisMatrix;
double gamma1, theta1;
gsl_vector *gammas;
double priorCons = 0;
double omega = 0;
int outputLevel = 0;

int main(int argc, char *argv[]){
  gsl_matrix *alphas, *betas, *posts, *logPosts;
	gsl_vector  *vits, *locs, *include, *include2;
	gsl_matrix  *emisLogs;
	double score1, sum, pFor;
	HMM *h;
	FILE *outfile;
	int i, j, **elements, elementCount, transitionPowerFlag;

	

	printf("segSiteHMM v0.3 -- thank you for playing thermonuclear war\n");
	printf("getting parameters/data\n");
	getParameters(argc,argv);




	//find max sampleSize
	maxSampleSize = 0;
	for(j = 0; j < snpNumber; j++){
	  if (data[j].n > maxSampleSize){
	    maxSampleSize = data[j].n;
	  }
	}
	//transition gaps?
	transitionPowerFlag = 0;

	//make sampleSize vector and summarize sfs
	sampleSizes = sampleSizeVector(data,snpNumber,maxSampleSize);
	sfsCounts = summarizeSFSBool(maxSampleSize, data, snpNumber);
	

	
	//was smoothing used?
	if(smoothFlag == 'f'){
	  //init with equal starting probs
	  sum = 1.0 / states;
	  for (i = 0; i < states; i++){
	    gsl_vector_set(piStart, i, sum);
	    gsl_vector_set(piEnd, i, sum);
	  }
	}
	else if(smoothFlag == 's'){
	  //set start and end probs based on smoothing
	  gsl_vector_set(piStart, 1, priorCons);
	  gsl_vector_set(piStart, 0, 1.0 - gsl_vector_get(piStart, 1));
	  gsl_vector_set(piEnd, 0, gsl_vector_get(piStart,0));
	  gsl_vector_set(piEnd, 1, gsl_vector_get(piStart,1));
	  
	  //define transition probs based on omega and gamma
	  gsl_matrix_set(m, 1, 0, 1.0 / omega);
	  gsl_matrix_set(m, 1, 1, 1.0 - gsl_matrix_get(m, 1, 0));
	  gsl_matrix_set(m, 0, 1, (gsl_matrix_get(m, 1, 0) / gsl_vector_get(piStart, 0)) \
			 -  gsl_matrix_get(m, 1,0));
	  gsl_matrix_set(m, 0, 0, 1.0 - gsl_matrix_get(m, 0, 1));
	}
	  
	  
	printf("path is %d SNPs long\n",snpNumber);
	locs = gsl_vector_alloc(snpNumber);
	//fill up locs vector for use in calculating transition powers
	for(i = 0; i < snpNumber; i++){
	  gsl_vector_set(locs, i, data[i].pos);
	}



	//initialize HMM
	printf("initializing model\n");
	h = newHMM(m, piStart,piEnd,locs);
	initHMM(h);
	printHMM(h);
	
	switch(modeFlag){
	case 't' :
	  printf("training mode\nEM iterations follow:\n");
	  fflush(stdout);
	  emAlg(h, gammas,snpNumber, data, transitionPowerFlag);
	  //output selection coefficient
	  strcat(outfileName,".alphas");
	  outfile = fopen(outfileName, "w");
	  gsl_vector_fprintf(outfile,gammas,"%f");
	  fclose(outfile);
	  //output transition matrix
	  strcat(outfileName,".trainedTransMat");
	  printTransitions(h, outfileName);
	  break;
	case 'b' :
	  printf("train and decode mode\nEM iterations follow:\n");
	  fflush(stdout);
	  emAlg(h, gammas,snpNumber, data, transitionPowerFlag);
	  	if(outputLevel>0){
		//output selection coefficient
		  strcat(outfileName,".alphas");
		  outfile = fopen(outfileName, "w");
		  fprintf(outfile,"%f",gsl_vector_get(gammas,1));
		  fclose(outfile);
		  //output transition matrix
		  strcat(outfileName,".trainedTransMat");
		  printTransitions(h, outfileName);
		}
	  strcat(outfileName,".decoding");
	  printf("Viterbi decoding mode\nelement predicition will be stored in a file called \"%s\"\nPosterior Probs and Viterbi calls follow:\n",outfileName);
	  vits = gsl_vector_alloc(snpNumber);
          if(!fixedEmisFlag)
          {
	    emisLogs = gsl_matrix_alloc(h->nstates, snpNumber);
	    calculateEmissions(h, emisLogs, gammas, theta1, snpNumber, data);
          }
	  alphas = gsl_matrix_alloc(h->nstates,snpNumber);
	  betas = gsl_matrix_alloc(h->nstates,snpNumber);  
	  posts = gsl_matrix_alloc(h->nstates, snpNumber);
	  logPosts =  gsl_matrix_alloc(h->nstates, snpNumber);
	  pFor = forwardAlg(h, alphas, emisLogs, snpNumber, transitionPowerFlag);
	  backwardAlg(h, betas, emisLogs, snpNumber, transitionPowerFlag);
	  posteriorProbsReduced(h, posts, pFor, alphas, betas, emisLogs, snpNumber, transitionPowerFlag);
	  viterbi(h, vits, emisLogs, snpNumber, transitionPowerFlag);

	  //print out likelihood
	  printf("Likelihood = %lf\n",pFor);
	
	if (outputLevel > 0){	
		//output posteriors and viterbi calls
		outfile = fopen(outfileName, "w");
		for(i = 0; i < emisLogs->size2; i++){
			fprintf(outfile,"%d",(int)gsl_vector_get(h->otherData,i));
			for(j = 0; j < h->nstates; j++){
				fprintf(outfile,"\t%lf",gsl_matrix_get(posts, j, i));
			}
			fprintf(outfile,"\t%d\n",(int)gsl_vector_get(vits, i));
		}
		fclose(outfile);
  	}
	  //score elements
	  strcat(outfileName,".elements");
	  scoreElements(h, emisLogs, vits, outfileName, transitionPowerFlag);
	  break;
	case 'v' :
	  printf("Viterbi decoding mode\nelement predicition will be stored in a file called \"%s\"\nPosterior Probs and Viterbi calls follow:\n",outfileName);
	  vits = gsl_vector_alloc(snpNumber);
	  emisLogs = gsl_matrix_alloc(h->nstates, snpNumber);
          if (fixedEmisFlag)
          {
              initializeEmissions(emisLogs, emisMatrix, snpNumber, data);
          }
          else
          {
	      calculateEmissions(h, emisLogs, gammas, theta1, snpNumber, data);
          }
	  alphas = gsl_matrix_alloc(h->nstates,snpNumber);
	  betas = gsl_matrix_alloc(h->nstates,snpNumber);  
	  posts = gsl_matrix_alloc(h->nstates, snpNumber);
	  logPosts =  gsl_matrix_alloc(h->nstates, snpNumber);

	  pFor = forwardAlg(h, alphas, emisLogs, snpNumber, transitionPowerFlag);
	  backwardAlg(h, betas, emisLogs, snpNumber, transitionPowerFlag);
	  posteriorProbsReduced(h, posts, pFor, alphas, betas, emisLogs, snpNumber, transitionPowerFlag);
	  viterbi(h, vits, emisLogs, snpNumber, transitionPowerFlag);
	  
	  //print out likelihood
	  printf("Likelihood = %lf\n",pFor);

	if(outputLevel > 0){
	//output posteriors and viterbi calls
		outfile = fopen("snpPosts", "w");
		for(i = 0; i < emisLogs->size2; i++){
			fprintf(outfile,"%d",(int)gsl_vector_get(h->otherData,i));
			for(j = 0; j < h->nstates; j++){
				fprintf(outfile,"\t%lf",gsl_matrix_get(posts, j, i));
			}
			fprintf(outfile,"\t%d\n",(int)gsl_vector_get(vits, i));
		}
		fclose(outfile);
		scoreElements(h, emisLogs, vits, outfileName, transitionPowerFlag);
		break;
		}
	}
	
	freeHMM(h);
	return 0;
}



void usage(){
    fprintf(stderr,"usage: pgHMM - population genetic HMMs\nrun modes:\n\t-t train mode -- EM training for model\n\t-v viterbi decoding mode\n\t-b both train and decode\n\ninput control\n\t-S stateNumber\n\t-d dataFile\n\t-p transitionMatrix\n\t-s startProbs\n\t-e endProbs\n\t-g prior probability of being in selected state\n\t-E emissionMatrixFile sampleSize -- specifies a file containing an SFS for each state (including frequencies of monomorphic sites); entires must sum to 1, and data must have constant sampleSize\n\t-w omega smoothing parameter (note either transitionMatrix or gamma and omega need to be defined)\n\t-a beta1 beta2 -- value for selected classes\n\t-u 4Nu initial value\n\t-o outfile name (trained matrix or elements)\n\t-V verbose\n\t-O output level (default= 0)\n");
    exit(1);
}

/* sets the parameters/run mode and then opens data */
void getParameters(int argc, char *argv[]){
  FILE *infile;
  int  i, j, n, args;
  long int pos;
         
  if (argc < 2){
    usage();
  }
  else{
    args = 1;
    while(args < argc){
      switch(argv[args][1]){
      case 't' :
	/* train mode? */
	modeFlag = 't';
	smoothFlag = 'f';
	break;
      case 'v' :
	/* viterbi decode */
	modeFlag = 'v';
	smoothFlag = 'f';
	break;
      case 'b' :
	/*both train and decode */
	modeFlag = 'b';
	smoothFlag = 'f';
	break;
      case 'a' :
	/*set alphas for selected class */
	for(j = 1; j < states; j++){
	  gsl_vector_set(gammas,j,atof(argv[++args]));
	}
	break;
      case 'u' :
	//set theta 
	theta1 = atof(argv[++args]);
	break;
	case 'O' :
	//set outputLevel 
	outputLevel = atoi(argv[++args]);
	break;
      case 'S' :
	//set state number
	states = atoi(argv[++args]);
	m = gsl_matrix_alloc(states,states);
	gsl_matrix_set_zero(m);
	//initialize selection params and allocate stuff
	gammas = gsl_vector_alloc(states);
	piStart = gsl_vector_alloc(states);
	piEnd = gsl_vector_alloc(states);
	break;
      case 'd' :
      	/* open data file, errors? */
	infile = fopen(argv[++args], "r");
	printf("%s is datafile\n",argv[args]);
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get snp info; expect pos, i, n */
	snpNumber = 0;
	maxSampleSize = 0;
	while (fscanf(infile, "%ld\t%d\t%d", &pos, &i, &n) != EOF){
		data[snpNumber].pos = pos;
		data[snpNumber].i = i;
		data[snpNumber].n = n;
		snpNumber++;
		if(n > maxSampleSize){
		  maxSampleSize = n;
		}
	}
	fclose(infile); 
	break;
	//either whole transition matrix needs to be defined or gamma and omega
      case 'p':
	//get transition matrix
	infile = fopen(argv[++args], "r");
	if (infile == NULL){
	  fprintf(stderr,"Error opening matrixfile! ARRRRR!!!!\n");
	  exit(1);
	}
	printf("%s is matrixfile\n",argv[args]);
	gsl_matrix_fscanf(infile, m);
	//gsl_matrix_fprintf(stdout,m,"%f");
	fclose(infile);
	break;
      case's':
	//get start vector
		infile = fopen(argv[++args], "r");
	if (infile == NULL){
	  fprintf(stderr,"Error opening start prob. file! ARRRRR!!!!\n");
	  exit(1);
	}
	printf("%s is start prob. file\n",argv[args]);
	gsl_vector_fscanf(infile, piStart);
	fclose(infile);
	smoothFlag = 'piSet';
	break;
      case'e':
	//get end vector
		infile = fopen(argv[++args], "r");
	if (infile == NULL){
	  fprintf(stderr,"Error opening end prob. file! ARRRRR!!!!\n");
	  exit(1);
	}
	printf("%s is end prob. file\n",argv[args]);
	gsl_vector_fscanf(infile, piEnd);
	smoothFlag = 'piSet';
	fclose(infile);
	break;
      case 'E':
        //get emission matrix
        infile = fopen(argv[++args],"r");
        sampleSize = atoi(argv[++args]);
        if (infile == NULL){
          fprintf(stderr,"Error opening end emission prob. file! ARRRRR!!!!\n");
          exit(1);
        }
        printf("%s is emission prob. file; sampleSize=%d\n",argv[args-1],sampleSize);
	emisMatrix = gsl_matrix_alloc(states,sampleSize);
        gsl_matrix_fscanf(infile, emisMatrix);
        printf("finished reading emission prob. matrix\n");
        for (i=0;i<states;i++)
        {
            for (j=0;j<sampleSize;j++)
            {
                printf("%f ",gsl_matrix_get(emisMatrix,i,j));
            }
            printf("\n");
        }
        fixedEmisFlag = 1;
        fclose(infile);
        break;

      case 'g':
	//get gamma for transitions
	priorCons = atof(argv[++args]);
	printf("prior probability of conservation defined\n");
	break;
      case 'w':
	//get omega for transitions
	omega =  atof(argv[++args]);
	smoothFlag = 's';
	break;
      case 'o':
	//get outfile name
	outfileName = argv[++args];
	break;
      case 'V':
	//verbose training
	verboseFlag = 1;
	break;
      }
      args++;
    }
  }
  if (omega == 0 && gsl_matrix_get(m, 0, 0) == 0){
    usage();
    exit(1);
  }
}

/* calculateEmissions- this fills up a matrix of Nstates * L observations where the elements are
   the computed allele frequency probabilities drawn from the stationary distribution of Wright's genic
   selection model. The population is assumed to be at equilibrium and all snps are considered to be
    independent draws from the same distributions */
void calculateEmissions(HMM *hmm, gsl_matrix *emisLogs, gsl_vector *gammas, double theta,\
			int snpNumber, struct snp theData[MAXSNPS]){
  int i, j;
  gsl_matrix *probs;
  gsl_vector *pVector;
  
  //error check
  assert(emisLogs->size1 == hmm->nstates);
  assert(emisLogs->size2 == snpNumber);
  //alloc pVector
  pVector = gsl_vector_alloc(hmm->nstates);
  probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize); 
  //calculate emissions and fill up matrix
  for(i = 0; i < emisLogs->size1; i++){
    gsl_vector_set(pVector,0,gsl_vector_get(gammas,i));
    gsl_vector_set(pVector,1,theta);
    //make prob matrix (note this is in logs)
    siteProbMatrix(probs,pVector,maxSampleSize, sampleSizes);
    //go through snps, calculate probs
    for(j=0; j < emisLogs->size2; j++){
      gsl_matrix_set(emisLogs, i, j, gsl_matrix_get(probs, data[j].n, data[j].i));
    }
  }
  gsl_vector_free(pVector);
  gsl_matrix_free(probs);
}

void initializeEmissions(gsl_matrix *emisLogs, gsl_matrix *emisMatrix, int snpNumber, struct snp theData[MAXSNPS]){
  int i, j;
  int freq;
  double logProb;

  //calculate emissions and fill up matrix
  for(i = 0; i < emisLogs->size1; i++){
    for(j=0; j < emisLogs->size2; j++){
      //need to get the allele frequency here, and then extract the probs from the corresponding emisMatrix elements
      freq = data[j].i;
      logProb = log(gsl_matrix_get(emisMatrix, i, freq));
      //if (emisLogs->size2 - j <= 10)
      //{
      //    fprintf(stderr,"emission probability for state %d, position %d: %f; log=%f\n",i,j,gsl_matrix_get(emisMatrix, i, freq),logProb);
      //}
      gsl_matrix_set(emisLogs, i, j, logProb);
    }
  }
}

/* emAlg- this is mean to train the emissions (the parameter gamma) using the EM algorithm. 
 it returns the trained up gamma value  */
double emAlg(HMM *hmm, gsl_vector *startGammas, int snpNumber, struct snp theData[MAXSNPS], int transitionPowerFlag){
  gsl_matrix *alphas, *betas, *posts, *logPosts, *emisLogs, *tempTrans, *powTrans, **jpMatArray;
  gsl_vector *rowSums, *estimateVector;
  double  lastLik, newLik, delta, dummy, val, pFor, pBack;
  gsl_function likFunc;
   struct my_jointLik_params paramSet;
  int i, j,  k, count;
	
  //initialize and allocate
  emisLogs = gsl_matrix_alloc(hmm->nstates,snpNumber);
  if(fixedEmisFlag)
  {
      initializeEmissions(emisLogs, emisMatrix, snpNumber, theData);
  }
  else
  {
      calculateEmissions(hmm, emisLogs, startGammas, theta1, snpNumber, theData);
  }
  posts = gsl_matrix_alloc(hmm->nstates,snpNumber);
  logPosts =  gsl_matrix_alloc(hmm->nstates,snpNumber);
  alphas = gsl_matrix_alloc(hmm->nstates,snpNumber);
  betas = gsl_matrix_alloc(hmm->nstates,snpNumber);  
  rowSums = gsl_vector_alloc(hmm->nstates);
  delta = 1e-4;
  newLik = 1 ;
  lastLik = 0;
  tempTrans = gsl_matrix_alloc(hmm->nstates,hmm->nstates);
  
  //make jpMatArray
  jpMatArray = malloc(sizeof(gsl_matrix) * snpNumber);
  for(i = 0; i < snpNumber; i++){
    //each is nstate x nstate transition matrix
    jpMatArray[i] = gsl_matrix_alloc(hmm->nstates,hmm->nstates);
    gsl_matrix_set_all(jpMatArray[i],0.0);
  }

				     //setup parameter estimate vector
  estimateVector = gsl_vector_alloc(hmm->nstates); // 1 mutation parameter, n-1 selection parameters
  gsl_vector_set(estimateVector,0,theta1);
  for(i = 1; i < hmm->nstates; i++){
    gsl_vector_set(estimateVector,i,gsl_vector_get(startGammas,i));
  }

  //initializing joint lik params  
  paramSet.data = theData;
  paramSet.snpNumber = snpNumber;
  paramSet.nstates = hmm->nstates;
  paramSet.maxSampleSize = maxSampleSize;
  paramSet.sampleSizeVector = sampleSizes;
  paramSet.sfsBools = sfsCounts;
  paramSet.startingPoint = estimateVector;

  //main loop
  count = 0;
  while ( fabs(newLik - lastLik) >  delta){
    transitionPowerFlag = 0;
    
    //E step
    //set lastLik to newLik
    dummy = newLik;
    lastLik = dummy;
		
    //recalculate emissions
    if (!fixedEmisFlag)
    {
        calculateEmissions(hmm, emisLogs, startGammas, theta1, snpNumber, theData);
    }
		
    //run forwards/backwards and calculate posterior probs
    pFor = forwardAlg(hmm, alphas, emisLogs, snpNumber, transitionPowerFlag);
    pBack = backwardAlg(hmm, betas, emisLogs, snpNumber,transitionPowerFlag);
    //print current status
    printf("%d\tlik= %f\t\ttheta_hat= %f", count,pFor,theta1);
    for(i = 1; i < hmm->nstates; i++){
      printf("\tgamma_hat= %f",gsl_vector_get(startGammas, i));
    }
    printf("\n");
    fflush(stdout);
    //check probs; calculate posterior probs; set "weights"
    assert(fabs(pFor - pBack) < 0.1);
    posteriorProbsReduced(hmm, posts, pFor, alphas, betas, emisLogs, snpNumber,transitionPowerFlag);
    posteriorProbsReducedLog(hmm, logPosts, pFor, alphas, betas, emisLogs, snpNumber,transitionPowerFlag);

    //reestimate the new shiny way
    reestimateTransitions(hmm, posts, jpMatArray,  pFor, alphas, betas, emisLogs, snpNumber,transitionPowerFlag);

    //update start and end vectors based on posteriors
    for(j = 0; j < hmm->nstates; j++){
      gsl_vector_set(hmm->piStart,j,gsl_matrix_get(posts,j,0));
      gsl_vector_set(hmm->piEnd,j,gsl_matrix_get(posts,j,snpNumber - 1));
      gsl_vector_set(hmm->piStart_scores, j, log(gsl_vector_get(hmm->piStart,j)));
      gsl_vector_set(hmm->piEnd_scores, j, log(gsl_vector_get(hmm->piEnd,j)));
    } 
          
    // M step for other params
    if (!fixedEmisFlag)
    {
        paramSet.posts = posts;
        paramSet.startingPoint = estimateVector;
        likFunc.params = &paramSet;
        jointMLEstNState(estimateVector, &newLik, &paramSet);
        
        //maximization
        theta1 = gsl_vector_get(estimateVector,0);
        for(j = 1; j < hmm->nstates; j++){
          gsl_vector_set(startGammas,j, gsl_vector_get(estimateVector, j));
        }
    }
    newLik = pFor;
    count += 1;
    if (verboseFlag) printHMM(hmm);
  }

 //free jpMatArray
  for(i = 0; i < snpNumber; i++){
    //each is nstate x nstate transition matrix
    gsl_matrix_free(jpMatArray[i]);
  }
  //free(jpMatArray);

  gsl_matrix_free(emisLogs);
  gsl_matrix_free(posts);
  gsl_matrix_free(alphas);
  gsl_matrix_free(betas);
  gsl_vector_free(rowSums);
  gsl_vector_free(estimateVector);
  printHMM(hmm);
  return(gsl_vector_get(startGammas,1));
}	

/*This returns an array of arrays of ints which indicate the Viterbi decoded elements from the model.
  The arrays returned are in the format startIndex, stopIndex */
int **elementArray(gsl_vector *vits, double state, int *pCount){
  int i, **array, elementCount, size, last;
 
  //allocate storage for element arrays
  array = malloc(vits->size * sizeof(int *));
  for(i = 0; i < vits->size; i++){
    array[i] = malloc(3 * sizeof(int));
  }
  //set last to initial obs
  elementCount = 0;
  size = state; //this is the first one we're about to find!
  //find first occurence of state
  for(last = 0; last < vits->size; last++){
    if (gsl_vector_get(vits,last) == state){
      break;
    }
  }

  //iterate through Viterbi's, find runs
  for(i = last+1; i < vits->size; i++){
    if (gsl_vector_get(vits,last) == gsl_vector_get(vits,i)){
      if (size == 0){ //new one?
	last = i;
      }
      size++;
      if ( i == vits->size - 1 && size > 1){
	array[elementCount][0] = last;
	array[elementCount][1] = i;
	elementCount++;
      }
    }
    else{
      if (size > 1){
	array[elementCount][0] = last;
	array[elementCount][1] = i - 1;
	array[elementCount][2] = size;
	elementCount++;
	size = 0;
	last = i - 1 ;
      }
      else{
	size = 0;
      }
    }
  }
  *pCount = elementCount;
  return(array);
}

void scoreElements(HMM *h, gsl_matrix *emisLogs, gsl_matrix *vits, char *outfileName, int transitionPowerFlag){
  gsl_vector *include, *include2;
  FILE *outfile;
  int elementCount, **elements, i, j;
  double score1;

  //open file to print to
  outfile = fopen(outfileName, "w");
  if (outfile == NULL){
    fprintf(stderr,"Error opening outfile! ARRRRR!!!!\n");
    exit(1);
  }
  //print header
  fprintf(outfile,"state\tstart\tstop\tlog_odds_score\tlength\n");

  //alloc state inclusion vectors
  include = gsl_vector_alloc(h->nstates);
  include2 = gsl_vector_alloc(h->nstates);

  //iterate across states
  for(j = 0; j < h->nstates; j++){
    //get state elements from Viterbi
    elementCount = 0;
    elements = elementArray(vits, j, &elementCount);
    //set up state sets for element scoring
    gsl_vector_set_all(include,0);
    gsl_vector_set(include, j, 1);
    gsl_vector_set_all(include2,1);
    gsl_vector_set(include2, j, 0);
    //score elements
    for(i=0; i < elementCount; i++){
      score1 = hmm_subpath_log_odds(h, emisLogs, include,include2, elements[i][0], \
				    elements[i][1] - elements[i][0] + 1 , transitionPowerFlag);
      fprintf(outfile,"%d\t%d\t%d\t%f\t%d\n",j,(int) gsl_vector_get(h->otherData,elements[i][0]), \
	      (int) gsl_vector_get(h->otherData,elements[i][1]),score1, (int) elements[i][2]);
      free(elements[i]);
    }
    free(elements);
  }
  fclose(outfile);
}

