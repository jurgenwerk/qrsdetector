#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wfdb/wfdb.h>
#include <wfdb/ecgmap.h>

//Author: Matic Jurglic

#define MAX_LEADS 3
#define MAX_SIGS  5
#define MAX_LEN   31000000
#define MAX_BUF   100000
#define MAX_GAP   1000     

static WFDB_Siginfo s[MAX_LEADS];

unsigned stklen = 40000;
long nrSamps;
double *sig0;
double *sig1;
double *sig2;

int M;
int N;
int openLeads;
long samples;
double normCnst = 32;

int
openRecord(char *record){
  return (isigopen(record, s, MAX_LEADS));
}

long ReadBuffer(long annot){
  long i;
  long j;
  WFDB_Sample vec[MAX_LEADS];

  for (i=0;i<nrSamps;i++){
    if (getvec(vec)<0){
      fprintf(stderr,"Error\n");
      break;
    }
    sig0[i]=vec[0];
  }
  return(nrSamps);
}

void 
writeQRS(char *record, int chan){
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  int i;
  FILE *out;

  annIFO.name = "qrs"; annIFO.stat = WFDB_WRITE;
  if (annopen(record, &annIFO, 1) < 0){
    fprintf(stderr,"Error opening QRS file\n");
    return;
  }
  annot.subtyp = annot.chan = annot.num = 0; annot.aux = NULL;
  for (i=0;i<samples;i++){
    if (sig1[i]!=0){
      annot.anntyp = sig1[i];
      annot.time=i;
      if (putann(0, &annot) < 0) break;
    }
  }
}

void 
readQRS(char *record, char *ant, int chan){
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  long i;
  FILE *out;

  annIFO.name = ant; annIFO.stat = WFDB_READ;
  if (annopen(record, &annIFO, 1) < 0){
    fprintf(stderr,"Error opening QRS file\n");
    return;
  }
  annot.subtyp = annot.chan = annot.num = 0; annot.aux = NULL;
  for (i=0;i<samples;i++) sig1[i]=0;
  i=0;
  while (getann(0, &annot)==0) {
    if (annot.time>samples) {
      fprintf(stderr,"Error reading annotation times\n");
      return;
    }
    sig1[annot.time]=annot.anntyp;
  }
}
double movingAverage(int sample){
  long i;  
  double temp = 0;
  for(i=0; i<M; i++)
		temp += sig0[sample-i];
  temp /= M;
  return temp;
}
double summation(int sample){
  long i;
  double temp = 0;
  
	for (i=N; i>0 ;i--) 
		temp += sig1[sample - i];    
  return temp;
  
}
void findPeaks(){
  
  int startedClimbing = 0;
  long i,j;
  long start=0;
  long end=0;
  long peakIdx;
  double peak=0;
  double temp; 
  double descentFactor = 1;
  double climbingFactor = 1.6;
  
  
  for (i=0;i<samples;i++) 
		sig1[i] = 0;

	for(i=N; i<samples; i++){
		if(startedClimbing == 0){
			if((sig2[i]/sig2[i-1]) > climbingFactor){ //starts climbing
				start = i;
				startedClimbing = 1;
			}				
		}
		else if((sig2[i]/sig2[i-1]) < descentFactor){  //starts descending
			startedClimbing = 0;			
			end = i;		
			for(j = start; j<end; j++){ //finding peak between climbing and descending index
				temp = fabs(sig0[j]);				
				if(temp > peak){
					peak = temp;
					peakIdx = j;
				}				
			}
			sig1[peakIdx] = peak;
			peak = 0;
		}
	}
}
double decisionMaking(double alpha, double gamma){
  
  long i;
  int j;
  double threshold = 0;
  double temp = 0;

	for(i=0; i<nrSamps; i++){
		if(sig1[i]>0){
			if(threshold == 0){ //init threshold
				for (j = 0; j < 200; j++)
					temp+=sig1[j];
				threshold = temp/200.0;
			  sig1[i] = 1;
			}
			else if (sig1[i]>threshold){
				threshold = alpha*gamma*sig1[i] + (1 - alpha)*threshold;
				sig1[i]=1;
			}
			else sig1[i]=0;
		}
	}
  return 0;
}
double detectQRS(){
  long length;
  long i,j,k, max;  
	
  M=9;
  N=38;

  double alpha=0.06;
  double gamma=0.18;
  length = samples;

  //step 1: Linear High pass filtering stage
  
	//step 1.1: M-point moving average 
  for (i=M;i<length;i++) 
		sig1[i] = movingAverage(i);
  
  //step 1.2: delay system with a group delay of (M+1)/2 samples and subtracting y2 - y1 
  for (i=(M+1)/2;i<length;i++) 
		sig1[i] = sig0[i-((M+1)/2)]-sig1[i];

  //step 2: Nonlinear Low pass filtering stage
  
	//step 2.1: point by point squaring
  for (i=0;i<length;i++) 
  	sig1[i] = sig1[i]*sig1[i];

  //step 2.2: moving summation window of size N
  for (i=N;i<length;i++) 
		sig2[i] = summation(i);
   
	//step 3: find peaks and make a decision 
  findPeaks();
  decisionMaking(alpha,gamma);

}

int
main(int argc, char *argv[]){
  long i, j;
  char *record = NULL;
  char *annotator=NULL;
  double thresh;	
  
  for (i=1;i<argc;i++){
    if (argv[i][0]!='-'){
      fprintf(stderr,"Error parsing command line\n");
      exit(1);
    }
    switch (argv[i][1]){
    case 'r':
      i++;
      record = (char *)malloc(sizeof(char)*(strlen(argv[i])+2));
      strcpy(record,argv[i]);
      break;
    case 'a': 
      i++;
      annotator = (char *) malloc(sizeof(char)*(strlen(argv[i])+2));
      strcpy(annotator, argv[i]);
      break;
    case 'n': 
      i++;
      normCnst = atof(argv[i]);
      break;
    default:
      fprintf(stderr,"Wrong switch\n");
      exit (2);
    }
  }
  if (record == NULL){
    fprintf(stderr,"No record was specified, exiting\n");
    exit (2);
  }
  if ((openLeads=openRecord(record))<0){
    fprintf(stderr,"Error opening record, exiting\n");
    exit (3);
  }
  nrSamps = s->nsamp;	
  sig0 = (double *)malloc(sizeof(double ) * nrSamps);
  sig1 = (double *)malloc(sizeof(double ) * nrSamps);
  sig2 = (double *)malloc(sizeof(double ) * nrSamps);
  
  if ((samples=ReadBuffer(0))<=0){
    fprintf(stderr,"Error opening record, exiting\n");
    exit (3);
  }

  detectQRS();
  writeQRS(record, 1);
  wfdbquit();
  free (sig0);
  free (sig1);
  free (sig2);
  return 0;
}




