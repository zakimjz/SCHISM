#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<fstream>
#include<algorithm>
#include<ext/algorithm> // for linux code
//#include<algorithm>     // for solaris/unix code
#include<utility>
#include<unistd.h>      // for getopt() function
using namespace std;
using namespace __gnu_cxx; // for linux code

enum Digit {MINUS_ONE=-1,ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,TEN};
enum DistMode {FULL_DIM=0,HYPER_RECT_UNIFORM,MVNORMAL,LINEAR_DEPENDENCE};
static const int RANGE = 1000;
static const int SAMPLE_SIZE = 1000000;
static const int NUM_ARGS = TWO;
//input args
int iNoOfClusters=5,iNoOfPts=1000,iNoOfDims=400,iRealDims=50,iMode=2,iNoOfDatasets=1;
double dProbNoise=0.05,dStdDev=0.5,dMaxClusterSzRatio=4.0,dProbOverlap = 0.5,dProbConstrain=0.5,x = 0.5,dBias = 0.05,dWidth=20;//0.25*RANGE;
char caOutFile[300];
bool bCentresOnly = false;

// see Probability Methods for Computer Science - Sheldon Ross (Chpt 9)
double genStdNormal() {   
  while(1) {
    double dY1 = -log(1.0 - (rand()%RANGE)/(double)RANGE);
    double dY2 = -log(1.0 - (rand()%RANGE)/(double)RANGE);
    double dY = (dY2-(dY1-1)*(dY1-1))/2.0;
    if(dY >= 0.0) return ((rand()%RANGE > RANGE/TWO) ? dY1 : -dY1);
  }
}

double estimateHashProb(double dStdDev,int iDiv=10) {
  int iNoOfHits = 0;
  double dDiv = 1.0/iDiv;
  for(int i=0;i<SAMPLE_SIZE;++i) {
    double dY = double(rand()%(int(dDiv*RANGE)))/RANGE;
    double dSample = genStdNormal();
    if(dSample <= (-dY+dDiv)/dStdDev && dSample >= -dY/dStdDev) 
      ++iNoOfHits;
  }
  return double(iNoOfHits)/SAMPLE_SIZE;
}

void createClusterSizes(int* iaCoord,double dBias,int iNoOfClusters,double v,double dProbNoise,ofstream& fOut,int RANGE=1000) {
  int iaRange[RANGE];
  for(int i = MINUS_ONE;++i < RANGE;) iaRange[i] = i;
  // generate P(i \in jth cluster) randomly with uniform bias dBias
  int iMinMinAlpha = int(max(dBias,(1.0-dProbNoise)/(ONE+v*(iNoOfClusters-ONE)))*RANGE);
  int iMaxMinAlpha = int(RANGE*(1.0-dProbNoise)/(v+iNoOfClusters-ONE));
  if(iMaxMinAlpha > iMinMinAlpha) {
    iaCoord[ZERO]=rand()%(iMaxMinAlpha-iMinMinAlpha)+iMinMinAlpha;
  } else iaCoord[ZERO]=iMinMinAlpha;
  int iRangeSize = int(RANGE*(1.0-dProbNoise)-iaCoord[ZERO]*(iNoOfClusters-TWO+v+ONE));
  if(iRangeSize<iNoOfClusters-TWO) {
    cerr << "ERROR: UNSUITABLE PARAMETERS: Decrease -b " << dBias << " or -a " << v << endl;
    exit(0);
  }
  random_sample_n(iaRange,iaRange+iRangeSize,iaCoord+ONE,iNoOfClusters-TWO);
  for(int j=ZERO;++j < iNoOfClusters-TWO;) iaCoord[j]=iaCoord[j+1];
  iaCoord[iNoOfClusters-TWO]=iRangeSize;
  sort(iaCoord+ONE,iaCoord+iNoOfClusters-ONE);
  for(int i=ZERO;++i<iNoOfClusters-ONE;)iaCoord[i] += ((i+1)*iaCoord[ZERO]);
  iaCoord[iNoOfClusters-ONE]=int((1.0-dProbNoise)*RANGE);
  int iTmp = ZERO;
  fOut << iNoOfClusters << ' ' << dProbNoise << ' ';
  for(int i=0;i < iNoOfClusters;++i) {
    fOut << double(iaCoord[i]-iTmp)/RANGE << ' ';
    iTmp=iaCoord[i];
  }
  fOut << endl;
}

void createClusterParams(int* iaMean,double *daStdDev,int iNoOfDims,int iNoOfClusters,DistMode distMode,int iDims,double dProbConstrain=0.5,double dProbOverlap=0.5) {
  int iaRange[RANGE];
  double r=2;
  int s=2;
  for(int iCntr = MINUS_ONE;++iCntr < iNoOfClusters;) {
    for(int d = MINUS_ONE;++d < iNoOfDims;) {
/*      int iScaleFactor = (rand()%s) + ONE;
      if(distMode==MVNORMAL || distMode==HYPER_RECT_UNIFORM)
      daStdDev[iNoOfDims*iCntr+d]=iScaleFactor*iScaleFactor*r*r;
      daStdDev[iNoOfDims*iCntr+d]=double(10+rand()%20);*/
      daStdDev[iNoOfDims*iCntr+d]=dWidth;
    }
  }
  for(int i = MINUS_ONE;++i < RANGE;) iaRange[i] = i;
  double x = 0.8;
  int iaCoord[iNoOfClusters];

  // generate the random full-dimensional centres
  switch(distMode) {
    case FULL_DIM:
      for(int d = MINUS_ONE;++d < iNoOfDims;) {
	// choose iNoOfClusters random points in range of dimension as centres
        random_sample_n(iaRange,iaRange+RANGE,iaCoord,iNoOfClusters);
        random_shuffle(iaCoord,iaCoord+iNoOfClusters);
        for(int iCntr = MINUS_ONE;++iCntr < iNoOfClusters;) {
          iaMean[iNoOfDims * iCntr + d] = iaCoord[iCntr];
        }
      }
      break;

    case HYPER_RECT_UNIFORM:
      for(int iCntr = MINUS_ONE;++iCntr < iNoOfClusters;) {
	double c = dProbConstrain*((iNoOfDatasets > 1)?(x+ (rand()%int((1.0-x)*1000))/500.0):1);
	double o = dProbOverlap*((iNoOfDatasets > 1)?(x+ (rand()%int((1.0-x)*1000))/500.0):1);
        double dO = ((iCntr > ZERO) ? (1.0-o)/(1.0-c) : 1.0); //~
        for(int d = MINUS_ONE;++d < iNoOfDims;) {
          if(iCntr > ZERO && iaMean[iNoOfDims*(iCntr-1)+d]!=MINUS_ONE) {
            if(rand()%100 < o*100)
              iaMean[iNoOfDims*iCntr+d]=iaMean[iNoOfDims*(iCntr-1)+d]+2*int(daStdDev[iNoOfDims*(iCntr-1)+d]*(((rand()%2)==0)?-1:1));//2 \sigma from prev mean
          } else {
            if(rand()%100 < c*dO*100)
              iaMean[iNoOfDims*iCntr+d]=rand()%RANGE;
          }
        }
      }
      break;

    case MVNORMAL:
      for(int iCntr = MINUS_ONE;++iCntr < iNoOfClusters;) {
	double c = dProbConstrain*((iNoOfDatasets > 1)?(x+ (rand()%int((1.0-x)*1000))/500.0):1);
	double o = dProbOverlap*((iNoOfDatasets > 1)?(x+ (rand()%int((1.0-x)*1000))/500.0):1);
        double dO = ((iCntr > ZERO) ? (1.0-o)/(1.0-c) : 1.0);
        for(int d = MINUS_ONE;++d < iNoOfDims;) {
          if(iCntr > ZERO && iaMean[iNoOfDims*(iCntr-1)+d]!=MINUS_ONE) {
            if(rand()%100 < o*100)
              iaMean[iNoOfDims*iCntr+d]=iaMean[iNoOfDims*(iCntr-1)+d]+2*int(daStdDev[iNoOfDims*(iCntr-1)+d]*(((rand()%2)==0)?-1:1));//2 \sigma from prev mean
//              iaMean[iNoOfDims*iCntr+d]=iaMean[iNoOfDims*(iCntr-1)+d]+(rand()%(2*int(daStdDev[iNoOfDims*(iCntr-1)+d])))*((rand()%2)?-1:1);
          } else {
            if(rand()%100 < c*dO*100)
              iaMean[iNoOfDims*iCntr+d]=rand()%RANGE;
          }
        }
      }
      break;

    default:
      cerr << "ERROR: UNDEFINED DATA GENERATION MODE " << distMode << endl;
      exit(0);
  }
#ifdef OLD_DEBUG
  for(int i=MINUS_ONE;++i < iNoOfClusters;) {
    for(int d=MINUS_ONE;++d < iDims;) cerr << iaMean[iNoOfDims*i+d] << ' ';
    cerr << endl;
  }
  for(int i=MINUS_ONE;++i < iNoOfClusters;) {
    for(int d=MINUS_ONE;++d < iDims;) cerr << daStdDev[iNoOfDims*i+d] << ' ';
    cerr << endl;
  }
#endif
}

void writeNoisyPoint(ofstream& fOut,int iNoOfDims,int iDims) {
  // for each dimension
  for(int d = MINUS_ONE;++d < iNoOfDims;) {
    int iRand = rand();
    if(d >= iDims) continue;
    fOut << iRand%RANGE;
    if(d != iDims - ONE) fOut << ' ';
  }
}

void writePoint(ofstream& fOut,int iNoOfDims,int iDims,int* iaMean,double *daStdDev,DistMode distMode) {
  double dRand;
  int iRand = 100;
  // for each dimension
  for(int d = MINUS_ONE;++d < iNoOfDims;) {
    dRand = genStdNormal();
    iRand=rand();
    if(d >= iDims) continue;
    if(iaMean[d]==MINUS_ONE) {
      fOut << iRand%RANGE;
    } else {
      switch(distMode) {
        case FULL_DIM:
          fOut << max(0,min(int(iaMean[d]+daStdDev[d]*dRand),RANGE));
          break;

        case HYPER_RECT_UNIFORM:
          fOut << max(0,min(RANGE,iaMean[d]-int(.5*daStdDev[d])+iRand%(int(daStdDev[d]))));
          break;

        case MVNORMAL:
          fOut << max(0,min(int(iaMean[d]+daStdDev[d]*dRand),RANGE));
          break;

        default:
          cerr << "ERROR: UNDEFINED DATA GENERATION MODE " << distMode << endl;
          exit(0);
      }
    }
    if(d != iDims - ONE) fOut << ' ';
  }
}

void writePt(ofstream& fOut,int iNoOfDims,int iDims,int* iaMean,double *daStdDev,int* iaCoord,DistMode distMode) {
  // choose a cluster
  int iRand = rand()%RANGE;
  int iClusterId = ZERO;
  while(iaCoord[iClusterId++] <= iRand);
  --iClusterId;
  writePoint(fOut,iNoOfDims,iDims,iaMean+iClusterId*iNoOfDims,daStdDev+iClusterId*iNoOfDims,distMode);
}

void parseArguments(int argc,char **argv) {

  // read args
  if(argc < NUM_ARGS) {
    cerr << "USAGE: dataGen\n";
    cerr << "\t-a<Max cluster size ratio   default="<<dMaxClusterSzRatio<<">\n";
    cerr << "\t-b<Min cluster size frac    default=" << dBias << ">\n";
    cerr << "\t-c<P(dim is bounded)        default=" << dProbConstrain << ">\n";
    cerr << "\t-d<Real Dims                default=" << iRealDims << ">\n";
    cerr << "\t-D<#(Dims)                  default=" << iNoOfDims << ">\n";
    cerr << "\t-e<generate subspaces only, no points>\n";
    cerr << "\t-m<mode                     default=" << iMode << ">\n";
    cerr << "\t-k<#(Clusters)              default=" << iNoOfClusters << ">\n";
    cerr << "\t-n<#(pts)                   default=" << iNoOfPts << ">\n";
    cerr << "\t-o<out file>\n";
    cerr << "\t-O<P(dim bounding overlaps) default=" << dProbOverlap << ">\n";
    cerr << "\t-r<noise                    default=" << dProbNoise << ">\n";
    cerr << "\t-R<range of variation       default=" << x << ">\n";
    cerr << "\t-s<#(datasets)              default=" << iNoOfDatasets << ">\n";
    cerr << "\t-w<width of clu in uni dist default=" << dWidth << ">\n";
    exit(0);
  }
  int c;
  while((c=getopt(argc,argv,"a:b:c:d:D:e:k:m:n:o:O:r:R:s:w:")) != -1) {
    switch(c) {
      case 'a':
        dMaxClusterSzRatio = atof(optarg);
        break;

      case 'b':
        dBias = atof(optarg);
        break;

      case 'c':
        dProbConstrain = atof(optarg);
        break;

      case 'd':
        iRealDims = atoi(optarg);
        break;

      case 'D':
        iNoOfDims = atoi(optarg);
        break;

      case 'e':
        bCentresOnly = true;
        break;

      case 'k':
        iNoOfClusters = atoi(optarg);
        break;

      case 'm':
	iMode = atoi(optarg);
        break;

      case 'n':
        iNoOfPts = atoi(optarg);
        break;

      case 'o':
	sprintf(caOutFile,"%s",optarg);
        break;

      case 'O':
	dProbOverlap = atof(optarg);
        break;

      case 'r':
        dProbNoise = atof(optarg);
        break;

      case 'R':
        x = atof(optarg);
        break;

      case 's':
        iNoOfDatasets = atoi(optarg);
        break;

      case 'w':
        dWidth = atof(optarg);
        break;

      default:
	cerr << "ERROR: OPTION -" << c << " NOT SUPPORTED\n";
	exit(0);
    }
  }
}

int main(int argc,char **argv) {
  parseArguments(argc,argv);
  DistMode distMode = DistMode(iMode);

  // test validity of parameters
  srand(20000);
  if(iRealDims > iNoOfDims) {
    cerr << " ERROR: INCORRECT DIMENSIONS : NoOfDims=" << iNoOfDims << " RealDims=" << iRealDims << endl;
    exit(0);
  }

  int k = (iNoOfDatasets > 1)? int((1.0+x)*iNoOfClusters):iNoOfClusters;
  int iaMean[k*iNoOfDims],iaMeanTmp[iNoOfDims];
  double daStdDev[k*iNoOfDims],daStdDevTmp[iNoOfDims];
  for(int i=ZERO;i<k*iNoOfDims;++i) {
    iaMean[i]=MINUS_ONE;
    daStdDev[i]=dStdDev;
  }
  createClusterParams(iaMean,daStdDev,iNoOfDims,k,distMode,iRealDims,dProbConstrain,dProbOverlap);
  for(int iDatasetId = ZERO;iDatasetId < iNoOfDatasets;++iDatasetId) {
    // open file to write to
    char caFile[1000],caOut[1000],caCentre[1000];
    memset(caOut,'\0',1000);
    if(iNoOfDatasets>ONE) {
      sprintf(caFile,"%.*ss%.2d%s",strlen(caOutFile)-3,caOutFile,iDatasetId,caOutFile+strlen(caOutFile)-3);
    } else sprintf(caFile,"%s",caOutFile);
    ofstream fHa(caFile,ios::out);
    if(fHa.bad()) {
      cerr << " ERROR: OPENING FILE " << caFile << endl;
      exit(0);
    }
    strncpy(caOut,caFile,strlen(caFile)-3);
    ofstream fOut(caOut,ios::out);
    if(fOut.bad()) {
      cerr << " ERROR: OPENING FILE " << caOut << endl;
      exit(0);
    }
    sprintf(caCentre,"%s_out",caOut);
    ofstream fCentre(caCentre,ios::out);
    if(fCentre.bad()) {
      cerr << " ERROR: OPENING FILE " << caCentre << endl;
      exit(0);
    }
    // perturb number of points and interesting subspaces in dataset
    int n = int(iNoOfPts*((iNoOfDatasets > 1)?x+ (rand()%int(1000*(1.0-x)))/500.0:1));
    int k = int(iNoOfClusters*((iNoOfDatasets > 1)?x+ (rand()%int(1000*(1.0-x)))/500.0:1));
    int iDims = iRealDims;//int(iRealDims*((iNoOfDatasets > 1)?x+ (rand()%int(1000*(1.0-x)))/500.0:1));
    int iaCoord[k];
    fCentre << iDims << endl;
    if(!bCentresOnly) {
      fOut << n << ' ' << iDims << ' ';
      fHa << iDims << endl;
      createClusterSizes(iaCoord,dBias,k,dMaxClusterSzRatio,dProbNoise,fOut);
    }
    int iPtsWritten = ZERO,iTmp = ZERO;
    for(int iCluId=ZERO;iCluId<k;++iCluId) { // for each cluster
      // determine number of points ascribed to that cluster
      int iNoOfCluPts = int(n*double(iaCoord[iCluId]-iTmp)/RANGE);
      iTmp = iaCoord[iCluId];
      for(int d=ZERO;d<iNoOfDims;++d) {
	// perturb the centre and deviation of the planted subspace
	if(iaMean[iNoOfDims*iCluId+d] < ZERO) {
	  iaMeanTmp[d]=MINUS_ONE;
          daStdDevTmp[d]=dStdDev;
	  continue;
	}
	double x = 0.02;
	iaMeanTmp[d]=max(0,min(RANGE,int(((iNoOfDatasets > 1)?iaMean[iCluId*iNoOfDims+d]-x*RANGE + 2*(rand()%int(RANGE*x)):iaMean[iCluId*iNoOfDims+d]))));
	daStdDevTmp[d]=daStdDev[iCluId*iNoOfDims+d]*((iNoOfDatasets > 1)?(1.0-x) + (rand()%int(1000*x))/500.0:1);
      }
      for(int d=MINUS_ONE;++d < iDims;) fCentre << iaMeanTmp[d] << ' ';
      fCentre << endl;
#ifdef DEBUG
//  for(int d=MINUS_ONE;++d < iDims;) fCentre << daStdDevTmp[d] << ' ';
//  fCentre << endl;
#endif
      if(!bCentresOnly) {
        for(int j=ZERO;j<iNoOfCluPts;++j) { // for each point ascribed to it
          if((iCluId+j) != ZERO) fHa << endl;
          writePoint(fHa,iNoOfDims,iDims,iaMeanTmp,daStdDevTmp,distMode);
        }
        iPtsWritten += iNoOfCluPts;
      }
    }
    if(!bCentresOnly) {
      while(++iPtsWritten <= n) {
        fHa << endl;
        writeNoisyPoint(fHa,iNoOfDims,iDims);
      }
      fHa << endl;
      fHa.flush();
    }
  }
  return 0;
}
