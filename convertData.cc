#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <numeric>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include<unistd.h>
#include <string>
#include "Point.h"
#include "Dim.h"
#include "Centre.h"
#include "convertData.h"

using namespace std;
using namespace __gnu_cxx;

char caNewCentreFile[1000],caCentreFile[1000],caNewFile[1000], caDataFile[1000];
bool bCentresOnly = false;
double *daSUtotal,*daSUtmp,*daSLtotal,*daSLtmp; // for maintaining bounds
Point *pMax, *pMin;


void parseArgs(int argc,char **argv) {
  if(argc < NUM_ARGS) {
    cerr << "USAGE: convertData\n";
    cerr << "\t-i<PtsFileName>\n";
    cerr << "\t-l<Delimiter                                      default=" << DELIMITER << ">\n";
    cerr << "\t-d<# dimensions                                   default=" << NO_OF_DIMENSIONS << ">\n";
    cerr << "\t-m<Mode MAFIA=-6,LDR=-5,WEKA=-4,IBM=-3,ASC_BIN=-2 default=" << ((int)IBM_DATA_FORMAT) << ">\n";
    cerr << "\t-M<Max value of dimension                         default=" << MAX << ">\n";
    cerr << "\t-L<Min value of dimension                         default=" << MIN << ">\n";
    cerr << "\t-o<NewFileName>                              default=<PtsFlNmPrefix>.ibm\n";
    cerr << "\t-n<# points                                       default=" << NO_OF_POINTS << ">\n";
    cerr << "\t-c<Centre file name>                         default=<PtsFlNmPrefix>_out\n";
    cerr << "\t-C<new centre file name>                     default=<PtsFlNmPrefix>_out.ibm\n";
    cerr << "\t-s<#(datasets)                                    default=" << iNoOfDatasets << ">\n";
    cerr << "\t-S<Start Dimension                                default=" << START_DIMENSION << ">\n";
    cerr << "\t-x<Divisions/dimension                            default=" << NORM << ">\n";
    exit(0);
  }
  int c;
  while((c=getopt(argc,argv,"c:C:d:e:i:l:L:m:M:n:o:s:S:x:")) != -1) {
    switch(c) {
      case 'c':
        sprintf(caCentreFile,"%s",optarg);
        break;

      case 'C':
        sprintf(caNewCentreFile,"%s",optarg);
        break;

      case 'd':
        NO_OF_DIMENSIONS = atoi(optarg);
        break;

      case 'e':
        bCentresOnly = true;
        break;

      case 'i':
        sprintf(caDataFile,"%s",optarg);
        break;

      case 'l':
        DELIMITER = *(optarg);
        break;

      case 'L':
        MIN = atof(optarg);
        break;

      case 'm':
        iMode = atoi(optarg);
        break;

      case 'M':
        MAX = atof(optarg);
        break;

      case 'n':
        NO_OF_POINTS = atoi(optarg);
        break;

      case 'o':
        sprintf(caNewFile,"%s",optarg);
        break;

      case 's':
        iNoOfDatasets = atoi(optarg);
        break;

      case 'S':
        START_DIMENSION = atoi(optarg);
        break;

      case 'x':
        NORM = atoi(optarg);
        break;

      default:
        cerr << "ERROR: OPTION -" << c << " NOT SUPPORTED\n";
        exit(0);
    }
  }
}

void setNames(char* caFileIn,char* caFileOut,char* caExt, int iDatasetId) {
  if(iNoOfDatasets>ONE) {
    sprintf(caFileIn,"%ss%.2d.ha",caDataFile,iDatasetId);
    sprintf(caFileOut,"%ss%.2d%s",caDataFile,iDatasetId,caExt);
    if(*(1 + strrchr(caFileIn,'/'))=='s') {
      sprintf(caCentreFile,"%ss%.2d_out",caDataFile,iDatasetId);
      sprintf(caNewCentreFile,"%s.ibm",caCentreFile);
    }
  } else {
    sprintf(caFileIn,"%s.ha",caDataFile);
    sprintf(caFileOut,"%s%s",caDataFile,caExt);
    if(*(1 + strrchr(caFileIn,'/'))=='s') {
      sprintf(caCentreFile,"%s_out",caDataFile);
      sprintf(caNewCentreFile,"%s.ibm",caCentreFile);
    }
  }
}


void readDimension(ifstream& fData,int iNoOfPts,Dim& dimCurr,double* daMagnitude,double dLpNorm=NORM) {
  double dDim[iNoOfPts];
  fData.read((char*)dDim,sizeof(double)*iNoOfPts);
  if(NO_OF_POINTS > iNoOfPts || (int)fData.gcount() != int(sizeof(double)*iNoOfPts)) {
    if(NO_OF_POINTS > iNoOfPts)
      cerr << " ERROR: NOT ENOUGH POINTS : " << iNoOfPts << endl;
    else if (fData.gcount() != (unsigned int)(sizeof(double)*iNoOfPts))
      cerr << " ERROR: ENOUGH POINTS NOT READ " << endl;
    exit(0);
  }
  dimCurr.setSize(NO_OF_POINTS);
  for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
    dimCurr.setCoordinate(iPt,dDim[iPt]);
    if(dLpNorm < ZERO) {
      daMagnitude[iPt] += dDim[iPt]*dDim[iPt];
    }
  }
  if(dLpNorm >= ZERO) {
    dimCurr.updateMin();
    dimCurr.updateMax();
    dimCurr.normalise();
    dimCurr.updateSum();
  }
}

void readFileVertically(ifstream& fData,Dim* dimArr,char cDelim=DELIMITER,ClusterMode cMode=VCLUSTER,double dLpNorm=NORM) {
  fData.read((char*)&iNoOfPts,sizeof(int));
  if(NO_OF_POINTS > iNoOfPts) {
    cerr << " ERROR: NOT ENOUGH POINTS : " << iNoOfPts << endl;
    exit(0);
  }
  if(cMode > IMP_VCLUSTER) {
    daSUtotal = new double[NO_OF_POINTS];
    daSUtmp = new double[NO_OF_POINTS];
    daSLtotal = new double[NO_OF_POINTS];
    daSLtmp = new double[NO_OF_POINTS];
    memset(daSUtotal,0,sizeof(double)*NO_OF_POINTS);
    memset(daSLtotal,0,sizeof(double)*NO_OF_POINTS);
  }
  double* daMagnitude;
  if(dLpNorm == MINUS_ONE) {
    daMagnitude = new double[NO_OF_POINTS];
    memset(daMagnitude,0,sizeof(double)*NO_OF_POINTS);
  }
  for(int d = MINUS_ONE;++d < NO_OF_DIMENSIONS;) {
    readDimension(fData,iNoOfPts,dimArr[d],daMagnitude,dLpNorm);
    if(cMode < OPT_VCLUSTER) continue;
    for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
      if(dLpNorm < ZERO) {
        daSUtotal[iPt] += dimArr[d].getCoordinate(iPt);   // T(q^+),used in upper bound computations
      } else {
        daSUtotal[iPt] += pow(max(1.0-dimArr[d].getCoordinate(iPt),dimArr[d].getCoordinate(iPt)),dLpNorm);
        daSLtotal[iPt] += dimArr[d].getCoordinate(iPt);   // T(q^+)
      }
    }
  }
  if(dLpNorm < ZERO) {    // convert to unit vector for vertical DB
    for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
      daMagnitude[iPt] = sqrt(daMagnitude[iPt]);
    }
    for(int d = MINUS_ONE;++d < NO_OF_DIMENSIONS;) {
      for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
        dimArr[d].setCoordinate(iPt,dimArr[d].getCoordinate(iPt)/daMagnitude[iPt]);
        if(cMode > IMP_VCLUSTER) daSLtotal[iPt] += dimArr[d].getCoordinate(iPt);
      }
    }
  }
//  sort(dimArr,dimArr+NO_OF_DIMENSIONS);    MUST EXPLORE THIS
}

void getMax(Point const& ptSrc,list<Centre *> const& listPt,Point& pMax) {
  for(list<Centre *>::const_iterator lIter=listPt.begin();lIter!=listPt.end();++lIter) {
    for(int d = MINUS_ONE; ++d < pMax.getCoordSize();) {
      pMax.setCoordinate(d,max(pMax.getCoordinate(d),(*lIter)->getCoordinate(d)));
    }
  }
  for(int d = MINUS_ONE; ++d < pMax.getCoordSize();) {
    pMax.setCoordinate(d,max(pMax.getCoordinate(d),ptSrc.getCoordinate(d)));
  }
}

void getMax(Point const* ptArr,Point& pMax) {
  for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
    for(int d = MINUS_ONE; ++d < pMax.getCoordSize();) {
      pMax.setCoordinate(d,max(pMax.getCoordinate(d),ptArr[iPt].getCoordinate(d)));
    }
  }
}

void getMin(Point const* ptArr,Point& pMin) {
  for(int iPt=MINUS_ONE;++iPt < NO_OF_POINTS;) {
    for(int d = MINUS_ONE; ++d < pMin.getCoordSize();) {
      pMin.setCoordinate(d,min(pMin.getCoordinate(d),ptArr[iPt].getCoordinate(d)));
    }
  }
}

bool readPtHorizontally(ifstream& fData,Point& pNew,char cDelim=DELIMITER,double dLpNorm=NORM,int iNoOfDims=NO_OF_DIMENSIONS) {
  if(cDelim == 'b') {    // while reading from binary files
    double daCoord[iNoOfDims];
    fData.read((char*)daCoord,sizeof(double)*iNoOfDims);
    if(iNoOfDims < NO_OF_DIMENSIONS) {
      cerr << " ERROR: iNoOfDims = " << iNoOfDims << " NO_OF_DIMENSIONS = " << NO_OF_DIMENSIONS << endl;
      exit(0);
    }
    for(int iDimensionId=MINUS_ONE;++iDimensionId < iNoOfDims;) {
      pNew.setCoordinate(iDimensionId,daCoord[iDimensionId]);
    }
    return true;
  }
  char caPoint[MAX_LINE_SIZE];
  if(fData.getline(caPoint,MAX_LINE_SIZE).eof()) return false;
  if(strchr(caPoint,cDelim) == (char *)NULL) {
    cerr << " ERROR: DID NOT FIND DELIMITER (" << cDelim << ")" << endl;
    exit(0);
  }
  int iOffset = ZERO;
  char caCoord[MAX_COORD_SIZE];
  for(int iDimensionId=MINUS_ONE;++iDimensionId < START_DIMENSION;) {
    iOffset += (strchr(caPoint + iOffset,cDelim) - (caPoint + iOffset) + ONE);
  }

  // \forall d \in [START_DIMENSION,NO_OF_DIMENSIONS]
  for(int iDimensionId=MINUS_ONE;++iDimensionId < pNew.getCoordSize();) {
    char* caOffset = caPoint + iOffset;
    if(strchr(caOffset,cDelim) == (char *)NULL) {
      pNew.setCoordinate(pNew.getCoordSize()-ONE,atof(caPoint+iOffset));
    } else {
      int iLen = strchr(caOffset,cDelim)-caOffset;
//      if(iLen > ZERO) {
        memset(caCoord,'\0',MAX_COORD_SIZE);
        strncpy(caCoord,caOffset,iLen);
        pNew.setCoordinate(iDimensionId,atof(caCoord));
//      } else {
//        pNew.setCoordinate(pNew.getCoordSize()-ONE,atof(caPoint+iOffset));
//        pNew.setCoordinate(iDimensionId,0);
//      }
      iOffset += (iLen + ONE);
    }
  }
  if(pMax != (Point *)NULL)pNew.normalise(*pMax,*pMin,dLpNorm);
  return true;
}

void readFileHorizontally(ifstream& fData,Point * ptArr,char cDelim=DELIMITER,double dLpNorm = NORM) {
  int iNoOfDims;
  if(cDelim=='b') fData.read((char*)&iNoOfDims,sizeof(int));
  for(int iPt = MINUS_ONE;++iPt < NO_OF_POINTS;) {
    if(!readPtHorizontally(fData,ptArr[iPt],cDelim,dLpNorm,iNoOfDims)) break;
  }
  iNoOfPts = NO_OF_POINTS;
  if(cDelim=='b') return;
  if(dLpNorm < ZERO) {
    for(int i=ZERO;i < iNoOfPts;i++) ptArr[i].normalise(*pMax,*pMin,dLpNorm);
    return;
  }
  if(pMax == (Point *)NULL) {              // if we're not reading Centre file
    pMax = new Point(ptArr[ZERO]);
    pMin = new Point(ptArr[ZERO]);
    getMax(ptArr,*pMax);
    getMin(ptArr,*pMin);
    for(int i=ZERO;i < iNoOfPts;i++) ptArr[i].normalise(*pMax,*pMin,dLpNorm);
  }
}

void readData(char const* caDataFile,Point* ptArr,Dim* dimArr,ClusterMode cMode=KCLUSTER,char cDelim=DELIMITER,double dLpNorm=NORM) {
  // open data file
  ifstream fData(caDataFile);
  if(fData.bad()) {
    cerr << "ERROR IN DATA FILE " << caDataFile << endl;
    exit(0);
  }
  if(cMode >= VCLUSTER)
    readFileVertically(fData,dimArr,cDelim,cMode,dLpNorm);
  else
    readFileHorizontally(fData,ptArr,cDelim,dLpNorm);
}

// converts horizontal DB to binary vertical DB
void horizontalToVertical(char const* caHorizontalFile,char const *caVerticalFile, char cDelimiter=DELIMITER, double dLpNorm=NORM) {
  int iStartDim = START_DIMENSION, iNumDim = START_DIMENSION+NO_OF_DIMENSIONS;
  NO_OF_DIMENSIONS = min(NO_OF_DIMENSIONS,25); // read in batches
  ofstream fVertical(caVerticalFile,ios::out|ios::binary);
  while(START_DIMENSION < iNumDim) {           // for each dimension to be read
    Point* ptArr = new Point[NO_OF_POINTS];
    Dim* dimArr;
    readData(caHorizontalFile,ptArr,dimArr,HORZ_TO_VERT,cDelimiter);
    if(START_DIMENSION == iStartDim) {
      fVertical.write((char *)(&iNoOfPts),sizeof(int));
      NO_OF_POINTS = min(NO_OF_POINTS,iNoOfPts);
    }
    for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS; ) {
      for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
        double dCoord = ptArr[iPt].getCoordinate(d);
        fVertical.write((char *)(&dCoord),sizeof(double));
      }
    }
    START_DIMENSION += NO_OF_DIMENSIONS;       // update batch parameters
    NO_OF_DIMENSIONS = min(iNumDim-START_DIMENSION,NO_OF_DIMENSIONS);
    delete[] ptArr;
  }
}

// converts horizontal DB to binary vertical DB
void asciiToBinary(char const* caHorizontalFile,char const *caVerticalFile, char cDelimiter=DELIMITER, double dLpNorm=NORM) {
  ofstream fVertical(caVerticalFile,ios::out|ios::binary);
  if(fVertical.bad()) {
    cerr << "ERROR IN OUTPUT FILE " << caVerticalFile << endl;
    exit(0);
  }
  fVertical.write((char *)(&NO_OF_DIMENSIONS),sizeof(int));
  Point* ptArr = new Point[NO_OF_POINTS];
  Dim* dimArr;
  readData(caHorizontalFile,ptArr,dimArr,ASC_TO_BIN,cDelimiter);
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
    for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS; ) {
      double dCoord = ptArr[iPt].getCoordinate(d);
      fVertical.write((char *)(&dCoord),sizeof(double));
    }
  }
  fVertical.close();
  delete[] ptArr;
}

// converts to ascii Weka data format. Use makebin after this to produce binary
void horzAsciiToWekaHash(char const*caHorizontalFile,char const*caWekaFile,char cDelimiter=DELIMITER,double dLpNorm=NORM) {
  ofstream fVertical(caWekaFile,ios::out);
  if(fVertical.bad()) {
    cerr << "ERROR IN OUTPUT FILE " << caWekaFile << endl;
    exit(0);
  }
  Point ptCurr;
  ifstream fData(caHorizontalFile);
  if(fData.bad()) {
    cerr << "ERROR IN DATA FILE " << caHorizontalFile << endl;
    exit(0);
  }
  if(pMax == (Point *)NULL) {              // if we're not reading Centre file
    pMax = new Point(ptCurr);
    pMin = new Point(ptCurr);
    for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS; ) {
      pMax->setCoordinate(d,MAX);
      pMin->setCoordinate(d,0.0);
    }
  }
  iNoOfPts = NO_OF_POINTS;
  fVertical << "@RELATION arff"<< endl;
  for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS; ) {
    fVertical << "@attribute a" << d << " numeric" << endl;
  }
  fVertical << "@DATA"<< endl;
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
    if(!readPtHorizontally(fData,ptCurr,cDelimiter,dLpNorm,NO_OF_DIMENSIONS)) break;
    for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS; ) {
      if(d != ZERO) fVertical << ',';
      if (ptCurr.getCoordinate(d)<ZERO) fVertical << "-1";
      else fVertical << min(int((d+1)*dLpNorm-1),int(d*dLpNorm)+int(floor((ptCurr.getCoordinate(d))*dLpNorm)));
    }
    fVertical << endl;
  }
  fVertical.close();
}

// converts to ascii LDR data format. Use makebin after this to produce binary
void horzAsciiToLDRHash(char const*caHorizontalFile,char const*caLDRFile,char cDelimiter=' ',double dLpNorm=NORM) {
  ofstream fLDR(caLDRFile,ios::out);
  if(fLDR.bad()) {
    cerr << "ERROR IN OUTPUT FILE " << caLDRFile << endl;
    exit(0);
  }
  Point* ptArr = new Point[NO_OF_POINTS];
  Dim* dimArr;
  readData(caHorizontalFile,ptArr,dimArr,HORZ_TO_VERT,cDelimiter);
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
    fLDR << (iPt+ONE);
    int d = MINUS_ONE;
    while(++d < START_DIMENSION) ;
    --d;
    while(++d < NO_OF_DIMENSIONS) {
      if (ptArr[iPt].getCoordinate(d)<ZERO) fLDR << " -1";
      else fLDR << ' ' << ptArr[iPt].getCoordinate(d);
    }
    fLDR << endl;
  }
  fLDR.close();
//  delete ptCurr;
}

// converts from ascii to horizontal binary (MAFIA) data format
void horzAsciiToBinHash(char const*caHorizontalFile,char const*caBinFile,char cDelimiter=' ',double dLpNorm=NORM) {
  /* ofstream fBin(caBinFile,ios::out||ios::binary); */
  ofstream fBin(caBinFile,ios::binary);
  if(fBin.bad()) {
    cerr << "ERROR IN OUTPUT FILE " << caBinFile << endl;
    exit(0);
  }
  ifstream fData(caHorizontalFile);
  if(fData.bad()) {
    cerr << "ERROR IN DATA FILE " << caHorizontalFile << endl;
    exit(0);
  }
  int iNoOfDims = NO_OF_DIMENSIONS;
  fData >> iNoOfDims;
  char caDummy[30];
  fData.getline(caDummy,30);
  if(NO_OF_DIMENSIONS != iNoOfDims) {
    pMax=NULL;
    pMin=NULL;
  }
  NO_OF_DIMENSIONS = iNoOfDims;
  Point ptCurr;
  if(pMax == (Point *)NULL) {              // if we're not reading Centre file
    pMax = new Point(ptCurr);
    pMin = new Point(ptCurr);
    for(int d = MINUS_ONE; ++d < iNoOfDims; ) {
      pMax->setCoordinate(d,MAX);
      pMin->setCoordinate(d,MIN);
    }
  }
  iNoOfPts = NO_OF_POINTS;
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
    if(!readPtHorizontally(fData,ptCurr,cDelimiter,dLpNorm,iNoOfDims)) break;
    for(int d = MINUS_ONE; ++d < iNoOfDims; ) {
      int iVal = (int)(ptCurr.getCoordinate(d));
      fBin.write((char *)&iVal,sizeof(int));
    }
  }
  fBin.close();
}

// converts to ascii IBM data format. Use makebin after this to produce binary
void horzAsciiToIBMHash(char const*caHorizontalFile,char const*caIBMFile,char cDelimiter=' ',double dLpNorm=NORM) {
  ofstream fIBM(caIBMFile,ios::out);
  if(fIBM.bad()) {
    cerr << "ERROR IN OUTPUT FILE " << caIBMFile << endl;
    exit(0);
  }
  ifstream fData(caHorizontalFile);
  if(fData.bad()) {
    cerr << "ERROR IN DATA FILE " << caHorizontalFile << endl;
    exit(0);
  }
  int iNoOfDims = NO_OF_DIMENSIONS;
  fData >> iNoOfDims;
  char caDummy[30];
  fData.getline(caDummy,30);
  if(NO_OF_DIMENSIONS != iNoOfDims) {
    pMax=NULL;
    pMin=NULL;
  }
  NO_OF_DIMENSIONS = iNoOfDims;
  Point ptCurr;
  if(pMax == (Point *)NULL) {              // if we're not reading Centre file
    pMax = new Point(ptCurr);
    pMin = new Point(ptCurr);
    for(int d = MINUS_ONE; ++d < iNoOfDims; ) {
      pMax->setCoordinate(d,MAX);
      pMin->setCoordinate(d,MIN);
    }
  }
  iNoOfPts = NO_OF_POINTS;
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) {
    if(!readPtHorizontally(fData,ptCurr,cDelimiter,dLpNorm,iNoOfDims)) break;
    fIBM << (iPt+ONE) << ' ' << (iPt+ONE) << ' ' << iNoOfDims;
    for(int d = MINUS_ONE; ++d < iNoOfDims; ) {
      if (ptCurr.getCoordinate(d)<ZERO) fIBM << " -1";
      else fIBM << ' ' << min(int((d+1)*dLpNorm-1),int(d*dLpNorm)+int(floor((ptCurr.getCoordinate(d))*dLpNorm)));
    }
    fIBM << endl;
  }
  fIBM.close();
}



void driver(int argc, char** argv) {
  parseArgs(argc,argv);
  ClusterMode cMode = ClusterMode(iMode);
  char caFileIn[1000],caFileOut[1000];
  switch(cMode) {
    case MAFIA_FORMAT:
      for(int iDatasetId = ZERO; iDatasetId < iNoOfDatasets;++iDatasetId) {
        setNames(caFileIn,caFileOut,".maf",iDatasetId);
        if(!bCentresOnly) horzAsciiToBinHash(caFileIn,caFileOut);
        if(*(1+strrchr(caFileIn,'/'))=='s')
          horzAsciiToIBMHash(caCentreFile,caNewCentreFile);
      }
      break;

    case LDR_FORMAT:
      for(int iDatasetId = ZERO; iDatasetId < iNoOfDatasets;++iDatasetId) {
        setNames(caFileIn,caFileOut,".ldr",iDatasetId);
        if(!bCentresOnly) horzAsciiToLDRHash(caFileIn,caFileOut);
        if(*(1+strrchr(caFileIn,'/'))=='s')
          horzAsciiToLDRHash(caCentreFile,caNewCentreFile);
      }
      break;

    case WEKA_FORMAT:
      for(int iDatasetId = ZERO; iDatasetId < iNoOfDatasets;++iDatasetId) {
        setNames(caFileIn,caFileOut,".weka",iDatasetId);
        if(!bCentresOnly) horzAsciiToWekaHash(caFileIn,caFileOut);
        if(*(1+strrchr(caFileIn,'/'))=='s')
          horzAsciiToWekaHash(caCentreFile,caNewCentreFile);
      }
      break;

    case IBM_DATA_FORMAT:
      for(int iDatasetId = ZERO; iDatasetId < iNoOfDatasets;++iDatasetId) {
        setNames(caFileIn,caFileOut,".ibm",iDatasetId);
        if(!bCentresOnly) horzAsciiToIBMHash(caFileIn,caFileOut);
        if(*(1+strrchr(caFileIn,'/'))=='s')
          horzAsciiToIBMHash(caCentreFile,caNewCentreFile);
      }
      break;

    case ASC_TO_BIN:
      asciiToBinary(caDataFile,caNewFile);
      break;

    case HORZ_TO_VERT:
      horizontalToVertical(caDataFile,caNewFile);
      break;

    default:
      cerr << "ERROR: NO SUCH CONVERSION MODE " << cMode << endl;
      exit(0);
  }
}


int main(int argc, char **argv) {
  driver(argc,argv);
  return 1;
}
