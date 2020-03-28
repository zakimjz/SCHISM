#include<cstdio>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<utility>
#include<vector>
#include<unistd.h>
using namespace std;

typedef vector<double> Subspace;
typedef vector<Subspace> Dataset;
static int NUM_ARGS,iNoOfDatasets,iNoOfDimensions,MAX_LINE_SIZE=16484,MAX_COORD_SIZE=6,RANGE=10;
static double dSupport;
char caEclatCall[1000],caDataFile[1000];
enum Digit {MINUS_ONE=-1,ZERO,ONE,TWO,THREE};

bool readSubspace(ifstream& fData,Subspace& subspaceNew,char cDelim=' ') {
  char caSubspace[MAX_LINE_SIZE];
  if(fData.getline(caSubspace,MAX_LINE_SIZE).eof()) return false;
  if(strchr(caSubspace,cDelim) == (char *)NULL) {
    cerr << " ERROR: DID NOT FIND DELIMITER (" << cDelim << ")" << endl;
    exit(0);
  }
  int iOffset = ZERO,iLen=ZERO;
  char caCoord[MAX_COORD_SIZE];
  char* caOffset = caSubspace + iOffset;
  while(1) {
    caOffset = caSubspace + iOffset;
    if(strchr(caSubspace,cDelim) == (char *)NULL) {
      subspaceNew.push_back(atof(caSubspace+iOffset));
      break;
    } else if(iLen > ZERO) {
      iLen = strchr(caOffset,cDelim)-caOffset;
      memset(caCoord,'\0',MAX_COORD_SIZE);
      strncpy(caCoord,caOffset,iLen);
      subspaceNew.push_back(atof(caCoord));
      if(subspaceNew[subspaceNew.size()-1]<0)
	subspaceNew[subspaceNew.size()-1]=RANGE/2;
    } else {
      subspaceNew.push_back(0);
    }
    iOffset += (iLen + ONE);
  }
  return true;
}

void parseArgs(int argc, char** argv) {
  if(argc < NUM_ARGS) {
    cerr << "USAGE: matcher\n";
    cerr << "\t-d<#(dimensions)>\n";
    cerr << "\t-S<support fraction>\n";
    cerr << "\t-s<#(datasets)>\n";
    cerr << "\t-i<input file name>\n";
    exit(0);
  }
  int c;
  while((c=getopt(argc,argv,"i:s:S:")) != -1) {
    switch(c) {
      case 'd': iNoOfDimensions = atoi(optarg);  break;
      case 'i': sprintf(caDataFile,"%s",optarg); break;
      case 's': iNoOfDatasets = atoi(optarg);    break;
      case 'S': dSupport = atof(optarg);         break;
      default:  cerr << "ERROR: OPTION -" << c << " NOT SUPPORTED\n"; exit(0);
    }
  }
}

bool converged(){return true;}

void gibbsAlign(vector<Dataset>& vDataset) {
  while(!converged()) {
    for(int i=ZERO;i < vDataset.size();++i) {
      for(int i=ZERO;i < vDataset.size();++i) {
      
      }
    }
  }
}

void alignSubspaces(int iMode,vector<Dataset>& vDataset) {
  switch(iMode) {
    case 1: gibbsAlign(vDataset); break;
    default:
      cerr << " ERROR: ALIGN MODE " << iMode << " NOT SUPPORTED" << endl;
      exit(0);
  }
}

void normaliseSubspace(Subspace& s1) {
  double dSum = 0;
  for(int d=-1;++d < s1.size();) dSum += s1[d]*s1[d];
  double dNorm = sqrt(dSum);
  for(int d=-1;++d < s1.size();) s1[d] /= dNorm;
}

void multiplySubspaces(Dataset dataset1,Dataset dataset2, vector<Subspace> & vSubspaceProd) {
  vSubspaceProd.resize(dataset1.size());
  for(int i=MINUS_ONE;++i < dataset2.size();) normaliseSubspace(dataset2[i]);
  for(int i=MINUS_ONE;++i < dataset1.size();) {
    normaliseSubspace(dataset1[i]);
    vSubspaceProd[i].resize(dataset2.size());
    for(int j=MINUS_ONE;++j < dataset2.size();) {
      vSubspaceProd[i][j] = 0;
      for(int d=MINUS_ONE;++d < dataset1[i].size();) {
	vSubspaceProd[i][j] += dataset1[i][d]*dataset2[j][d];
      }
    }
  }
}

void driver() {
  char caInFile[1000];
  vector<Dataset> vDataset;
  for(int iDatasetId=ZERO;iDatasetId < iNoOfDatasets;++iDatasetId) {
    sprintf(caInFile,"%s.%d.ibm",caDataFile,iDatasetId);
    sprintf(caEclatCall,"eclat -i %s -a 2 -d 1 -s %f -f 1 2>tmp",caInFile,dSupport);
    system(caEclatCall);
    ifstream fTmp("tmp");
    if(fTmp.bad()) {
      cerr << "ERROR IN Tmp FILE" << endl;
      exit(0);
    }
    Dataset datasetNew;
    while(!fTmp.eof()) {
      Subspace subspaceNew;
      readSubspace(fTmp,subspaceNew);
      datasetNew.push_back(subspaceNew);
    }
    vDataset.push_back(datasetNew);
  }
  alignSubspaces(1,vDataset);
}

int main(int argc,char** argv) {
//  parseArgs(argc,argv);
//  driver();
  srand(1000);
  int i1[4] = {1,2,4,7},i2[4]={3,8,4,6},i3[4]={4,1,7,9},i4[4]={2,5,8,2};
  Subspace s1,s2,s3,s4; s1.resize(4); s2.resize(4); s3.resize(4); s4.resize(4);
  Dataset d1,d2,d3;
  for(int i=-1;++i < 4;) {
    s1[i]=i1[i]+rand()%3; s2[i]=i2[i]+rand()%3; s3[i]=i3[i]+rand()%3-1; s4[i]=i4[i]+rand()%3-1;
  }
  d1.push_back(s1); d1.push_back(s2); d1.push_back(s3); d1.push_back(s4);
  for(int i=-1;++i < 4;) {
    s1[i]=i1[i]+rand()%3; s2[i]=i2[i]+rand()%3; s3[i]=i3[i]+rand()%3-1; s4[i]=i4[i]+rand()%3-1;
  }
  d2.push_back(s2); d2.push_back(s1); d2.push_back(s4); d2.push_back(s3);
  multiplySubspaces(d1,d2,d3);
  for(int i=0;i<d1.size();++i) {
    for(int j=0;j<d2.size();++j) {
      cerr << d3[i][j] << "\t";
    }
    cerr << endl;
  }
  return 1;
}
