#include<iostream>
#include<fstream>
#include<algorithm>
#include<utility>
#include<numeric>
#include<vector>
#include<list>
#include<cstdlib>
#include<cstdio>
#include<cmath>
using namespace std;

enum Digit {MINUS_ONE=-1,ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,TEN};
enum ClusterMode {ASC_TO_BIN=-2,HORZ_TO_VERT,OKCLUSTER,KCLUSTER,COMPARE,HCLUSTER,VHEIRARCHY,ECLUSTER,JCLUSTER,VCLUSTER,IMP_VCLUSTER,OPT_VCLUSTER};
const static int NUM_ARGS=NINE;
const static int MAX_LINE_SIZE = 16384;
const static int MAX_COORD_SIZE = TEN;
static char DELIMITER = ' ';
static int iPointCnt = ZERO;
static int iClusterCnt = ZERO;
static int NO_OF_DIMENSIONS = SIX;
static int START_DIMENSION = ZERO;
static int NO_OF_POINTS = 1000;
static double NORM = TWO;
static double dNoOfCalls=ZERO;
static bool bRead = false;
static int iNoOfPts;

#include "Point.h"
#include "Dim.h"
#include "Centre.h"
#include "Cluster.h"
#include "ClusterSet.h"
#include "Common.h"

void parseArgs(int argc,char **argv) {
  if(argc < NUM_ARGS) {
    cerr << "USAGE: verticalCluster <PointsFileName> <Delimiter> <# dimensions> <Norm> <Mode> <Start Dimension> <CentresFileName/k> <# points>" << endl;
    exit(0);
  }
}

void driver(int argc, char** argv) {
  parseArgs(argc,argv);
  int iArgNo = ZERO;
  char* caDataFile = argv[++iArgNo];
  DELIMITER = *(argv[++iArgNo]);
  NO_OF_DIMENSIONS = atoi(argv[++iArgNo]);
  NORM = atof(argv[++iArgNo]);
  ClusterMode cMode = ClusterMode(atoi(argv[++iArgNo]));
  if(cMode > 10) {
    cMode = ClusterMode(int(cMode) - 10);
    bRead = true;
  }
  START_DIMENSION = atoi(argv[++iArgNo]);
  if(cMode > 0) {
    Point* ptArr;
    Dim* dimArr;
    vector<Cluster> vecCluster;
    char* caNumClusters = argv[++iArgNo];
    NO_OF_POINTS = atoi(argv[++iArgNo]);
    if(cMode >= VCLUSTER) dimArr = new Dim[NO_OF_DIMENSIONS];
    else                  ptArr = new Point[NO_OF_POINTS];
    readData(caDataFile,ptArr,dimArr,cMode);
    createClusters(ptArr,vecCluster,dimArr,caNumClusters,cMode);
    if(bRead == true) {
      if(cMode >= VCLUSTER) delete[] dimArr;
      else                  delete[] ptArr;
    }
    clusterData(caDataFile,vecCluster,ptArr,dimArr,cMode,DELIMITER,atoi(caNumClusters));
  } else {
    char* caNewFile = argv[++iArgNo];
    NO_OF_POINTS = atoi(argv[++iArgNo]);
    if(cMode == HORZ_TO_VERT) {
      horizontalToVertical(caDataFile,caNewFile);
    } else {
      asciiToBinary(caDataFile,caNewFile);
    }
  }
}

int main(int argc, char **argv) {
  driver(argc,argv);
  return 0;
}
