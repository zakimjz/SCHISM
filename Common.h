#define MAXDOUBLE   1.79769313486231570e+308

double dTotalCnt = ZERO;  // number of pow/mul operations during execution
double *daMp,*daU,*daL; // for maintaining bounds
double const DIFF = 0.00000000005; // to account for FP precision
int iIter = ZERO;

double distance(Point const& ptSrc,Point const& ptDest,double const dLpNorm=NORM,ClusterMode cMode=KCLUSTER) {
  double dDist = ZERO;

  // \forall d \in [1,NO_OF_DIMENSIONS]
  for(int iDimensionId=MINUS_ONE;++iDimensionId < ptSrc.getCoordSize();) {
    if(dLpNorm < ZERO) {
      dDist += ptSrc.getCoordinate(iDimensionId)*ptDest.getCoordinate(iDimensionId);
    } else {
      dDist += pow(fabs(ptSrc.getCoordinate(iDimensionId)-ptDest.getCoordinate(iDimensionId)),dLpNorm);
    }
  }
  dNoOfCalls += ONE;
  if(dLpNorm < ZERO) {
    if(dLpNorm == MINUS_ONE) dDist = acos(dDist);
  } else if(cMode < VCLUSTER) dDist = pow(dDist,1.0/dLpNorm);
  return dDist;
}

// idea from Jain paper - maximum movement
void updateMaxMoved(vector<Cluster>& vecCluSrc,ClusterMode cMode = KCLUSTER) {
  if(cMode == JCLUSTER) {
    // find cluster who's centre moved the most and second most
    int iMaxCluId=MINUS_ONE,iNumClusters=vecCluSrc.size();
    double dMaxMoved=ZERO;
    for(int j = MINUS_ONE; ++j < iNumClusters;) {
      if(dMaxMoved < vecCluSrc[j].getMoved()) {
        dMaxMoved = vecCluSrc[j].getMoved();
        iMaxCluId = j;
      }
    }
    double dMaxNextMoved=ZERO;
    for(int j = MINUS_ONE; ++j < iNumClusters;) {
      if(j == iMaxCluId) continue;
      if(dMaxNextMoved < vecCluSrc[j].getMoved()) {
        dMaxNextMoved = vecCluSrc[j].getMoved();
      }
    }
    for(int j = MINUS_ONE; ++j < iNumClusters;) {
      if(j == iMaxCluId) vecCluSrc[j].setMaxMoved(dMaxNextMoved);
      else               vecCluSrc[j].setMaxMoved(dMaxMoved);
    }
  }
}

void updateInterClusterDist(vector<Cluster>& vecCluSrc,ClusterMode cMode=VCLUSTER,double const dLpNorm=NORM) {
  double dDist,dDist2=ZERO;
  int iNoOfClusters = vecCluSrc.size();

  // \forall (j,l) \in [1,k] X [1,k]
  for(int j = MINUS_ONE; ++j < iNoOfClusters;) {
    int iFirst = j;
    if(cMode > VCLUSTER) iFirst += iIter*iNoOfClusters;
    for(int l = j; ++l < iNoOfClusters;) {
      int iSecond = l;
      switch(cMode) {
        case KCLUSTER:
          return;

        case OKCLUSTER:
          dDist = distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[l].getOtherCentre()),dLpNorm);
          dDist2 = distance(Point(vecCluSrc[l].getCentre()),Point(vecCluSrc[j].getOtherCentre()),dLpNorm);
          break;

        case ECLUSTER:
          dDist = distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[l].getCentre()),dLpNorm);
          break;

        case VCLUSTER:
          dDist = distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[l].getOtherCentre()),dLpNorm);
          dDist2 = distance(Point(vecCluSrc[l].getCentre()),Point(vecCluSrc[j].getOtherCentre()),dLpNorm);
          break;

        case IMP_VCLUSTER:
          dDist = distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[l].getOtherCentre()),dLpNorm);
          dDist2 = distance(Point(vecCluSrc[l].getCentre()),Point(vecCluSrc[j].getOtherCentre()),dLpNorm);
          iSecond += iIter*iNoOfClusters;
          break;

        case OPT_VCLUSTER:
          dDist = distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[l].getOtherCentre()),dLpNorm);
          dDist2 = distance(Point(vecCluSrc[l].getCentre()),Point(vecCluSrc[j].getOtherCentre()),dLpNorm);
          iSecond += iIter*iNoOfClusters;
          break;

        default:
          dDist = distance(Point(vecCluSrc[j].getOtherCentre()),Point(vecCluSrc[l].getOtherCentre()),dLpNorm,cMode);
          break;
      }
      vecCluSrc[j].setVecDist(iSecond,dDist);
      if(cMode == OKCLUSTER || cMode >= VCLUSTER)
        vecCluSrc[l].setVecDist(iFirst,dDist2);
      else 
        vecCluSrc[l].setVecDist(iFirst,dDist);
      if(cMode == JCLUSTER || cMode == ECLUSTER) {
        if(dDist<vecCluSrc[j].getVecDist(vecCluSrc[j].getNearestClusterId())) {
          vecCluSrc[j].setNearestClusterId(l);
        }
      }
    }
    if(cMode == JCLUSTER || cMode == ECLUSTER || cMode >= VCLUSTER) {
      dDist=distance(Point(vecCluSrc[j].getCentre()),Point(vecCluSrc[j].getOtherCentre()),dLpNorm);
      vecCluSrc[j].setMoved(dDist);
      vecCluSrc[j].setVecDist(iFirst,dDist);
    }
  }
  if(cMode == ECLUSTER) return;
  updateMaxMoved(vecCluSrc,cMode);
  for(int j = MINUS_ONE; ++j < iNoOfClusters;) {
    if(cMode < VCLUSTER) vecCluSrc[j].setCentre(vecCluSrc[j].getOtherCentre());
  }
}

// subspace clustering-related function
double comb(int n, int k) {
  if(k > n-k) {
    k = n - k;
  }

  double dNum = ONE;
  int i=n-k;
  while(++i<=n) {
    dNum *= i;
  }

  double dDenom = ONE;
  i=ONE;
  while(++i<=k) {
    dDenom *= i;
  }
  return dNum/dDenom;
}

// for vertical databases
void readDimension(ifstream& fData,int iNoOfPts,Dim& dimCurr,double* daMagnitude,double dLpNorm=NORM) {
  double dDim[iNoOfPts];
  fData.read((char*)dDim,sizeof(double)*iNoOfPts);
  if(NO_OF_POINTS > iNoOfPts || (int)fData.gcount() != int(sizeof(double)*iNoOfPts)) {
    if(NO_OF_POINTS > iNoOfPts)
      cerr << " ERROR: NOT ENOUGH POINTS : " << iNoOfPts << endl;
    else if (fData.gcount() != int(sizeof(double)*iNoOfPts))
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
    double daCoord[NO_OF_DIMENSIONS];
    fData.read((char*)daCoord,sizeof(double)*iNoOfDims);
    if(iNoOfDims < NO_OF_DIMENSIONS) {
      cerr << " ERROR: iNoOfDims = " << iNoOfDims << " NO_OF_DIMENSIONS = " << NO_OF_DIMENSIONS << endl;
      exit(0);
    }
    for(int iDimensionId=MINUS_ONE;++iDimensionId < NO_OF_DIMENSIONS;) {
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
    int iLen = strchr(caOffset,cDelim)-caOffset;
    if(abs(iLen) > MAX_COORD_SIZE) {
      pNew.setCoordinate(pNew.getCoordSize()-ONE,atof(caPoint+iOffset));
    } else if(iLen > ZERO) {
      memset(caCoord,'\0',MAX_COORD_SIZE);
      strncpy(caCoord,caOffset,iLen);
      pNew.setCoordinate(iDimensionId,atof(caCoord));
    } else {
      pNew.setCoordinate(iDimensionId,0);
    }
    iOffset += (iLen + ONE);
  }
  if(pMax != (Point *)NULL)pNew.normalise(*pMax,*pMin,dLpNorm);
  return true;
}

void readFileHorizontally(ifstream& fData,Point * ptArr,char cDelim=DELIMITER,ClusterMode cMode=KCLUSTER,double dLpNorm = NORM) {
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
  pMax = new Point(ptArr[ZERO]);
  pMin = new Point(ptArr[ZERO]);
  getMax(ptArr,*pMax);
  getMin(ptArr,*pMin);
  for(int i=ZERO;i < iNoOfPts;i++) ptArr[i].normalise(*pMax,*pMin,dLpNorm);
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
    readFileHorizontally(fData,ptArr,cDelim,cMode,dLpNorm);
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

*/#include "convertData.h"
#include "Algorithms.h"

// subspace clustering related function
bool prefixMatch(vector<int> const& vInt1,vector<int> const & vInt2,int d) {
  if(vInt1[d]!=vInt2[d] && equal(vInt1.begin(),vInt1.begin()+d,vInt2.begin()))
    return false;
  else
    return true;
}

// subspace clustering related function
bool mergeClusterSets(vector<vector<ClusterSet *> *>& vvCluSet,vector<ClusterSet *>* vecCluSet) {
  vector<ClusterSet *> *vecCluNew = new vector<ClusterSet *>();
  vector<int> vPoints;
  int d=(*vecCluSet)[ZERO]->getDimensions().size()-ONE;
//  cerr << endl;

  // for all ClusterSets in the current layer
  for(unsigned int i=ZERO;i < vecCluSet->size();++i) {
    int iMerged=ZERO;// cerr << endl << " i=" << i;

    // for all other ClusterSets in the current layer
    for(unsigned int j=i;++j < vecCluSet->size();) {
      if(prefixMatch((*vecCluSet)[i]->getDimensions(),(*vecCluSet)[j]->getDimensions(),d))
        continue;
      vPoints.resize((*vecCluSet)[i]->getPoints().size());
      vector<int>::iterator vIter=set_intersection((*vecCluSet)[i]->getPoints().begin(),(*vecCluSet)[i]->getPoints().end(),(*vecCluSet)[j]->getPoints().begin(),(*vecCluSet)[j]->getPoints().end(),vPoints.begin());
      vPoints.resize(vIter-vPoints.begin());
      if(vPoints.size()<SIX) continue;

      vector<int> vDim((*vecCluSet)[i]->getDimensions());
      vDim.push_back(((*vecCluSet)[j]->getDimensions())[vvCluSet.size()-ONE]);
      vector<ClusterSet*> vBase((*vecCluSet)[i]->getBase());
      vBase.push_back(((*vecCluSet)[j]->getBase())[vvCluSet.size()-ONE]);
      ClusterSet *csNew = new ClusterSet();
      csNew->setDimensions(vDim);
      csNew->setBase(vBase);
      csNew->setPoints(vPoints);
/*  cerr << " i= " << i << " j=" << j << ' ';
    for(int m=-1;++m < (*vecCluSet)[j]->getPoints().size();)
      cerr << ((*vecCluSet)[j]->getPoints())[m] << ' ';
    cout << endl;
    for(int m=-1;++m < vPoints.size();)
      cout << vPoints[m] << ' ';
    cout << endl;*/
      csNew->updateProb(*(vvCluSet[ZERO]),vvCluSet.size());
//      double dProb = 1 - pow((1-csNew->getProb()),comb(iPointCnt,vPoints.size()));
//      if(vBase.size() > 3 && csNew->getProb()>0.005) { delete csNew; continue; }
      ((*vecCluSet)[i])->addChild(csNew);
      ((*vecCluSet)[j])->addChild(csNew);
      vecCluNew->push_back(csNew);
      ++iMerged;
    }
  }
  if(vecCluNew->size() > ZERO) {
    vvCluSet.push_back(vecCluNew);
    return false;
  } else return true;
}

// subspace clustering related function
void produceClusters(vector<vector<ClusterSet *> *>& vvCluSet,vector<Cluster>& vecCluSrc) {
  int iNumber=ZERO;
  int i=vvCluSet.size();
  while(--i>MINUS_ONE) {
    for(unsigned int j=ZERO;j < vvCluSet[i]->size();++j) {
      ClusterSet& csCurr = *(*(vvCluSet[i]->begin()+j));
      if(csCurr.getChildren().size()==ZERO);
        csCurr.print();
      ++iNumber;
    }
  }
  cerr << " iNumber = " << iNumber << endl;
}

// subspace clustering related function
int hClusterData(int* iaCluster,vector<Cluster>& vecCluSrc,Point* ptArr) {
  Point pMax;
  getMax(ptArr,pMax);

  vector<vector<ClusterSet *> *> vvCluSet;
  vector<ClusterSet *>* vecCluSet = new vector<ClusterSet *>();
  vvCluSet.push_back(vecCluSet);
  int const iPointId = iPointCnt;
  int const iClusterId = iClusterCnt;
  Point pOne(ONE);
  Point* vecNewPt = new Point[NO_OF_POINTS];  // MUST FIX THIS TOO
  for(int d = MINUS_ONE; ++d < pMax.getCoordSize();) {
    createClusters(d,iaCluster,ptArr,vvCluSet[ZERO],vecNewPt,pMax,HCLUSTER);
  }
  iClusterCnt = iClusterId;
  iPointCnt = iPointId;
  int iLayer=ZERO;
  bool bConverged = false;
  while(!bConverged) {
    bConverged = mergeClusterSets(vvCluSet,vvCluSet[iLayer++]);
  }
  produceClusters(vvCluSet,vecCluSrc);
  return ZERO;
}

// incorporate pruning techniques from horizontal DB-based papers
void preprocess(int* iaIter,bool *baTest,vector<Cluster>& vecCluSrc,int* iaCluster,double* daSMinus,bool* babaTest,ClusterMode cMode,double dLpNorm=NORM) {
  int iNumClusters = vecCluSrc.size();
  double dDenom = 1.0/dLpNorm;
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) {
    int iCnt = ZERO,iOffset = iPt*iNumClusters;
    if(cMode > IMP_VCLUSTER) {
      daSUtmp[iPt]=daSUtotal[iPt];            // S(v^+,q^+)
      daSLtmp[iPt]=daSLtotal[iPt];            // T(q^+)
    }
    baTest[iPt]=false;
    if(iIter > ZERO) {
      double dTmp = daSMinus[iOffset+iaCluster[iPt]];
      ++iCnt;                               // min_iCntr (d(iCntr,iPt))
      dTmp = (dLpNorm > ZERO)?pow(dTmp,dDenom):acos(dTmp);              // find S^{1/p}
      double dLHS = dTmp + vecCluSrc[iaCluster[iPt]].getVecDist(iaCluster[iPt]+(iaIter[iPt])*iNumClusters);                   // LHS=max S(iPt,iaCluster[iPt])
      for(int iIterId = iaIter[iPt];++iIterId < iIter;) {
        dLHS += vecCluSrc[iaCluster[iPt]].getVecDist(iaCluster[iPt]+iIterId*iNumClusters);
      }
      int iCand = ZERO;
      for(int iCntr = MINUS_ONE; ++iCntr < iNumClusters;) {
//        double dRHS = max(fabs(dTmp-vecCluSrc[iCntr].getVecDist(iCntr)),fabs(dTmp - vecCluSrc[iaCluster[iPt]].getVecDist(iCntr+(iIter)*iNumClusters)));
        double dRHS = fabs(dTmp - vecCluSrc[iaCluster[iPt]].getVecDist(iCntr+(iaIter[iPt])*iNumClusters));
        for(int iIterId = iaIter[iPt];++iIterId < iIter;) {
          dRHS -= vecCluSrc[iCntr].getVecDist(iCntr+iIterId*iNumClusters);
        }
        if(dLHS < dRHS) {                     // if current iCntr isn't closest
          babaTest[iOffset+iCntr] = false;    // do not evaluate S^+ completely
        } else {                              // if other centre is a candidate
          ++iCand;
          if(iCntr != iaCluster[iPt]) {
            daSMinus[iOffset+iCntr] = ZERO;   // initialise S^-(iPt,iCntr) = 0
          }
          babaTest[iOffset+iCntr] = true;     // recalculate daSMinus
        }
      }
      if(iCand > ONE) {                       // if #(candidate centre) > 1
        baTest[iPt]=true;                     // evaluate candidates later
        daSMinus[iOffset+iaCluster[iPt]]=ZERO;// initialise S^-(iPt) = 0
      } else if (iCand  == ZERO) cerr << "ERROR: UNASSIGNED PT " << iPt << endl;
    } else {                                  // for first iteration
      baTest[iPt]=true;                       // all centres are candidates
      for(int iCntr = MINUS_ONE; ++iCntr < iNumClusters;)
        babaTest[iOffset+iCntr]=true;
    }
    dTotalCnt += iCnt;
  }
  if(iIter==ZERO)
    memset(daSMinus,0,sizeof(double)*NO_OF_POINTS*iNumClusters);
  return;
}

// check for convergence
int checkChanges(int *iaCluSize,int* iaCluster,bool* babaTest,bool* baTest,int iNumClusters,int* iaIter,double* daSMinus,ClusterMode cMode,double dLpNorm) {
  int iChanged = ZERO;
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) {
    if(baTest[iPt]==false) continue;
    // for each point, find closest centre
    int iMinCntr = MINUS_ONE;
    double dMinDist = MAXDOUBLE;
    int iOffset = iNumClusters*iPt;
    for(int iCntr = MINUS_ONE; ++iCntr < iNumClusters;) {
      if(babaTest[iOffset+iCntr]==false) continue;
      if(cMode < OPT_VCLUSTER) {
        double dTmp=daSMinus[iOffset+iCntr];
        if(dLpNorm==MINUS_ONE) dTmp = acos(dTmp); // SHOULD WE COUNT THIS?
        if(dMinDist > dTmp) {
          dMinDist = dTmp;
          iaIter[iPt]=iIter;
          iMinCntr = iCntr;
        }
      } else {
        iaIter[iPt]=iIter;
        iMinCntr = iCntr;
        break;
      }
    }
    if(iMinCntr == MINUS_ONE) {
      cerr << "ERROR : IN checkChanges() " << iPt << endl;
    }
    // reassign point to that centre
    if(iMinCntr != iaCluster[iPt]) {
      if(iaCluster[iPt] != MINUS_ONE)iaCluSize[iaCluster[iPt]]--;
      iaCluster[iPt]=iMinCntr;
      iaCluSize[iaCluster[iPt]]++;
      ++iChanged;
    }
#ifdef DEBUG
    cerr << iPt << "->" << iaCluster[iPt] << ' ';
#endif
  }
  return iChanged;
}

int reassignPts(ifstream& fData,Dim* dimCurr,Dim* dimArr,vector<Cluster>& vecCluSrc,int* iaIter,int *iaCluster,int* iaCluSize,double* daSMinus, bool* babaTest,ClusterMode cMode=VCLUSTER,double const dLpNorm=NORM) {
  int iChanged = ZERO;
  int iNumClusters=vecCluSrc.size();
  bool *baTest = new bool[NO_OF_POINTS];
  if(cMode > VCLUSTER) {
    preprocess(iaIter,baTest,vecCluSrc,iaCluster,daSMinus,babaTest,cMode,dLpNorm);
  } else {
    for(int i=MINUS_ONE;++i < iNumClusters*NO_OF_POINTS;) {
      if(dLpNorm < ZERO) daSMinus[i] = ZERO;
      else daSMinus[i] = ZERO;
      babaTest[i]=true;
    }
    for(int i=MINUS_ONE;++i < NO_OF_POINTS;)baTest[i]=true;
  }
  if(bRead == true) fData.read((char*)&iNoOfPts,sizeof(int));
  double dDenom2 = MAXDOUBLE;
  double *daMagnitude = new double[NO_OF_POINTS];
  memset(daMagnitude,0,sizeof(double)*NO_OF_POINTS);

  // for each dimension
  for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS;) {
    double daSmin[iNumClusters];
    double daSmax[iNumClusters];
    int iCnt = ZERO;
    if(bRead == true)readDimension(fData,iNoOfPts,*dimCurr,daMagnitude,dLpNorm);
    else             dimCurr = &(dimArr[d]);
    if(cMode > IMP_VCLUSTER) {
      if(dLpNorm > 0) {
        dDenom2=pow(NO_OF_DIMENSIONS-d-ONE,dLpNorm-ONE); ++iCnt;
      }
      for(int j = MINUS_ONE; ++j < iNumClusters;)
        vecCluSrc[j].getOtherCentre().setPlusMax(vecCluSrc[j].getOtherCentre().getPlusMax()-vecCluSrc[j].getOtherCentre().getCoordinate(d));         // T(v^+,d)
    }
    double dSPlusMax = MAXDOUBLE;

    // for each point in the dimension
    for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) {
      if(baTest[iPt]==false) continue;
      double dMinSmax = NO_OF_DIMENSIONS, dMaxSmin = ZERO;
      int iOffset = iPt *iNumClusters;
      if(cMode > IMP_VCLUSTER) {
        if(dLpNorm == MINUS_ONE) {
          daSUtmp[iPt] -= dimCurr->getCoordinate(iPt);   // T(q^+,d)
        } else {
          daSUtmp[iPt] -= pow(max(1.0-dimCurr->getCoordinate(iPt),dimCurr->getCoordinate(iPt)),dLpNorm);
          dSPlusMax = daSUtmp[iPt];
          daSLtmp[iPt] -= dimCurr->getCoordinate(iPt);   // T(q^+,d)
        }
        ++iCnt;
      }

      // for all cluster centres
      for(int iCntr = MINUS_ONE; ++iCntr < iNumClusters;) {
        if(babaTest[iOffset+iCntr]==false) continue;

        // S^-(iPt,iCntr,d) += ||D(d,iPt) - C(iCntr,d)||^p
        if(dLpNorm == MINUS_ONE) {
          daSMinus[iOffset+iCntr] += (dimCurr->getCoordinate(iPt)*vecCluSrc[iCntr].getOtherCentre().getCoordinate(d));
        } else {
          daSMinus[iOffset+iCntr] += pow(fabs(dimCurr->getCoordinate(iPt)-vecCluSrc[iCntr].getOtherCentre().getCoordinate(d)),dLpNorm);
        }
        ++iCnt;

        if(cMode > IMP_VCLUSTER) {

          // S_max(iPt,iCntr) = S^-(iPt,iCntr,d) + S^+_{max}(iPt,iCntr,d)
          daSmax[iCntr] = daSMinus[iOffset+iCntr];
          if(dLpNorm == MINUS_ONE) {   // my bounds for cosine similarity
            if(vecCluSrc[iCntr].getOtherCentre().getPlusMax() >= 1)
              daSmax[iCntr] += (daSUtmp[iPt] >= 1) ? ONE: daSUtmp[iPt];
            else daSmax[iCntr] += (daSUtmp[iPt]>=1)?vecCluSrc[iCntr].getOtherCentre().getPlusMax():vecCluSrc[iCntr].getOtherCentre().getPlusMax()*daSUtmp[iPt];
          } else daSmax[iCntr] += dSPlusMax;

          // S_min(iPt,iCntr) = S^-(iPt,iCntr,d) + S^+_{min}(iPt,iCntr,d)
          daSmin[iCntr] = daSMinus[iOffset+iCntr];
          if(dLpNorm != MINUS_ONE) {
            if(NO_OF_DIMENSIONS - d - ONE > ZERO) {
              daSmin[iCntr] += pow(fabs(vecCluSrc[iCntr].getOtherCentre().getPlusMax()-daSLtmp[iPt]),dLpNorm)/dDenom2; ++iCnt;
            }
          }
          if(dLpNorm == MINUS_ONE) {    // see SIGMOD BOND paper
            dMaxSmin = max(dMaxSmin,daSmin[iCntr]);
            if(daSmax[iCntr]+DIFF<dMaxSmin) babaTest[iOffset+iCntr]=false;
          } else {
            dMinSmax = min(dMinSmax,daSmax[iCntr]);
            if(daSmin[iCntr]>dMinSmax+DIFF) babaTest[iOffset+iCntr]=false;
          }
        }
      }
      if(cMode < OPT_VCLUSTER) continue;
      int iCand= ZERO;
      for(int iCntr = MINUS_ONE; ++iCntr < iNumClusters;) {
        if(babaTest[iOffset+iCntr]==false) continue;
        if(dLpNorm == MINUS_ONE) {
          if(daSmax[iCntr]+DIFF<dMaxSmin) babaTest[iOffset+iCntr]=false;
          else                               ++iCand;
        } else {
          if(daSmin[iCntr]>dMinSmax+DIFF) babaTest[iOffset+iCntr]=false;
          else                               ++iCand;
        }
      }
      if (iCand  == ZERO) cerr << "ERROR NO CENTRE CANDIDATES FOR " << iPt << ' ' << d << endl;
    }
    dTotalCnt += iCnt;
  }
  iChanged = checkChanges(iaCluSize,iaCluster,babaTest, baTest,iNumClusters,iaIter,daSMinus,cMode,dLpNorm);
  for(int j = MINUS_ONE; ++j < iNumClusters;) {
    vecCluSrc[j].setCentre(vecCluSrc[j].getOtherCentre());
  }
  return iChanged;
}

int vClusterData(char const* caDataFile,int* iaCluster,vector<Cluster>& vecCluSrc,Point* ptArr, Dim* dimArr,int iNumClusters,ClusterMode cMode=VCLUSTER, double const dLpNorm = NORM) {
  int iChanged = ONE;
  int *iaIter = new int[NO_OF_POINTS];
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) iaIter[iPt] = MINUS_ONE;
  int iaCluSize[iNumClusters];
  for(int j = MINUS_ONE; ++j < iNumClusters;) {
    iaCluSize[j] = ZERO;
  }
  double *daSMinus = new double[NO_OF_POINTS*iNumClusters];
  double *daMagnitude = new double[NO_OF_POINTS];
  memset(daMagnitude,0,sizeof(double)*NO_OF_POINTS);
  bool *babaTest = new bool[NO_OF_POINTS*iNumClusters];
  Dim *dimCurr;
  if(bRead == true) dimCurr = new Dim(NO_OF_POINTS);
  
  while(iChanged) {
    ifstream fData(caDataFile);
    iChanged = reassignPts(fData,dimCurr,dimArr,vecCluSrc,iaIter,iaCluster,iaCluSize,daSMinus,babaTest,cMode);
    for(int j = MINUS_ONE; ++j < iNumClusters;) { // initialise clusters
      vecCluSrc[j].getOtherCentre().setPlusMax(ZERO);
      for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS;) 
        vecCluSrc[j].getOtherCentre().setCoordinate(d,ZERO);
    }
    ifstream fData2(caDataFile);                  // recalculate centres
    if(bRead == true) fData2.read((char*)&iNoOfPts,sizeof(int));
    for(int d = MINUS_ONE; ++d < NO_OF_DIMENSIONS;) {
      if(bRead == true)readDimension(fData2,iNoOfPts,*dimCurr,daMagnitude,dLpNorm);
      else             dimCurr = &(dimArr[d]);
      for(int iPtId = MINUS_ONE; ++iPtId < NO_OF_POINTS;) {
        vecCluSrc[iaCluster[iPtId]].getOtherCentre().setCoordinate(d,vecCluSrc[iaCluster[iPtId]].getOtherCentre().getCoordinate(d)+dimCurr->getCoordinate(iPtId));
      }
      for(int j = MINUS_ONE; ++j < iNumClusters;) {
        if(iaCluSize[j] > ZERO) {
          double dCoord=vecCluSrc[j].getOtherCentre().getCoordinate(d)/iaCluSize[j];
          vecCluSrc[j].getOtherCentre().setPlusMax(vecCluSrc[j].getOtherCentre().getPlusMax()+dCoord);
          vecCluSrc[j].getOtherCentre().setCoordinate(d,dCoord);
        }
      }
    }
#ifdef DEBUG
    printClusters(vecCluSrc);
#endif
    if(iChanged == ZERO) break;
    if(dLpNorm < ZERO) {
      for(int iCntr = MINUS_ONE;++iCntr < iNumClusters;) {
        double dMagnitude = ZERO;
        for(int d= MINUS_ONE;++d < vecCluSrc[iCntr].getCentre().getCoordSize();)
          dMagnitude += vecCluSrc[iCntr].getCentre().getCoordinate(d)*vecCluSrc[iCntr].getCentre().getCoordinate(d);
        dMagnitude = sqrt(dMagnitude);
        for(int d= MINUS_ONE;++d < vecCluSrc[iCntr].getCentre().getCoordSize();)
          vecCluSrc[iCntr].getCentre().setCoordinate(d,vecCluSrc[iCntr].getCentre().getCoordinate(d)/dMagnitude);
      }
    }
    for(int j = MINUS_ONE; ++j < iNumClusters;)
      vecCluSrc[j].updateVecDist(iNumClusters*(iIter+ONE));
    updateInterClusterDist(vecCluSrc,cMode,dLpNorm);
    ++iIter;
  }
  delete [] daSMinus;
  delete [] babaTest;
  delete [] iaIter;
  return iIter;
}

void clusterData(char const* caDataFile,vector<Cluster>& vecCluSrc,Point* ptArr, Dim* dimArr,ClusterMode cMode=COMPARE, char cDelim=DELIMITER,int const iNumClusters = FIVE,double const dLpNorm=NORM) {

  // until convergence
  int iChanged = ONE,iIter=ZERO;
  int *iaCluster = new int[NO_OF_POINTS];
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) iaCluster[iPt] = MINUS_ONE;
  double *daMinDist = new double[NO_OF_POINTS];
  memset(daMinDist,0,sizeof(double)*NO_OF_POINTS);
  while(iChanged) {
    ifstream fData(caDataFile);
    switch(cMode) {
      case OKCLUSTER:    // obsolete
        iChanged = oKClusterData(daMinDist,iaCluster,vecCluSrc,ptArr,KCLUSTER);
        break;
      case KCLUSTER:    // regular K-means
        iChanged = kClusterData(fData,daMinDist,iaCluster,vecCluSrc,ptArr,KCLUSTER,cDelim);
        break;
      case COMPARE:     // sanity check
        iChanged = oKClusterData(daMinDist,iaCluster,vecCluSrc,ptArr);
/*        iChanged += kClusterData(fData,daMinDist,iaCluster,vecCluSrc2,ptArr,cDelim);
        for(int i = MINUS_ONE; ++i < NO_OF_POINTS;) {
          if(ptArr[i].getClusterId() != ptArr[i].getClusterId()) {
            cerr << " ERROR: Pt " << i << " assigned to clusters " << ptArr[i].getClusterId() << " and " << ptArr[i].getClusterId() << " in iter " << iIter << endl;
            exit(0);
          }
        }*/
        break;
      case HCLUSTER:     // subspace clustering
        iChanged = hClusterData(iaCluster,vecCluSrc,ptArr);
        break;
      case ECLUSTER:    // Elkan k-means clustering (not completely implemented)
        if(iIter == ZERO) {
          daU = new double[NO_OF_POINTS];
          daL = new double[NO_OF_POINTS*iNumClusters];
          memset(daL,0,sizeof(double)*NO_OF_POINTS*iNumClusters);
        }
        iChanged = kClusterData(fData,daMinDist,iaCluster,vecCluSrc,ptArr,ECLUSTER,cDelim);
        break;
      case JCLUSTER:    // Jain k-means clustering (not completely implemented)
        if(iIter == ZERO) {
          daMp = new double[NO_OF_POINTS];
          memset(daMp,0,sizeof(double)*NO_OF_POINTS);
        }
        iChanged = kClusterData(fData,daMinDist,iaCluster,vecCluSrc,ptArr,JCLUSTER,cDelim);
        break;
      case VCLUSTER:   // henceforth all vertical k-means algorithms
        iChanged = ZERO;
        iIter = vClusterData(caDataFile,iaCluster,vecCluSrc,ptArr,dimArr,iNumClusters,cMode);
        dNoOfCalls += (NO_OF_POINTS +int((iIter+dTotalCnt)/NO_OF_DIMENSIONS));
        break;
      case IMP_VCLUSTER:  // vertical K-means using pruning techniques
        iChanged = ZERO;
        iIter = vClusterData(caDataFile,iaCluster,vecCluSrc,ptArr,dimArr,iNumClusters,cMode);
        dNoOfCalls += (NO_OF_POINTS +int((iIter+dTotalCnt)/NO_OF_DIMENSIONS));
        break;
      case OPT_VCLUSTER:  // vertical K-means using pruning and BOND techniques
        iChanged = ZERO;
        iIter = vClusterData(caDataFile,iaCluster,vecCluSrc,ptArr,dimArr,iNumClusters,cMode);
        dNoOfCalls += (NO_OF_POINTS +int((iIter+dTotalCnt)/NO_OF_DIMENSIONS));
        break;
      default:
        cerr << " ERROR : INVALID MODE " << cMode;
        exit(0);
    }
    ++iIter;
  }
  delete[] iaCluster;
  delete[] daMinDist;
  if(cMode != COMPARE) {
    int iMode = int(cMode);
    if(bRead == true) iMode += 10;
    cout << NO_OF_POINTS << '\t' << iIter << "\t" << iMode << '\t' << NO_OF_DIMENSIONS << '\t' << vecCluSrc.size() << '\t' << dNoOfCalls << endl;
#ifdef DEBUG2
    printClusters(vecCluSrc);
#endif
  }
}
