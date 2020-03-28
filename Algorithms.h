void selectCentres(Point* ptArr,Dim* dimArr, vector<Point>& vecCntrPt,int k,ClusterMode cMode,double dLpNorm=NORM) {
  srand(1000);
#ifdef RANDOM                  // Centres are randomly selected points in DB
  int iaCentre[k];
  for(int i = MINUS_ONE;++i < k;) {
    iaCentre[i] = rand()%NO_OF_POINTS;
  }
  if(cMode < VCLUSTER) {
    for(int i=MINUS_ONE;++i<k;)
      vecCntrPt.push_back(ptArr[iaCentre[i]]);
  } else {
    for(int i=MINUS_ONE;++i<k;) {
      Point pNew;
      for(int d = MINUS_ONE; ++d < pNew.getCoordSize();) {
        pNew.setCoordinate(d,dimArr[d].getCoordinate(iaCentre[i]));
      }
      vecCntrPt.push_back(pNew);
    }
  }
#else                        // idea from Elkan pg 3 2nd last para
  Point pMean;               // first centre is mean of all points in DB
  for(int d = MINUS_ONE; ++d < pMean.getCoordSize(); ) {
    pMean.setCoordinate(d,ZERO);
    if(cMode < VCLUSTER) {          // horizontal DB
      for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; )
        pMean.setCoordinate(d,pMean.getCoordinate(d)+ptArr[iPt].getCoordinate(d));
      pMean.setCoordinate(d,pMean.getCoordinate(d)/NO_OF_POINTS);
    } else {
      pMean.setCoordinate(d,dimArr[d].getSum()/NO_OF_POINTS);
    }
  }
  vecCntrPt.push_back(pMean);
  vecCntrPt[0].setPointId(NO_OF_POINTS);
  vecCntrPt[0].normalise(*pMax,*pMin,dLpNorm);

  for(int i=ZERO;++i<k;) {       // to select the remaining centres
    double dMaxDist = ZERO,dMinDist=ZERO,dDist=ZERO;
    int iCntrPt = MINUS_ONE;
    Point pNew;
    for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS; ) { // for each point
      if(cMode < VCLUSTER) {    // find d(iPt, iCntr)
        dMinDist = distance(vecCntrPt[ZERO],ptArr[iPt],NORM);
      } else {
        for(int d = MINUS_ONE; ++d < pNew.getCoordSize();) {
          pNew.setCoordinate(d,dimArr[d].getCoordinate(iPt));
        }
        dMinDist = distance(vecCntrPt[ZERO],pNew,NORM);
      }
      for(unsigned int iCntr = ZERO; ++iCntr < vecCntrPt.size();) {
        if(cMode < VCLUSTER) {                     // horizontal DB
          dDist = distance(vecCntrPt[iCntr],ptArr[iPt],NORM);
        } else {                                   // vertical DB
          dDist = distance(vecCntrPt[iCntr],pNew,NORM);
        }
        dMinDist = min(dMinDist,dDist);   // find dist from nearest cntr
      }
      if(dMaxDist < dMinDist) {
        dMaxDist = dMinDist;
        iCntrPt = iPt;  // new centre = argmax_{iPt} min_{iCntr} d(iCntr,iPt)
      }
    }// cerr << iCntrPt << ' ';
    if(cMode < VCLUSTER) {
      vecCntrPt.push_back(ptArr[iCntrPt]);
    } else {
      for(int d = MINUS_ONE; ++d < pNew.getCoordSize();) {
        pNew.setCoordinate(d,dimArr[d].getCoordinate(iCntrPt));
      }
      vecCntrPt.push_back(pNew);
    }
    vecCntrPt[i].setPointId(NO_OF_POINTS+i);
  }
#endif
}

void createClusters(Point* ptArr,vector<Cluster>& vecCluster,Dim* dimArr,char* caCentreFile,ClusterMode cMode) {
  vector<Point> vecCntr;
  iClusterCnt=ZERO;
  if(strlen(caCentreFile)>THREE) {  // if file with initial centres is provided
//    readData(caCentreFile,vecCntr,dimArr,KCLUSTER);  MUST FIX THIS CASE
  } else {
    selectCentres(ptArr,dimArr,vecCntr,atoi(caCentreFile),cMode);
  }
  vector<Centre> vecCentre;
  int iNumClusters = vecCntr.size();
  for(int i=MINUS_ONE;++i < iNumClusters;) {
    vecCentre.push_back(Centre(vecCntr[i]));
  }
  // \forall j \in [1,k]
  for(int j = MINUS_ONE;++j < iNumClusters;) {
    Cluster cNew(vecCntr[j]);
    vecCluster.push_back(cNew);
  }
}

// this is a sub-space clustering related function
void createClusters(int d,int* iaCluster,Point* ptArr,vector<ClusterSet *> * vecCluSet,Point* ptArrNew,Point & pMax,ClusterMode cMode) {
  int iNoOfPts = NO_OF_POINTS;
  for(int i=MINUS_ONE;++i<iNoOfPts;) {
    ptArrNew[i].setPointId(i);
    ptArrNew[i].setCoordinate(ZERO,ptArr[i].getCoordinate(d)/pMax.getCoordinate(d));
  }
  vector<Cluster> vecCluster,vecCluster2;
  char caK[FIVE];
  sprintf(caK,"%d",int(NORM));
  Dim* dimArr;// = (Dim *)NULL;
  createClusters(ptArrNew,vecCluster,dimArr,caK,cMode);
  for(int i=MINUS_ONE;++i<iNoOfPts;) {
    vecCluster[int(ptArrNew[i].getCoordinate(d)*vecCluster.size())].add(ptArr[i],iaCluster);
  }
  vector<int> vDim(ONE);
  vDim[ZERO]=d;
  int iNumClusters = vecCluster.size();
  for(int i=MINUS_ONE; ++i < iNumClusters;) {

    // prune under-represented ClusterSets
//    if(vecCluster[i].getSize()<log(NO_OF_POINTS))
//      continue;
    ClusterSet *csNew = new ClusterSet();
    csNew->setDimensions(vDim);
    for(int j=MINUS_ONE; ++j < vecCluster[i].getSize();) {
      csNew->addPoint(vecCluster[i].getPtId(j));
    }
    sort(csNew->getPoints().begin(),csNew->getPoints().end());
    csNew->setProb(csNew->getPoints().size()/(1.0*NO_OF_POINTS));
    vecCluSet->push_back(csNew);
    ((*vecCluSet)[vecCluSet->size()-1]->getBase()).push_back((*vecCluSet)[vecCluSet->size()-1]);
  }
}

// nearest neighbour implementation
pair<double,int> findNearestCluster(int iPt,double* daMinDist,int* iaCluster,vector<Cluster> const & vecCluSrc,Point & ptSrc,list<Centre *> const& listCntr,ClusterMode cMode=KCLUSTER,double const dLpNorm=NORM) {
  double dMinDist = MAXDOUBLE,dDist,dNextMinDist = MAXDOUBLE;
  int iMinDistId = ZERO,iNumClusters = vecCluSrc.size(), iId = ZERO;
  double daDist[listCntr.size()];

  switch(cMode) {
    case KCLUSTER:
      for(list<Centre *>::const_iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
        dDist = distance(ptSrc,*(*lIter),dLpNorm);
        if(dMinDist > dDist) {
          dMinDist = dDist;
          iMinDistId = (*lIter)->getClusterId();
        }
      }
      break;

    case ECLUSTER:      // ECLUSTER uses techniques mentioned in the Elkan paper
      iMinDistId = iaCluster[ptSrc.getPointId()];   // c(x)
      if(iMinDistId == MINUS_ONE) {                 // find nearest centre
        for(list<Centre *>::const_iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
          if(dMinDist <= vecCluSrc[(*lIter)->getClusterId()].getVecDist(iMinDistId)/2) continue;                               // using Lemma 1
          dDist = distance(ptSrc,*(*lIter),dLpNorm);// d(x,c)
          daL[iPt*iNumClusters+(*lIter)->getClusterId()]=dDist;//ptSrc.setL((*lIter)->getClusterId(),dDist);
          if(dMinDist > dDist) {
            dMinDist = dDist;
            iMinDistId = (*lIter)->getClusterId();
          }
        }
        daU[iPt] = dMinDist; //ptSrc.setU(dMinDist);
        break;
      }
      if(daU[iPt] <= vecCluSrc[iMinDistId].getVecDist(vecCluSrc[iMinDistId].getNearestClusterId())/2) break;               // u(x) > d(c(x),c)/2
      for(list<Centre *>::const_iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
        if((*lIter)->getClusterId()==iMinDistId) continue; // c != c(x)
        if(ptSrc.tested() == true) {                // if r(x)
          dMinDist=distance(ptSrc,const_cast<Cluster&>(vecCluSrc[iMinDistId]).getCentre(),dLpNorm);
          ptSrc.setTested(false);                   // compute d(x,c(x))
          daMinDist[ptSrc.getPointId()]=dMinDist;
          daU[iPt] = dMinDist; //ptSrc.setU(dMinDist);
          daL[iPt*iNumClusters+iMinDistId]=dMinDist;//ptSrc.setL(iMinDistId,dMinDist);
        }
        if(daMinDist[ptSrc.getPointId()] > min(daL[iPt*iNumClusters+(*lIter)->getClusterId()],vecCluSrc[iMinDistId].getVecDist((*lIter)->getClusterId())/2)) {
          dDist = distance(ptSrc,const_cast<Cluster&>(vecCluSrc[(*lIter)->getClusterId()]).getCentre(),dLpNorm);
          daL[iPt*iNumClusters+(*lIter)->getClusterId()]=dDist;//ptSrc.setL((*lIter)->getClusterId(),dDist);
          if(dDist < dMinDist) {                    // if(d(x,c)<d(x,c(x))
            iMinDistId = (*lIter)->getClusterId();
            dMinDist = dDist;
            daU[iPt] = dMinDist; //ptSrc.setU(dMinDist);
            daMinDist[ptSrc.getPointId()]=dMinDist;
          }
        }
      }
      break;

    case JCLUSTER:      // JCLUSTER uses techniques mentioned in the Jain paper
      if(iaCluster[ptSrc.getPointId()]!=MINUS_ONE) {
        iMinDistId = iaCluster[ptSrc.getPointId()];
        daMp[iPt] -= (const_cast<Cluster &>(vecCluSrc[iMinDistId]).getMoved()+const_cast<Cluster &>(vecCluSrc[iMinDistId]).getMaxMoved());
        if(daMp[iPt] > 0) {
          dMinDist = distance(ptSrc,const_cast<Cluster &>(vecCluSrc[iMinDistId]).getCentre(),dLpNorm);
          break;
        }
      }
      for(list<Centre *>::const_iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
        if(dMinDist <= vecCluSrc[(*lIter)->getClusterId()].getVecDist(iMinDistId)/2) continue;
        daDist[++iId] = distance(ptSrc,*(*lIter),dLpNorm);
        if(dMinDist > daDist[iId]) {
          dMinDist = daDist[iId];
          iMinDistId = (*lIter)->getClusterId();
        }
      }
      iId = ZERO;
      for(list<Centre *>::const_iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
        ++iId;
        if((*lIter)->getClusterId() == iMinDistId) continue;
        if(dNextMinDist > daDist[iId]) {
          dNextMinDist = daDist[iId];
        }
      }
      daMp[iPt] = dNextMinDist-dMinDist;
      break;

    default:
      cerr << " Invalid mode for function to be called: " << cMode << endl;
      exit(0);
  }
  if(dMinDist >= MAXDOUBLE) {     // sanity check
    cerr << ptSrc.getPointId() << " HAS INFINITE DISTANCE TO NEAREST CLUSTER " << dMinDist << endl;
  }
  return pair<double,int>(dMinDist,iMinDistId);
}

// redundant function
pair<double,int> findNearestCentre(Point const& ptSrc,list<Centre *> listCntr,double const dLpNorm = NORM) {
  double dMinDist = MAXDOUBLE;
  double daMax[ptSrc.getCoordSize()];
  int iMinDistId = ZERO;
  Point pMax;
  getMax(ptSrc,listCntr,pMax);
  for(int d = MINUS_ONE; ++d < ptSrc.getCoordSize();) {
    daMax[d] = pow(max(ptSrc.getCoordinate(d),pMax.getCoordinate(d)-ptSrc.getCoordinate(d)),dLpNorm);
    dTotalCnt += ONE;
  }
  double dPlusMax = accumulate(daMax,daMax+ptSrc.getCoordSize(),0.0);

  // \forall d \in [1,NO_OF_DIMENSIONS]
  for(int d = MINUS_ONE; ++d < ptSrc.getCoordSize();) {
    double dKmax = MAXDOUBLE;
    dPlusMax -= daMax[d];

    // for each plausible centre j
    for(list<Centre *>::iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
      // update D^-(j), D^+_{max}(j)
      (*lIter)->setMinus((*lIter)->getMinus()+pow(abs(ptSrc.getCoordinate(d)-(*lIter)->getCoordinate(d)),dLpNorm));
      dTotalCnt += ONE;
      (*lIter)->setPlusMax(dPlusMax);

      // K_{max} = \min_{j}\{D^-(j)+D^+_{max}(j)\}
      double ddMax = (*lIter)->getMinus() + (*lIter)->getPlusMax();
      if(dKmax > ddMax) {
        dKmax = ddMax;
      }
    }
    if(listCntr.size() == ONE)
      continue;

    // for each plausible centre j
    for(list<Centre *>::iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {

      // if(D^-(j)+D^+_{min}(j) \geq K_{max}, vecCentre.erase(j)
      if((*lIter)->getMinus() + (*lIter)->getPlusMin() > (dKmax+0.05)) {
        if(lIter==listCntr.begin()) {
          listCntr.erase(lIter);
          lIter = listCntr.begin();
        } else {
          list<Centre *>::iterator lIter1=lIter--;
          listCntr.erase(lIter1);
        }
      }
    }
  }
  if(listCntr.size() >= ONE) {
    for(list<Centre *>::iterator lIter=listCntr.begin();lIter!=listCntr.end();++lIter) {
      if((*lIter)->getMinus() < dMinDist) {
        dMinDist = (*lIter)->getMinus();
        iMinDistId = (*lIter)->getClusterId();
      }
    }
  }
  return pair<double,int>(dMinDist,iMinDistId);
}

int assignPoint(int iPt,double* daMinDist,int* iaCluster,Point& ptSrc,vector<Cluster>& vecCluSrc,list<Centre *> const& listCntr,ClusterMode cMode=KCLUSTER,double const dLpNorm=NORM) {
  int iOldCluId = iaCluster[ptSrc.getPointId()];
  pair<double,int> pNearest = findNearestCluster(iPt,daMinDist,iaCluster,vecCluSrc,ptSrc,listCntr,cMode,dLpNorm);
  daMinDist[ptSrc.getPointId()]=pNearest.first;
  if(pNearest.second != iOldCluId) {
    if(iOldCluId != MINUS_ONE)
      vecCluSrc[iOldCluId].remove(ptSrc,iaCluster);
    vecCluSrc[pNearest.second].add(ptSrc,iaCluster);
    return ONE;
  }
  return ZERO;
}

void printClusters(vector<Cluster> & vecCluSrc) {
  for(unsigned int j = ZERO; j < vecCluSrc.size();j++) {
    vecCluSrc[j].print();
  }
}

int kClusterData(ifstream& fData,double* daMinDist,int* iaCluster,vector<Cluster>& vecCluSrc,Point* ptArr,ClusterMode cMode=KCLUSTER,char cDelim=DELIMITER,double const dLpNorm=NORM) {
  int iChanged = ZERO,iNumClusters=vecCluSrc.size();
  list<Centre *> listCntr;
  for(int j=MINUS_ONE;++j<iNumClusters;) {
    listCntr.push_back(&(vecCluSrc[j].getCentre()));
    vecCluSrc[j].updateVecDist(iNumClusters);
  }
  if(bRead) iPointCnt = ZERO;

  // \forall i \in [1,iNoOfPts]
  for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) {
    if(bRead) {
      Point pNew;
      if(!readPtHorizontally(fData,pNew,cDelim)) break;
      iChanged += assignPoint(iPt,daMinDist,iaCluster,pNew,vecCluSrc,listCntr,cMode,dLpNorm);
    } else {
      iChanged += assignPoint(iPt,daMinDist,iaCluster,ptArr[iPt],vecCluSrc,listCntr,cMode,dLpNorm);
    }
  }
  if(iChanged == ZERO) return ZERO;
  if(dLpNorm < ZERO) {     // convert to unit vector for spherical K-means
    for(int iCntr = MINUS_ONE;++iCntr < iNumClusters;) {
      vecCluSrc[iCntr].getCentre().normalise(*pMax,*pMin,dLpNorm);
    }
  }
  for(int j = MINUS_ONE; ++j < iNumClusters;) {
    vecCluSrc[j].updateOtherCentre();
    if(cMode == KCLUSTER) {
      vecCluSrc[j].setCentre(vecCluSrc[j].getOtherCentre());
    }
#ifdef DEBUG
    vecCluSrc[j].print();
#endif
  }
  if(cMode == ECLUSTER && iIter == ZERO) {
    for(int iPt = MINUS_ONE; ++iPt < NO_OF_POINTS;) {
      int iOffset = iPt*iNumClusters;
      for(int k = MINUS_ONE; ++k < iNumClusters;) {
        daL[iOffset+k]=max(daL[iOffset+k]-vecCluSrc[k].getMoved(),0.0);
//        ptArr[iPt].setL(k,max(ptArr[iPt].getL(k)-vecCluSrc[k].getMoved(),0.0));                                 // l(x,c)=max{l(x,c)-d(c,m(c)),0}
      }
//      ptArr[iPt].setU(ptArr[iPt].getU()+vecCluSrc[ptArr[iPt].getClusterId()].getMoved());
      daU[iPt] += vecCluSrc[iaCluster[ptArr[iPt].getPointId()]].getMoved();
      ptArr[iPt].setTested(true);
    }
    for(int j = MINUS_ONE; ++j < iNumClusters;)
      vecCluSrc[j].setCentre(vecCluSrc[j].getOtherCentre());
  }
  updateInterClusterDist(vecCluSrc,cMode,dLpNorm);
  return iChanged;
}

// redundant function
void getListOfPossibleCentres(double* daMinDist,vector<Cluster> & vecCluSrc,int j,int m,list<Centre *>& listCntr,Point* ptArr) {
  for(unsigned int l=ZERO; l < vecCluSrc.size();l++) {
    if(daMinDist[ptArr[vecCluSrc[j].getPtId(m)].getPointId()] >= (vecCluSrc[j].getVecDist(l)-vecCluSrc[j].getVecDist(j))/2) {
      listCntr.push_back(&(vecCluSrc[l].getOtherCentre()));
    }
  }
}

// redundant function
int oKClusterData(double* daMinDist,int* iaCluster,vector<Cluster>& vecCluSrc,Point* ptArr,ClusterMode cMode=KCLUSTER,double const dLpNorm=NORM) {
  int iChanged = ZERO,iNumClusters = vecCluSrc.size();

  // \forall j \in [1,k]
  if(iaCluster[ptArr[0].getPointId()]==MINUS_ONE)  {
    for(int j = MINUS_ONE; ++j < iNumClusters;) {
      vecCluSrc[j].updateVecDist(iNumClusters);
    }
    list<Centre  *> listCntr;
    for(int i=MINUS_ONE;++i<iNumClusters;)
      listCntr.push_back(&(vecCluSrc[i].getCentre()));

    // \forall i \in [1,ptArr.size()]
    for(int i = MINUS_ONE; ++i < NO_OF_POINTS;) 
      iChanged += assignPoint(i,daMinDist,iaCluster,ptArr[i],vecCluSrc,listCntr,cMode,dLpNorm);
  } else {
    for(int j = MINUS_ONE; ++j < iNumClusters;) {

      // \forall m \in [1,vecCluSrc[j].getMembers().size()]
      for(int m = MINUS_ONE; ++m < vecCluSrc[j].getSize();) {
        if(ptArr[vecCluSrc[j].getPtId(m)].tested()) continue;
        ptArr[vecCluSrc[j].getPtId(m)].setTested(true);
        list<Centre *> listCntr;
        getListOfPossibleCentres(daMinDist,vecCluSrc,j,m,listCntr,ptArr);
        if(listCntr.size()<TWO) continue;
        if(assignPoint(m,daMinDist,iaCluster,ptArr[vecCluSrc[j].getPtId(m)],vecCluSrc,listCntr,cMode,dLpNorm)) {
          ++iChanged;
          --m;
        }
      }
    }
  }
  if(iChanged == ZERO) return ZERO;
  for(int j = MINUS_ONE; ++j < iNumClusters;) {
    for(int m = MINUS_ONE; ++m < vecCluSrc[j].getSize();) {
      ptArr[vecCluSrc[j].getPtId(m)].setTested(false);
    }
    vecCluSrc[j].updateOtherCentre();
  }
  updateInterClusterDist(vecCluSrc,cMode);
  return iChanged;
}
