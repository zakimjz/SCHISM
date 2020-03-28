class Cluster {

  // private members
  Centre pCurrent;
  Centre pNew;
  Centre pNext;
  int iClusterId;
  int iNearestClusterId;    // used for JCLUSTER - spheres of G assignment
  double dMoved;            // used for JCLUSTER - max movement effect
  double dMaxMoved;         // used for JCLUSTER - max movement effect
  vector<int> vecIds;     // size = # points assigned to cluster
  vector<double> vecDist;   // size = # clusters

  public:

  // constructors
  Cluster(Centre pCntrSrc) {
    iClusterId = iClusterCnt++;
    iNearestClusterId = ZERO;
    if(iClusterId==ZERO)iNearestClusterId = ONE;
    pNext = pCurrent = pCntrSrc;
    pNext.setClusterId(iClusterId);
    pCurrent.setClusterId(iClusterId);
  }

  // public methods

  Centre & getCentre() { return pCurrent; }
  Centre & getOtherCentre() { return pNext; }
  int getPtId(int iIndex) { return vecIds[iIndex];}
  int getSize() const { return vecIds.size();}
  double getMoved() const { return dMoved; }
  double getMaxMoved() const { return dMaxMoved; }
  double getVecDist(int iId) const { return vecDist[iId]; }
  int getNearestClusterId() const { return iNearestClusterId; }

  void setCentre(Centre& pCntrSrc) { pCurrent = pCntrSrc;}
  void setMoved(double dMovedSrc) { dMoved = dMovedSrc; }
  void setMaxMoved(double dMaxMovedSrc) { dMaxMoved = dMaxMovedSrc; }
  void setVecDist(int iId,double dDist) { vecDist[iId]=dDist; }
  void setNearestClusterId(int iIdSrc) { iNearestClusterId = iIdSrc; }
  void updateOtherCentre();
  void updateVecDist(int iSize) { vecDist.resize(iSize); }

  double minVecDist() { 
    double dMin1 = *min(vecDist.begin(),vecDist.begin()+iClusterId);
    double dMin2 = *min(vecDist.begin()+iClusterId+ONE,vecDist.end());
    return (dMin1 < dMin2) ? dMin1 : dMin2;
  }
  void add(Point& pSrc,int *iaCluster);
  void remove(Point& pSrc,int *iaCluster);
  void print() const;
};

void Cluster::add(Point& pSrc,int* iaCluster) {
  iaCluster[pSrc.getPointId()] = iClusterId;
  vecIds.push_back(pSrc.getPointId());
  for(int i = MINUS_ONE; ++i < pSrc.getCoordSize();) {
    pNew.setCoordinate(i,pNew.getCoordinate(i)+pSrc.getCoordinate(i));
  }
}

void Cluster::remove(Point& pSrc,int* iaCluster) {
  for(int i = MINUS_ONE; ++i < pSrc.getCoordSize();) {
    pNew.setCoordinate(i,pNew.getCoordinate(i)-pSrc.getCoordinate(i));
  }
  iaCluster[pSrc.getPointId()] = MINUS_ONE;
  vector<int>::iterator iIter=vecIds.begin();
  if((iIter=find(vecIds.begin(),vecIds.end(),pSrc.getPointId()))!=vecIds.end()) {
    vecIds.erase(iIter);
  } else {
    cerr << " Point " << pSrc.getPointId() << " NOT FOUND " << endl;
  }
}

void Cluster::print() const {
  cout << " Cluster ID = " << iClusterId << " PointIds: ";
#ifdef DEBUG
  for(unsigned int i = ZERO; i < vecIds.size();++i) {
    cout << vecIds[i] << ' ';
  }
  cout << endl;
#endif
  pNext.print();
}

void Cluster::updateOtherCentre() {
  if(vecIds.size()==ZERO) return;
  for(int i = MINUS_ONE; ++i < pNew.getCoordSize();) {
    pNext.setCoordinate(i,pNew.getCoordinate(i)/vecIds.size());
  }
}
