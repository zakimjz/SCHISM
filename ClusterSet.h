// subspace clustering-related class
class ClusterSet {

  // private members
  vector<int> vDim;
  vector<ClusterSet*> vBase;
  vector<ClusterSet*> vChildren;
//  ClusterSet* csParent1;
//  ClusterSet* csParent2;
  vector<int> vPoints;
  double dProb;

  public:

  // constructors
  ClusterSet() { vDim.resize(ONE);}
  ClusterSet(ClusterSet & csSrc) { 
    setDimensions(csSrc.getDimensions());
    setBase(csSrc.getBase());
    setChildren(csSrc.getChildren());
    setPoints(csSrc.getPoints());
    setProb(csSrc.getProb());
  }
 
  // methods
  vector<int>& getDimensions() { return vDim; }
  vector<ClusterSet*>& getBase() { return vBase; }
  vector<ClusterSet*>& getChildren() { return vChildren; }
//  ClusterSet& getParent1() { return *csParent1; }
//  ClusterSet& getParent2() { return *csParent2; }
  vector<int>& getPoints() { return vPoints; }
  double getProb() const { return dProb; }

  void setDimensions(vector<int> & vDimen) { 
    vDim.resize(vDimen.size());
    vDim = vDimen;
  }
  void setBase(vector<ClusterSet*> & vBaseSrc) { 
    vBase.resize(vBaseSrc.size());
    vBase = vBaseSrc;
  }
  void setChildren(vector<ClusterSet *> & vChildren) { 
    vChildren.resize(vChildren.size());
    vChildren = vChildren;
  }
//  void setParent1(ClusterSet* csParent) { csParent1 = csParent; }
//  void setParent2(ClusterSet* csParent) { csParent2 = csParent; }
  void setPoints(vector<int> const& vPts) { 
    vPoints.resize(vPts.size());
    vPoints = vPts;
  }
  void setProb(double dSrc) { dProb=dSrc; }

  double updateProb(vector<ClusterSet *>& vecCluSet,int iNoDimXamined) {
    double dProb = ONE;
    double dSum = ZERO;
    int k = ZERO;

    // \for_all dimensions 
    for(int i=MINUS_ONE;++i < iNoDimXamined;) {

      // in which they are they are identical
      if(vDim[k]==i) {
        dSum = vBase[k++]->getProb();
      } else {
        dSum = ONE;

        // \for_all ClusterSets corresponding to the dimension
//        for(int j=i*int(NORM)-ONE;++j < (i+1)*int(NORM);) {
//          dSum -= pow(vecCluSet[j]->getProb(),1.0*vPoints.size());
//        }
      }
      dProb *= dSum;
    }
    dProb *= (iNoDimXamined*1.0/vPoints.size());
    return dProb;
  }
  void addChild(ClusterSet* csSrc) { vChildren.push_back(csSrc); }
  void addPoint(int iPt) { vPoints.push_back(iPt); }
  void print() const {
    cout << " Count: " << vPoints.size() << " dProb: " << dProb << " Dim: " ;
    for(unsigned int i=ZERO;i < vDim.size();) cout << vDim[i++] << ' ';
//    cout << endl << " Points: ";
//    for(int i=MINUS_ONE;++i < vPoints.size();) cout << vPoints[i] << ' ';
    cout << endl;
  }
};
