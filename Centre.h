class Centre: public Point {

  // private members
  int iClusterId;
  double dPlusMax;   // used for the Jain paper techniques and vertical spherical K-means
  double dPlusMin;
  double dMinus;

  public:

  // constructors
  Centre(Point & pSrc) : iClusterId(MINUS_ONE) {
    setPointId(pSrc.getPointId());
    getCoord().resize(pSrc.getCoordSize());
    dPlusMax = ZERO;
    for(int d=MINUS_ONE;++d < pSrc.getCoordSize();) {
      setCoordinate(d,pSrc.getCoordinate(d));
      dPlusMax += pSrc.getCoordinate(d);   // T(v^+)
    }
    dPlusMin = ZERO;
  }

  Centre(int iDim=NO_OF_DIMENSIONS) : iClusterId(MINUS_ONE) {
    setPointId(iPointCnt++);
    getCoord().resize(iDim);
    dPlusMin = ZERO;
  }

  // public methods
  double getPlusMax() const { return dPlusMax; }
  double getPlusMin() const { return dPlusMin; }
  double getMinus() const { return dMinus; }
  int getClusterId() { return iClusterId; }
  void setPlusMax(double dPlusMaxSrc) { dPlusMax = dPlusMaxSrc; }
  void setPlusMin(double dPlusMinSrc) { dPlusMax = dPlusMinSrc; }
  void setMinus(double dMinusSrc) { dMinus = dMinusSrc; }
  void setClusterId(int iSrc) { iClusterId = iSrc; }
};
