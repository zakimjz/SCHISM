class Dim {

  // private members
  private:
  vector<double> vecCoord;
  double dMax;
  double dMin;
  double dSum;

  public:

  // constructors
  Dim(int iNoOfPts = NO_OF_DIMENSIONS) {  // WHAT'S THIS?
    vecCoord.resize(iNoOfPts);
  }
  Dim(Dim const & dSrc,int iNoOfPts = NO_OF_DIMENSIONS) {
    dMax = dSrc.getMax();
    dMin = dSrc.getMin();
    dSum = dSrc.getSum();
    vecCoord.resize(iNoOfPts);
    vecCoord = dSrc.getCoord();
  }

  // public methods
  double getCoordinate(int iPtId) const {return vecCoord[iPtId]; }
  vector<double> const & getCoord() const { return vecCoord; }
  int getSize() const { return vecCoord.size(); }
  double getMax() const { return dMax; }
  double getMin() const { return dMin; }
  double const getSum() const { return dSum; }

  void setCoordinate(int iPtId,double dVal) { vecCoord[iPtId]=dVal; }
  void setSize(int iSize) { vecCoord.resize(iSize); }
  void setMax(double dMaxSrc) { dMax = dMaxSrc; }
  void setMin(double dMinSrc) { dMin = dMinSrc; }
  void setSum(double dSumSrc) { dSum = dSumSrc; }
  void updateMin() { dMin = *(min_element(vecCoord.begin(),vecCoord.end())); }
  void updateMax() { dMax = *(max_element(vecCoord.begin(),vecCoord.end())); }
  void updateSum() { dSum = accumulate(vecCoord.begin(),vecCoord.end(),0.0); }
  void normalise(int iId) { vecCoord[iId] = (vecCoord[iId]-dMin)/(dMax-dMin); }
  void denormalise(int iId) { vecCoord[iId] = vecCoord[iId]*(dMax-dMin)+dMin; }
  void normalise() {
    for(unsigned int i=ZERO;i < vecCoord.size();)normalise(i++);
  }
  void denormalise() {
    for(unsigned int i=ZERO;i < vecCoord.size();)denormalise(i++);
  }
  bool operator>(Dim const& DimSrc) const { return dSum > DimSrc.getSum(); }
  bool operator<(Dim const& DimSrc) const { return dSum < DimSrc.getSum(); }
};
