#include <vector>
#include <iostream>
#include "convertData.h"

using namespace std;

class Point {

  // private members
  private:
  int iPointId;
  vector<double> vecCoord;
  bool bTested;

  public:

  // constructors
  Point(unsigned int iNoOfDims = NO_OF_DIMENSIONS) {
    iPointId=iPointCnt++;
    vecCoord.resize(iNoOfDims);
    bTested = false;
  }
  Point(Point const& pSrc,int iNoOfDims = NO_OF_DIMENSIONS) {
    iPointId = pSrc.getPointId();
    for(int d=MINUS_ONE;++d < iNoOfDims;)
      vecCoord.push_back(pSrc.getCoordinate(d));
    bTested = pSrc.tested();
  }

  // public methods
  double getCoordinate(int iDimension) const {return vecCoord[iDimension]; }
  int getPointId() const { return iPointId; }
  vector<double>& getCoord() { return vecCoord; }
  bool tested() const { return bTested; }

  void setCoordinate(int iDimension,double dVal) { vecCoord[iDimension]=dVal; }
  void setPointId(int iPtIdSrc) { iPointId = iPtIdSrc; }
  void setTested(bool bFlag) { bTested = bFlag; }

  int getCoordSize() const { return vecCoord.size(); }
  void normalise(Point const& pMax,Point const& pMin,double dLpNorm=NORM) {
    if(dLpNorm < ZERO) {
      double dMagnitude = ZERO;
      for(int d = MINUS_ONE;++d < getCoordSize();) dMagnitude += (vecCoord[d]*vecCoord[d]);
      dMagnitude = sqrt(dMagnitude);
      for(int d = MINUS_ONE;++d < getCoordSize();) vecCoord[d] /= dMagnitude;
      return;
    }
    for(unsigned int d=ZERO; d < vecCoord.size(); ++d) {
      vecCoord[d] = (vecCoord[d]-pMin.getCoordinate(d))/(pMax.getCoordinate(d)-pMin.getCoordinate(d));
    }
  }
  void denormalise(Point const& pMax,Point const& pMin) {
    for(unsigned int d=ZERO; d < vecCoord.size(); ++d) {
      vecCoord[d] = vecCoord[d]*(pMax.getCoordinate(d)-pMin.getCoordinate(d))+pMin.getCoordinate(d);
    }
  }
  bool operator==(Point const& pSrc) const { return pSrc.getPointId()==iPointId;}
  void print() const;
};

void Point::print() const {
  cout << " PointId = " << iPointId;
  cout << " Coordinates = ";
  for(unsigned int d = ZERO;d < vecCoord.size(); ++d) {
    cout << vecCoord[d] << ' ';
  }
  cout << endl;
}
