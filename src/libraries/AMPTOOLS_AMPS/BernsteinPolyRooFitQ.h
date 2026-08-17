#ifndef BERNSTEINPOLYROOFITQ
#define BERNSTEINPOLYROOFITQ

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooAbsReal;
class RooArgList;

class BernsteinPolyRooFitQ : public RooAbsPdf {
public:
  BernsteinPolyRooFitQ() {} ;
  BernsteinPolyRooFitQ(const char *name, const char *title,
                RooAbsReal& _x,
                RooAbsReal& _massA,
                RooAbsReal& _massB,
                RooAbsReal& _qmin,
                RooAbsReal& _qmax,
                const RooArgList& _coefList);

  BernsteinPolyRooFitQ(const BernsteinPolyRooFitQ& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new BernsteinPolyRooFitQ(*this, newname); }
  inline virtual ~BernsteinPolyRooFitQ() { }

  Double_t breakupMomentum(double mass0, double mass1, double mass2) const;

protected:

  RooRealProxy x;
  RooRealProxy massA;
  RooRealProxy massB;
  RooRealProxy qmin;
  RooRealProxy qmax;
  RooListProxy coefList;

  Double_t evaluate() const;

private:

};

#endif