#include "Riostream.h"
#include "BernsteinPolyRooFitQ.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include <math.h>
#include "TMath.h"
#include "TIterator.h"

BernsteinPolyRooFitQ::BernsteinPolyRooFitQ(const char* name, const char* title,
                              RooAbsReal& _x,
                              RooAbsReal& _massA,
                              RooAbsReal& _massB,
                              RooAbsReal& _qmin,
                              RooAbsReal& _qmax,
                              const RooArgList& _coefList) :
  RooAbsPdf(name, title),
  x("x", "Dependent (mass)", this, _x),
  massA("massA", "Daughter mass A", this, _massA),
  massB("massB", "Daughter mass B", this, _massB),
  qmin("qmin", "Minimum breakup momentum", this, _qmin),
  qmax("qmax", "Maximum breakup momentum", this, _qmax),
  coefList("coefList", "List of coefficients", this)
{
  TIterator* coefIter = _coefList.createIterator();
  RooAbsArg* coef;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << "BernsteinPolyRooFitQ::ctor(" << GetName()
                << ") ERROR: coefficient " << coef->GetName()
                << " is not of type RooAbsReal" << std::endl;
      assert(0);
    }
    coefList.add(*coef);
  }
  delete coefIter;
}

BernsteinPolyRooFitQ::BernsteinPolyRooFitQ(const BernsteinPolyRooFitQ& other, const char* name) :
  RooAbsPdf(other, name),
  x("x", this, other.x),
  massA("massA", this, other.massA),
  massB("massB", this, other.massB),
  qmin("qmin", this, other.qmin),
  qmax("qmax", this, other.qmax),
  coefList("coefList", this, other.coefList)
{
}

Double_t BernsteinPolyRooFitQ::breakupMomentum(double mass0, double mass1, double mass2) const
{
  double q;
  q = sqrt( fabs(   mass0*mass0*mass0*mass0 +
                     mass1*mass1*mass1*mass1 +
                     mass2*mass2*mass2*mass2 -
                     2.0*mass0*mass0*mass1*mass1 -
                     2.0*mass0*mass0*mass2*mass2 -
                     2.0*mass1*mass1*mass2*mass2  ) ) / (2.0 * mass0);
  return q;
}

Double_t BernsteinPolyRooFitQ::evaluate() const
{
  Int_t degree = coefList.getSize() - 1;
  if (degree < 0) return 0.0;

  // Compute breakup momentum for this value of x (mass)
  Double_t q = breakupMomentum(x, massA, massB);

  // Rescale to [0, 1] using qmin/qmax
  Double_t t = (q - qmin) / (qmax - qmin);
  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;

  Double_t result = 0.0;

  // Bernstein basis: sum_i coef_i * C(n,i) * t^i * (1-t)^(n-i)
  for (Int_t i = 0; i <= degree; i++) {
    Double_t coef = ((RooAbsReal&)coefList[i]).getVal();
    Double_t binom = TMath::Binomial(degree, i);
    result += coef * binom * TMath::Power(t, i) * TMath::Power(1 - t, degree - i);
  }

  return result;
}