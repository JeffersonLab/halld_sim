#if !defined(BERNSTEINPOLY)
#define BERNSTEINPOLY

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPUBernsteinPoly_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                             GDouble xmin, GDouble xmax, int degree,
                             int daught1, int daught2,
                             const GDouble* coeffs );
#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

//
// BernsteinPoly: a Bernstein polynomial background amplitude evaluated on
// the invariant mass of two daughter particles specified by index.
//
// Config file usage:
//   define bern_xmin  0.6
//   define bern_xmax  2.5
//   define bern_deg   3
//
//   amplitude <reaction>::<sum>::<ampName> BernsteinPoly bern_xmin bern_xmax bern_deg <d1> <d2> [c0] [c1] ... [cN]
//
// Positional arguments:
//   args[0]             : xmin   (double)
//   args[1]             : xmax   (double)
//   args[2]             : degree (int)
//   args[3]             : daughter 1 index string (e.g. "2")
//   args[4]             : daughter 2 index string (e.g. "3")
//   args[5..5+degree]   : Bernstein coefficients [c0] ... [cN]  (AmpParameters)
//
// Total args = 5 + (degree + 1) = degree + 6
//
// calcAmplitude returns sqrt( |B(m)| ) so |amplitude|^2 = B(m),
// matching the AmpTools convention used by BreitWigner.
//
class BernsteinPoly : public UserAmplitude< BernsteinPoly >
{

public:

    BernsteinPoly() : UserAmplitude< BernsteinPoly >() {}
    BernsteinPoly( const vector< string >& args );

    ~BernsteinPoly() {}

    string name() const { return "BernsteinPoly"; }

    complex< GDouble > calcAmplitude( GDouble** pKin ) const;

    void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

    void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
    bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:

    GDouble m_xmin;
    GDouble m_xmax;
    int     m_degree;

    pair< string, string > m_daughters;

    // Bernstein coefficients c0 ... cN registered as free AmpParameters
    vector< AmpParameter > m_coeffs;
};

#endif // BERNSTEINPOLY

