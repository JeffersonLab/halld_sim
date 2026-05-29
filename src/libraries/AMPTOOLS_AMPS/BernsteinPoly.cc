#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <cmath>

#include "TLorentzVector.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/BernsteinPoly.h"

// ---------------------------------------------------------------------------
// Helper: integer binomial coefficient C(n, k)
// ---------------------------------------------------------------------------
static int binomialCoeff( int n, int k )
{
    if( k < 0 || k > n ) return 0;
    if( k == 0 || k == n ) return 1;
    int result = 1;
    for( int i = 0; i < k; ++i ){
        result *= ( n - i );
        result /= ( i + 1 );
    }
    return result;
}

// ---------------------------------------------------------------------------
// Constructor
//
// args[0]           : xmin
// args[1]           : xmax
// args[2]           : degree
// args[3]           : daughter 1 index string  (e.g. "2")
// args[4]           : daughter 2 index string  (e.g. "3")
// args[5..5+degree] : Bernstein coefficients [c0] ... [cN]
//
// Mirrors the BreitWigner pattern: daughters identified by pKin index strings,
// coefficients registered as AmpParameters.
// ---------------------------------------------------------------------------
BernsteinPoly::BernsteinPoly( const vector< string >& args )
  : UserAmplitude< BernsteinPoly >( args )
{
    // minimum: xmin xmax degree d1 d2 c0  --> 6 args for degree=0
    assert( args.size() >= 6 );

    m_xmin   = atof( args[0].c_str() );
    m_xmax   = atof( args[1].c_str() );
    m_degree = atoi( args[2].c_str() );

    m_daughters = pair< string, string >( args[3], args[4] );

    int nCoeffs = m_degree + 1;
    assert( (int)args.size() == nCoeffs + 5 );

    cout << "BernsteinPoly:"
         << "  xmin="    << m_xmin
         << "  xmax="    << m_xmax
         << "  degree="  << m_degree
         << "  d1="      << m_daughters.first
         << "  d2="      << m_daughters.second
         << "  nCoeffs=" << nCoeffs << endl;

    m_coeffs.reserve( nCoeffs );
    for( int i = 0; i < nCoeffs; ++i ){
        m_coeffs.push_back( AmpParameter( args[5 + i] ) );
        registerParameter( m_coeffs.back() );
        cout << "  c" << i << " initialised to " << (double)m_coeffs[i] << endl;
    }
}

// ---------------------------------------------------------------------------
// calcAmplitude
//
// Builds the invariant mass of the two daughters from pKin indices,
// exactly as BreitWigner does, then evaluates:
//
//   t    = ( m - xmin ) / ( xmax - xmin )   in [0, 1]
//   B(t) = sum_{i=0}^{N} c_i * C(N,i) * t^i * (1-t)^(N-i)
//
// Returns sqrt( |B(t)| ) so that |amplitude|^2 = B(t).
// Events outside [xmin, xmax] return 0.
// ---------------------------------------------------------------------------
complex< GDouble >
BernsteinPoly::calcAmplitude( GDouble** pKin ) const
{
    TLorentzVector P1, P2, Ptemp;

    // accumulate daughter 1 (supports multi-digit index strings like BreitWigner)
    for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
        string num; num += m_daughters.first[i];
        int index = atoi( num.c_str() );
        Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                          pKin[index][3], pKin[index][0] );
        P1 += Ptemp;
    }

    // accumulate daughter 2
    for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
        string num; num += m_daughters.second[i];
        int index = atoi( num.c_str() );
        Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                          pKin[index][3], pKin[index][0] );
        P2 += Ptemp;
    }

    GDouble mass = ( P1 + P2 ).M();

    // return zero outside the defined range
    if( mass < m_xmin || mass > m_xmax )
        return complex< GDouble >( 0.0, 0.0 );

    GDouble t     = ( mass - m_xmin ) / ( m_xmax - m_xmin );
    GDouble one_t = 1.0 - t;

    GDouble bern = 0.0;
    for( int i = 0; i <= m_degree; ++i ){
        int     binom = binomialCoeff( m_degree, i );
        GDouble term  = (double)m_coeffs[i]
                        * binom
                        * pow( t,     (GDouble)i )
                        * pow( one_t, (GDouble)( m_degree - i ) );
        bern += term;
    }

    return complex< GDouble >( sqrt( fabs( bern ) ), 0.0 );
}

// ---------------------------------------------------------------------------
// updatePar — called by the framework when a parameter value changes.
// No expensive cached quantities here, so nothing to do.
// ---------------------------------------------------------------------------
void
BernsteinPoly::updatePar( const AmpParameter& par )
{
    // nothing to recompute
}

// ---------------------------------------------------------------------------
// GPU stub — full implementation will live in GPUBernsteinPoly.cu
// ---------------------------------------------------------------------------
#ifdef GPU_ACCELERATION
void
BernsteinPoly::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const
{
    // Encode daughter index strings as integers (mirrors BreitWigner convention)
    int daught1 = atoi( m_daughters.first.c_str() );
    int daught2 = atoi( m_daughters.second.c_str() );

    // Flatten coefficients into a plain array for the GPU kernel
    vector< GDouble > coeffArr( m_degree + 1 );
    for( int i = 0; i <= m_degree; ++i )
        coeffArr[i] = (double)m_coeffs[i];

    GPUBernsteinPoly_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
                           m_xmin, m_xmax, m_degree,
                           daught1, daught2,
                           coeffArr.data() );
}
#endif // GPU_ACCELERATION

