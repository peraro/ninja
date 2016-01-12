// This -*- C++ -*- header file contains typedefs and wrappers used in
// the Ninja library.


#ifndef NINJA_TYPES_HH
#define NINJA_TYPES_HH

#include <ninja/ninja_config.h>

//quadninja//#define QUADNINJA_TYPES_HH_INSIDE 1

#if defined(NINJA_QUADRUPLE) || defined(QUADNINJA_TYPES_HH_INSIDE)
# include <ninja/quadruple.hh>
#else
# include <cmath>
# include <complex>
#endif

#include <iostream>

#include <ninja/static_arrays.hh>
#include <ninja/zero_float.hh>

#if defined(NINJA_NO_EXCEPTIONS)
# define NINJA_THROW(exception) (std::terminate())
#else
# define NINJA_THROW(exception) throw exception
#endif

#define NINJA_REAL(x) (ninja::Real(x))

#if !defined(NINJA_QUADRUPLE) && !defined(QUADNINJA_TYPES_HH_INSIDE)
# define NINJA_REAL_DEF(x) (ninja::Real(x))
# define NINJA_COMPLEX_DEF(r,i) (ninja::Complex(r,i))
#else
# define NINJA_TOKEN_QPASTE0(x,y) x##y
# define NINJA_REAL_DEF(x) (ninja::Real(NINJA_TOKEN_QPASTE0(x,q)))
# define NINJA_COMPLEX_DEF(r,i) (ninja::Complex(NINJA_REAL_DEF(r)),ninja::Complex(NINJA_REAL_DEF(i)))
#endif


namespace ninja {

  // typedefs for Real and Complex floating-point types
#if !defined(NINJA_QUADRUPLE) && !defined(QUADNINJA_TYPES_HH_INSIDE)
  typedef double Real;
  typedef std::complex<Real> Complex;
  const Real INFRARED_EPS = 1.0e-09;
#else
  typedef ninja::Quadruple Real;
  typedef ninja::ComplexQuadruple Complex;
  const Real INFRARED_EPS = 1.0e+07*FLT128_EPSILON;
#endif

  // typedef for zero, real and complex masses
  typedef Real RealMasses;
  typedef Complex ComplexMasses;
  typedef ZeroFloat Massless;


  // Imaginary unit
  const Complex I(Real(0.),Real(1.));


  // Put complex sqrt in ninja-namespace
  inline Complex sqrt(const Complex & z)
  {
    return std::sqrt(z);
  }
  // Put real in ninja-namespace
  inline Real real(const Complex & z)
  {
    return std::real(z);
  }
  // Put imag in ninja-namespace
  inline Real imag(const Complex & z)
  {
    return std::imag(z);
  }
  // Put conj in ninja-namespace
  inline Complex conj(const Complex & z)
  {
    return std::conj(z);
  }
  // Put abs in ninja-namespace
  inline Real abs(const Complex & z)
  {
    return std::abs(z);
  }

  // The taxicab norm (or Manhattan norm) in the complex plane
  // 
  //    ||z|| = |real(z)| + |imag(z)|
  //
  // Its computation should be faster than std::abs(z)
  inline Real taxicab_norm (const Complex & z)
  {
    return std::abs(real(z))+std::abs(imag(z));
  }

  // overrides taxicab_norm for real types
  inline Real taxicab_norm(const Real & x)
  {
    return std::abs(x);
  }


  // const pointer type
  template<typename X> struct const_pointer
  {
    typedef const X* type;
  };

  // Specialize ninja::const_pointer<ZeroFloat>
  template<> struct const_pointer<ZeroFloat> {
    typedef ZeroFloatArray type;
  };


  // Convert Zero-Floats to reals by asking for real part
  inline Real real (ZeroFloat)
  {
    return Real();
  }

  // Converting ZeroFloats to Complex
  inline Complex toCmplx (const Complex & z)
  {
    return z;
  }
  inline Complex toCmplx (const ZeroFloat &)
  {
    return Complex();
  }


  // some constants

  const Real ZERO = Real(0.0);
  const Real INV8 = Real(0.125);
  const Real INV4 = Real(0.25);
  const Real HALF = Real(0.5);
  const Real ONE = Real(1.);
  const Real ONEDOTFIVE = Real(1.5);
  const Real TWO = Real(2.);
  const Real THREE = Real(3.);
  const Real FOUR = Real(4.);
  const Real FIVE = Real(5.);
  const Real SIX = Real(6.);
  const Real EIGHT = Real(8.);
  const Real TWELVE = Real(12.);
  const Real SIXTEEN = Real(16.);
  const Real SQRT2 = std::sqrt(TWO);
  const Real SQRT3 = std::sqrt(THREE);

  const Real INVSQRT2 = HALF*SQRT2;
  const Real INVSQRT3 = 1./std::sqrt(THREE);
  const Real INVSQRT6 = 1./std::sqrt(Real(6.));

#if !defined(NINJA_QUADRUPLE) && !defined(QUADNINJA_TYPES_HH_INSIDE)
  const ninja::Real PI = M_PI;
#else
  const ninja::Real PI = M_PIq;
#endif

} // namespace ninja

#undef QUADNINJA_TYPES_HH_INSIDE

#endif // NINJA_TYPES_HH
