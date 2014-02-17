// Declarations for the interface of the ninja library to the Fortran
// code generated by GoSam.


#ifndef NINJA_GOSAM_INTERFACE_HH
#define NINJA_GOSAM_INTERFACE_HH

#include <ninja/ninja.hh>

namespace ninja {
  namespace gosam_interface {

	// Function which evaluates the numerator
	typedef void (*Numerator)(const int & ncut, const ComplexMomentum & Q,
							  const Complex & mu2, Complex & result);

	// Routine which returns the t-expansion
	typedef void (*Numerator_t3)(const int & ncut,
								 const ComplexMomentum & a,
								 const ComplexMomentum & b,
								 const ComplexMomentum & c,
								 const Complex & param,
								 const int & ndeg,
								 Complex * coeffs);

	// Routine which returns the t-expansion
	typedef void (*Numerator_t2)(const int & ncut,
								 const ComplexMomentum & a0,
								 const ComplexMomentum & a1,
								 const ComplexMomentum & b,
								 const ComplexMomentum & c,
								 const Complex * param,
								 const int & ndeg,
								 Complex * coeffs);

	// Routine which returns the mu-expansion
	typedef void (*Numerator_d)(const int & ncut,
								const ComplexMomentum a[],
								Complex * coeffs);

  } // namespace gosam_interface
} // namespace ninja

extern "C"{

  // These are the routines to be called by the Fortran code in GoSam

  // Selects a subset of propagators

  // Recution with real masses
  void ninjago_diag_rm(ninja::gosam_interface::Numerator numerator,
					   ninja::gosam_interface::Numerator_t3 numerator_t3,
					   ninja::gosam_interface::Numerator_t2 numerator_t2,
					   ninja::gosam_interface::Numerator_d numerator_d,
					   const int & nprops_group,
					   const int & nprops, const int & rk,
					   const int * indices,
					   const ninja::RealMomentum * Vi,
					   const ninja::Real * msq,
					   ninja::Real * g_mat_data,
					   const ninja::Real & scale2,
					   const int & istop,
					   ninja::Complex tot[3],
					   ninja::Complex & totr,
					   int & return_status);

  // Reduction with complex masses
  void ninjago_diag_cm(ninja::gosam_interface::Numerator numerator,
					   ninja::gosam_interface::Numerator_t3 numerator_t3,
					   ninja::gosam_interface::Numerator_t2 numerator_t2,
					   ninja::gosam_interface::Numerator_d numerator_d,
					   const int & nprops_group,
					   const int & nprops, const int & rk,
					   const int * indices,
					   const ninja::RealMomentum * Vi,
					   const ninja::Complex * msq,
					   ninja::Real * g_mat_data,
					   const ninja::Real & scale2,
					   const int & istop,
					   ninja::Complex tot[3],
					   ninja::Complex & totr,
					   int & return_status);

  // Reduction with massless loop-propagators
  void ninjago_diag_nm(ninja::gosam_interface::Numerator numerator,
					   ninja::gosam_interface::Numerator_t3 numerator_t3,
					   ninja::gosam_interface::Numerator_t2 numerator_t2,
					   ninja::gosam_interface::Numerator_d numerator_d,
					   const int & nprops_group,
					   const int & nprops, const int & rk,
					   const int * indices,
					   const ninja::RealMomentum * Vi,
					   ninja::Real * g_mat_data,
					   const ninja::Real & scale2,
					   const int & istop,
					   ninja::Complex tot[3],
					   ninja::Complex & totr,
					   int & return_status);


  // Clear the cache of the master integrals
  void ninjago_clear_integral_cache();

  // Set the verbosity
  void ninjago_set_verbosity(int val);

  // Set the the test flag
  void ninjago_set_test(int val);

  // Set the tolerance for the relative errors in the tests
  void ninjago_set_test_tolerance(const ninja::Real & val);

  // Set number of digits in the output
  void ninjago_set_output_precision(int val);

  // Select the integral library
  void ninjago_set_integral_library(int flag);

} // extern "C"

#endif // NINJA_GOSAM_INTERFACE_HH
