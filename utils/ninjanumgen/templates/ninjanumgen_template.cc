// This -*- C++ -*- file has been generated by the Python package
// NinjaNumGen, which is distributed with the Ninja library.

#ninja_includes

#include <ninja/num_defs.hh>
#include <ninja/ninjanumgen.hh>

#ninja_cpp_definition NINJA_NUM
#ninja_cpp_definition NINJA_RANK_MINUS_N
#ninja_cpp_definition NINJA_NUM_NAMESPACE

#if NINJA_NUM_NAMESPACE
namespace NINJA_NUM_NAMESPACE {
#endif

  ninja::Complex NINJA_NUM::evaluate(const ninja::ComplexMomentum & ninjaQ,
									 const ninja::Complex & ninjaMu2,
									 int ninja_cut,
									 const ninja::PartitionInt ninja_cutidx[])
  {
	ninja::Complex ninja_result;
#ninja_inline
#ninja_inline_evaluate
    (void)(ninjaQ);
	(void)(ninjaMu2);
	(void)(ninja_cut);
	(void)(ninja_cutidx);
#ninja_numerator_evaluate
	return ninja_result;
  }

  void NINJA_NUM::muExpansion(const ninja::ComplexMomentum ninjaV[],
							  const ninja::PartitionInt ninja_cutidx[],
							  ninja::Complex ninjaC[])
  {
#ninja_mu2expansion_indices
	ninja::ComplexMomentum ninjaA0 = ninjaV[0];
#if (NINJA_RANK_MINUS_N > 0)
	ninja::ComplexMomentum ninjaA1 = ninjaV[1];
	(void)(ninjaA1);
#endif
#ninja_inline
#ninja_inline_mu2expansion
	(void)(ninjaA0);
	(void)(ninja_cutidx);
	(void)(ninjaC);
#ninja_numerator_mu2expansion
  }

  void NINJA_NUM::t3Expansion(const ninja::ComplexMomentum & ninjaA,
							  const ninja::ComplexMomentum & ninjaE3,
							  const ninja::ComplexMomentum & ninjaE4,
							  const ninja::Complex & ninjaP,
							  int ninja_mindeg, int ninja_cut,
							  const ninja::PartitionInt ninja_cutidx[],
							  ninja::Complex ninjaC[])
  {
#ninja_t3expansion_indices
#ninja_inline
#ninja_inline_t3expansion
	(void)(ninjaA);
	(void)(ninjaE3);
	(void)(ninjaE4);
	(void)(ninjaP);
	(void)(ninja_mindeg);
	(void)(ninja_cut);
	(void)(ninja_cutidx);
	(void)(ninjaC);
	{
#ninja_numerator_t3expansion1
	}
	if (ninja_mindeg > 1+NINJA_RANK_MINUS_N) {
#ninja_numerator_t3expansion2
	}
  }

  void NINJA_NUM::t2Expansion(const ninja::ComplexMomentum & ninjaA0,
							  const ninja::ComplexMomentum & ninjaA1,
							  const ninja::ComplexMomentum & ninjaE3,
							  const ninja::ComplexMomentum & ninjaE4,
							  const ninja::Complex ninjaP[],
							  int ninja_mindeg, int ninja_cut,
							  const ninja::PartitionInt ninja_cutidx[],
							  ninja::Complex ninjaC[])
  {
#ninja_t2expansion_indices
#ninja_inline
#ninja_inline_t2expansion
	const ninja::Complex ninjaP0 = ninjaP[0];
	const ninja::Complex ninjaP1 = ninjaP[1];
	const ninja::Complex ninjaP2 = ninjaP[2];
	(void)(ninjaA0);
	(void)(ninjaA1);
	(void)(ninjaE3);
	(void)(ninjaE4);
	(void)(ninjaP0);
	(void)(ninjaP1);
	(void)(ninjaP2);
	(void)(ninja_mindeg);
	(void)(ninja_cut);
	(void)(ninja_cutidx);
	(void)(ninjaC);
	{
#ninja_numerator_t2expansion1
	}
	if (ninja_mindeg > 1+NINJA_RANK_MINUS_N) {
#ninja_numerator_t2expansion2
	}
  }

#if NINJA_NUM_NAMESPACE
} // namespace NINJA_NUM_NAMESPACE
#endif
