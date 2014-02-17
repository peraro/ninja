// Definition of the constructors of the Spinor class.


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ninja/spinors.hh>
using namespace ninja;

namespace ninja {
  
  Spinor::Spinor(const ComplexMomentum & p): ap(0.), ame(0.), bp(0.), bme(0.)
  {
	if (taxicab_norm(p[0] + p[3]) < SPINOR_TOLLERANCE) {
	  if (taxicab_norm(p[1]+I*p[2]) < SPINOR_TOLLERANCE) {
		ame = TWO*p[0];
		bme = ONE;
		ap = TWO*p[1];
	  } else if (taxicab_norm(p[0]) < SPINOR_TOLLERANCE) {
		bp = ONE;
		ame = TWO*p[1];
	  } else {
		ame = TWO*p[0];
		bme = ONE;
		bp = p[1]/p[0];
	  }
	} else {
	  ap = sqrt(p[0]+p[3]);  // sqrt(p^+)
	  ame = (p[1]+I*p[2])/ap;
	  bp = ap;  // sqrt(p^+)
	  bme= (p[1]-I*p[2])/ap;  // sqrt(p^-)
	}
  }


  Spinor::Spinor(const RealMomentum & p): ap(0.), ame(0.), bp(0.), bme(0.)
  {
	if (taxicab_norm(p[0] - p[3]) < SPINOR_TOLLERANCE) {
	  ap = bp = (p[0] < 0. ? I*std::sqrt(-TWO*p[0]) : std::sqrt(TWO*p[0]));
	} else if (taxicab_norm( p[0] + p[3]) < SPINOR_TOLLERANCE) {
	  ame = bme = (p[0] < 0. ? I*std::sqrt(-TWO*p[0]) : std::sqrt(TWO*p[0]));
	} else {
	  // sqrt(p^+)
	  ap = (p[0] < 0. ? I*std::sqrt(-p[0]-p[3])	: std::sqrt(p[0]+p[3]));
	  ame = (p[1]+I*p[2])/ap;
	  // sqrt(p^+)
	  bp = ap;
	  // sqrt(p^-)
	  bme = (p[1]-I*p[2])/ap;
	}
  }

} // namespace ninja
