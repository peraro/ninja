// ttbarh.cc

// diagram 266, helicity4

#include <iostream>
#include <ninja/ninja.hh>
#include <ninja/rambo.hh>
#include "ttbarh_num.hh"
using namespace ninja;

int main()
{
  // External legs of the loop
  const int N_LEGS = 5;

  // Centre of mass energy
  const Real ENERGY_CM = 500;

  // Invariant s
  const Real S = ENERGY_CM*ENERGY_CM;

  // Rank of the numerator
  const int RANK = 4;

  // Top and Higgs Mass
  const Real mT = TTbarHDiagram::mT;
  const Real mH = TTbarHDiagram::mH;

  // Declare an instance of the numerator
  TTbarHDiagram num;

  // Define external momenta
  RealMomentum k[N_LEGS];

  // Create an istance of Rambo phase-space generator
  Real external_masses[] = {0,0,mH,mT,mT};
  Rambo phase_space(S,N_LEGS, external_masses);

  // get a random phase-space point
  phase_space.getMomenta(k);

  // define the internal momenta of the loop
  num.init(k);
  //num.printAbbr();

  // create an amplitude object
  ninja::Amplitude<RealMasses> amp(N_LEGS,RANK,
                                   num.getInternalMomenta(),
                                   num.getInternalMasses());
  amp.setRenormalizationScale(mT*mT);
  setTest(Test::GLOBAL);
  setVerbosity(Verbose::ALL);

  // evaluate the integral
  amp.evaluate(num);

  std::cout.precision(15);
  for (int i=0; i<5; ++i) {
    std::cout << "vecs(" << i+1 << ",:) = (/";
    for (int j=0; j<3; ++j)
      std::cout << k[i][j] << "_ki, ";
    std::cout << k[i][3] << "_ki/)\n";
  }

  // print the result
  std::cout << "Finite:      " << amp[0] // or amp.eps0(),  finite part
            << std::endl
            << "Single pole: " << amp[1] // or amp.epsm1(), single-pole
            << std::endl
            << "Double pole: " << amp[2] // or amp.epsm2(), double-pole
            << std::endl;

  return 0;
}
