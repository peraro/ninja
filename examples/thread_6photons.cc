// 6photons.cc

#include <algorithm>
#include <pthread.h>
#include <cassert>

#include <ninja/ninja.hh>
#include <ninja/rambo.hh>

#include <ninja/thread_safe_integral_library.hh>

#ifdef NINJA_USE_ONELOOP
# include <ninja/avholo.hh>
#endif

#ifdef NINJA_USE_LOOPTOOLS
# include <ninja/looptools.hh>
#endif

// Choose an integral library to use as the base one
#if defined(NINJA_USE_ONELOOP)
typedef ninja::AvHOneLoop BaseIntegralLibrary;
#elif defined(NINJA_USE_LOOPTOOLS)
typedef ninja::LoopTools BaseIntegralLibrary;
#else
# error "Sorry, this example only works with OneLoop or LoopTools enabled."
#endif

#include "6photons_num.hh"
using namespace ninja;


namespace {


  struct ThreadData {
    const char * helicities;
    const RealMomentum * momenta;
    Complex result;
  };


  // this is the function passed to each thread
  void * getAmplitude(void * thread_data_in)
  {
    ThreadData * thread_data = static_cast<ThreadData *>(thread_data_in);
    const RealMomentum * k = thread_data->momenta;
    const char * helicities = thread_data->helicities;
    thread_data->result = 0;

    // here we do the same as in the non-threaded version, for one
    // phase space point
    RealMomentum k_in[6];
    char hel_in[6];
    int permutation[6] = {0,1,2,3,4,5};
    Amplitude<Massless> amp;
    SixPhotons diagram;

    amp.setN(6).setRank(6).setCutStop(3);
    amp.setSMatrix(diagram.getSMatrix());

    do {

      if (permutation[1]>permutation[5])
        continue;
    
      for (int i=0; i<6; ++i) {
        k_in[i] = k[permutation[i]];
        hel_in[i] = helicities[permutation[i]];
      }

      diagram.init(k_in,hel_in);
      amp.setKinematics(diagram.getInternalMomenta());
      amp.evaluate(diagram);
      amp.onlyCutConstructible();

    } while (std::next_permutation(permutation+1,permutation+6));

    thread_data->result = amp[0];

    return NULL;
  }

} // namespace






int main()
{
  const int EVENTS_PER_THREAD = 10;
  const int N_THREADS = 4;
  const Real CM_ENERGY = 100;

  // define a thread safe integral library
  ThreadSafeIntegralLibrary<BaseIntegralLibrary> my_integral_lib;
  setDefaultIntegralLibrary(my_integral_lib);

  // thread data: kinematics and helicities
  ThreadData thread_data[N_THREADS];
  const char helicities[] = "+--++-";
  RealMomentum k[N_THREADS][6];

  // phase space generation
  Rambo gen(CM_ENERGY*CM_ENERGY,6);
  gen.setSeed(125);

  // print the banner and avoid race conditions
  printBanner();

  for (int i=0; i<EVENTS_PER_THREAD; ++i) {

    pthread_t threads[N_THREADS];

    for (int it=0; it<N_THREADS; ++it) {

      gen.getMomenta(k[it]);
      for (int i=2; i<6; ++i)
        k[it][i] *= -1;

      thread_data[it].momenta = k[it];
      thread_data[it].helicities = helicities;

      int thread_creation = pthread_create(&threads[it], NULL,
                                         getAmplitude,
                                         static_cast<void*>(&thread_data[it]));
      assert(thread_creation == 0);

    }

    for (int it=0; it<N_THREADS; ++it) {

      int thread_status = pthread_join(threads[it], NULL);
      assert(thread_status == 0);

      // here thread_data[it].result gives the finite part of the
      // amplitude in the computed phase space point

    }

    my_integral_lib.clearIntegralCache();
  }

  // print the result of the last computed p.s.p.
  std::cout << "Finite part:  "
            << thread_data[N_THREADS-1].result
            << std::endl;

}
