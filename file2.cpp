#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

// Function that takes two double inputs (a and b) and generates a random double in the range (a,b)

double RandomReal(double a, double b){

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  //Polling computer clock and adding time in seconds and microseconds to seed RNG
  struct timeval CurrentTime;
  gettimeofday(&CurrentTime,0);
  unsigned long seed = CurrentTime.tv_sec + CurrentTime.tv_usec;

  //Creating random number generator using gsl default, seeded with computer time
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
 
  gsl_rng_set(r, seed);

  //Take double num from rng and adjust to fit range given by (a,b)
  double num = gsl_rng_uniform(r);
  double adj_num = num*(b-a)+a;
  
  gsl_rng_free(r);
 
  return adj_num;
}

