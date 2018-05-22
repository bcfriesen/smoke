#include <svrng.h>
#include "rng.hh"

// Generate cache of random numbers.
void gen_rng_cache(double* rng_cache,
                   svrng_engine_t& engine,
                   svrng_distribution_t& distr1) {
  size_t i;
  for (i = 0; i < rng_cache_sz; i++) {
    rng_cache[i] = svrng_generate_double( engine, distr1 );
  }
}

// Cached RNG. Consume the next RNG from a pre-populated array.
double get_rng(const bool cached,
               long int& cache_index,
               svrng_engine_t& engine,
               svrng_distribution_t& distr1) {
  if (cache_index == rng_cache_sz-1) {
    gen_rng_cache(rng_cache, engine, distr1);
    cache_index = 0;
  }
  return (rng_cache[cache_index++]);
}

// On-the-fly (traditional) RNG from Intel MKL.
double get_rng(svrng_engine_t& engine, svrng_distribution_t& distr1) {
  return svrng_generate_double( engine, distr1 );
}
