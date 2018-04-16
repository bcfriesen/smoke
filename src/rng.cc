#include <svrng.h>

// Cached RNG. Consume the next RNG from a pre-populated array.
double get_rng(const bool cached, long int cache_index) {
  extern double* rng_cache;
  return (rng_cache[cache_index++]);
}

// On-the-fly (traditional) RNG from Intel MKL.
double get_rng(svrng_engine_t& engine, svrng_distribution_t& distr1) {
  return svrng_generate_double( engine, distr1 );
}
