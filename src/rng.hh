#ifndef RNG_HH
#define RNG_HH

double get_rng(const bool cached = false, long int cache_index = -1);
double get_rng(svrng_distribution_t& distr1, svrng_engine_t& engine);

extern long int rng_count;

#endif
