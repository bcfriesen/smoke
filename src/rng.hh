#ifndef RNG_HH
#define RNG_HH

double get_rng(const bool cached);
void gen_rng_cache(double* rng_cache);

extern long int rng_count;
extern const long int rng_cache_sz;
extern double* rng_cache;
extern int num_times_regen; // how many times did we regenerate the RNG cache?
extern const bool use_cached_rng;
extern svrng_engine_t        engine;
extern svrng_distribution_t  distr1;

#endif
