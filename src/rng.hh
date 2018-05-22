#ifndef RNG_HH
#define RNG_HH

double get_rng(const bool cached = false, long int cache_index = -1);
double get_rng(svrng_distribution_t& distr1, svrng_engine_t& engine);
void gen_rng_cache(double* rng_cache,
                   svrng_engine_t& engine,
                   svrng_distribution_t& distr1);

extern long int rng_count;
extern const long int rng_cache_sz;
extern double* rng_cache;
extern int num_times_regen; // how many times did we regenerate the RNG cache?

#endif
