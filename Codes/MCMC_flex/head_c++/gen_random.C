#include <boost/random.hpp>
#include <ctime>
#include <random>

using namespace boost;
 
double gen_random_normal(double mean, double sigma)
{
    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
    static std::mt19937 rng(static_cast<unsigned> (std::time(0)));
 
    // select Gaussian probability distribution
    boost::normal_distribution<double> norm_dist(mean, sigma);
 
    // bind random number generator to distribution, forming a function
    variate_generator<std::mt19937&, boost::normal_distribution<double> >  normal_sampler(rng, norm_dist);
 
    // sample from the distribution
    return normal_sampler();
}

double gen_random_unif01(void)
{
  static std::mt19937 rng(static_cast<unsigned> (std::time(0)));
  static uniform_01<std::mt19937> unif_sampler(rng);
  return unif_sampler();
}

