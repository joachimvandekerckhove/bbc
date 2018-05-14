#include <config.h>
#include "DCT.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>
#include <JRmath.h>
#include <stdio.h>

using std::vector;
using std::log;
using std::min;
using std::max;
using std::string;

static inline double MU    (vector<double const*> const &par) { return *par[0]; }
static inline double LAMBDA(vector<double const*> const &par) { return *par[1]; }
static inline double DF    (vector<double const*> const &par) { return *par[2]; }
static inline double BIAS  (vector<double const*> const &par) { return *par[3]; }
static inline double SIGMA (vector<double const*> const &par) { return    1.0 ; }
static inline double TAU   (vector<double const*> const &par) { return    1.0 ; }
static inline double DELTA (vector<double const*> const &par) { return MU(par)/SIGMA(par); }

#define ALPHA 0.05

namespace jags {
	namespace ct {

		DCT::DCT() 
			: RScalarDist("dct", 4, DIST_UNBOUNDED)
		{}

//		DCT::DCT(string const &name, unsigned int npar) : ScalarDist(name, npar, DIST_UNBOUNDED)
//		{}

		bool DCT::checkParameterValue (vector<double const *> const &par) const
		{
			bool wat;
//			printf("In DCT::checkParameterValue\n");
//			printf("%5s: %f\n", "MU"   , MU   (par));
//			printf("%5s: %f\n", "TAU"  , TAU  (par));
//			printf("%5s: %f\n", "DF"   , DF   (par));
//			printf("%5s: %f\n", "BIAS" , BIAS (par));
//			printf("%5s: %f\n", "SIGMA", SIGMA(par));
//			printf("%5s: %f\n", "DELTA", DELTA(par));
			wat = (DF(par)   > 0 &&
			       BIAS(par) <= 1 &&
			       BIAS(par) >= 0 &&
			       fabs(DELTA(par)) <= 37.62 &&
			       TAU(par) > 0 &&
			       LAMBDA(par) >= 0);
			return wat;
		}

		double DCT::d(double x, PDFType type,
					  vector<double const *> const &par, bool give_log) const
		{
//			printf("In DCT::d\n");
//			printf("%5s: %f\n", "x"    , x         );
//			printf("%5s: %f\n", "MU"   , MU   (par));
//			printf("%5s: %f\n", "TAU"  , TAU  (par));
//			printf("%5s: %f\n", "DF"   , DF   (par));
//			printf("%5s: %f\n", "BIAS" , BIAS (par));
//			printf("%5s: %f\n", "SIGMA", SIGMA(par));
//			printf("%5s: %f\n", "DELTA", DELTA(par));

			x /= SIGMA(par);

			double h0lower = pnt(x, DF(par), 0.0, 1, 0);

//			printf("p would be %f\n", h0lower);

			double crlower, crupper;
			crlower = qnt(ALPHA / 2.0, DF(par), 0, 1, 0);
			crupper = -crlower;

//			printf("q would be (%f, %f)\n", crlower, crupper);
			double h1lower, h1upper, h1middle, area;

			if (DELTA(par)) {
				h1lower = pnt(crlower, DF(par), DELTA(par), 1, 0);
				h1upper = pnt(crupper, DF(par), DELTA(par), 0, 0);
				h1middle = 1 - h1lower - h1upper;
				area = h1lower + h1upper + BIAS(par) * h1middle;
			} else {
				h1lower = ALPHA / 2.0;
				h1upper = ALPHA / 2.0;
				h1middle = 1 - ALPHA;
				area = h1lower + h1upper + BIAS(par) * h1middle;
			}

//			printf("lower tail would be %f\n", h1lower);
//			printf("upper tail would be %f\n", h1upper);
//			printf("center body would be %f\n", h1middle);
//			printf("total area would be %f\n", area);

			bool in_tails = (h0lower <= ALPHA/2.0 || h0lower >= (1 - ALPHA/2.0) );

			double dv;

			if (give_log) {
				dv = dnt(x, DF(par), DELTA(par), 1) - log(SIGMA(par));
//				printf("raw log-density is %f\n", dv);
				dv -= log(area);
				if (!in_tails) {
//					printf("... but this is in the body\n");
					dv += log(BIAS(par));
				}
//				printf("corrected log-density is %f\n", dv);
			}
			else {
				dv = dnt(x, DF(par), DELTA(par), 0) / SIGMA(par);
//				printf("raw density is %f\n", dv);
				dv /= area;
				if (!in_tails) {
//					printf("... but this is in the body\n");
					dv *= BIAS(par);
				}
//				printf("corrected density is %f\n", dv);
			}
			return dv;
		}

		double DCT::p(double x, vector<double const *> const &par, bool lower,
					  bool use_log) const
		{
			printf("In DCT::p... NOT YET IMPLEMENTED\n");
			return pnt(x/SIGMA(par), DF(par), DELTA(par), lower, use_log);
		}

		double DCT::q(double p, vector<double const *> const &par, bool lower,
					  bool log_p) const
		{
			printf("In DCT::q... NOT YET IMPLEMENTED\n");
			return qnt(p, DF(par), DELTA(par), lower, log_p) * SIGMA(par);
		}

		double DCT::r(vector<double const *> const &par, RNG *rng) const
		{
			printf("In DCT::r... NOT YET IMPLEMENTED\n");
			double k = DF(par);
			return rnorm(MU(par), SIGMA(par), rng)/sqrt(rchisq(k, rng)/k);
		}

/*		double DCT::logDensity(double x, PDFType type,
			vector<double const *> const &parameters,
			double const *lower, double const *upper) const
		{
		return d(x, type, parameters, true);
		}

		double DCT::randomSample(vector<double const *> const &parameters,
				double const *lower, double const *upper, RNG *rng) const
		{
			return r(parameters, rng);
		}

		double DCT::typicalValue(vector<double const *> const &parameters,
				double const *lower, double const *upper) const
		{
			return DELTA(parameters);
		}
*/
	} //namespace ct
} //namespace jags
