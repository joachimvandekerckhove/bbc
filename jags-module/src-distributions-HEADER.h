/*  This class contains functions for the ct (diffusion) model */
#ifndef DCT_H_
#define DCT_H_

#include <distribution/RScalarDist.h>

namespace jags {
	namespace ct {

	/**
	 * Non-central t-distribution on k degrees of freedom, with
	 * location/non-centrality parameter parameter mu, and
	 * precision parameter tau.
	 *
	 * The non-central t has no simple closed-form expression for
	 * its density but it can be defined constructively in terms
	 * of underlying normal and chi-squared distributions.
	 *
	 * <pre>
         * X ~ dnt(mu, tau, k)
	 * X = U/sqrt(V/k))
	 * U ~ N(mu, tau)
	 * V ~ Chi2(k)
	 * </pre>
	 *
	 * The non-central t-distribution is normally only defined for
	 * tau=1, in which case mu is the non-centrality
	 * parameter. The 3-parameter form here is a scaled version of
	 * a standard non-central t with non-centrality parameter
	 * delta = mu * sqrt(tau).
	 *
	 * @short non-central t distribution
	 */

		class DCT : public RScalarDist
		{
		public:
			DCT();

    /*
     * logDensity, randomSample and typicalValue use the below defined
     * d,p,q,r functions
     *
    double logDensity(double x, PDFType type,
          std::vector<double const *> const &parameters,
          double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
       double const *lower, double const *upper,
       RNG *rng) const;
    double typicalValue(std::vector<double const *> const &parameters,
       double const *lower, double const *upper) const;
    /*
     * Checks that:
     * a > 0
     * w in intervall [0,1]
     * terr > 0
     */
			double d(double x, PDFType type,
				std::vector<double const *> const &parameters,
				bool give_log) const;
			double p(double q, std::vector<double const *> const &parameters, bool lower,
				bool give_log) const;
			double q(double p, std::vector<double const *> const &parameters, bool lower,
				bool log_p) const;
			double r(std::vector<double const *> const &parameters, RNG *rng) const;

			bool checkParameterValue(std::vector<double const *> const &parameters) const;

		};
	} //namespace ct
} //namespace jags

#endif /* DCT_H_ */
