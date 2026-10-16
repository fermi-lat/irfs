/** 
 * @file Aeff.cxx
 * @brief Implementation for the effective Area class.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "Aeff.h"

namespace {
   class TwoLinear {
   public:
      TwoLinear(double a, double b, double x0, double y0) 
         : m_a(a), m_b(b), m_x0(x0), m_y0(y0) {}
      auto operator()(double x) const -> double {
         double value;
         if (x < m_x0) {
            value = m_a*(x - m_x0) + m_y0;
         } else {
            value = m_b*(x - m_x0) + m_y0;
         }
         if (value < 0) {
            return 0;
         }
         return value;
      }
   private:
      double m_a, m_b, m_x0, m_y0;
   };
   auto profile(double mu) -> double {
      static double a(0.32), b(0.8), c(0.18);
      if (mu > b) {
         return (mu + a)/(1. + a);
      } else if (mu < c) {
         return 0;
      } else {
         return (mu - c)*(b + a)/(b - c)/(1. + a);
      }
   }
}

namespace dc1aResponse {

auto Aeff::value(double energy, 
                   const astro::SkyDir &srcDir, 
                   const astro::SkyDir &scZAxis,
                   const astro::SkyDir &scXAxis,
                   double time) const -> double {
   (void)(scXAxis);
   double theta = srcDir.difference(scZAxis)*180./M_PI;
   return value(energy, theta, 0, time);
}

auto Aeff::value(double energy, double theta, double phi,
                   double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "dc1aResponse::Aeff::value(energy, theta, phi):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double mu(cos(theta*M_PI/180.));

   if (m_params.size() == 3) {
      return ( ::profile(mu)*m_params[0]/
               (1. + exp(-m_params[2]*log(energy/m_params[1]))) );
   } else {
      const double & p2(m_params[8]);

      ::TwoLinear p0_func(m_params[0], m_params[1], m_params[2], m_params[3]);
      ::TwoLinear p1_func(m_params[4], m_params[5], m_params[6], m_params[7]);
      double mu(-cos(theta*M_PI/180.));

      double p0(p0_func(mu));
      double p1(p1_func(mu));
      if (p0 <= 0 || p1 <= 0) {
         return 0;
      }
      return p0/(1. + exp(-p2*log(energy/p1)));
   }
}

auto Aeff::upperLimit() const -> double {
   return m_params.at(0);
}

} // namespace dc1aResponse
