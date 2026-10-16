/** 
* @file Edisp.cxx
* @brief Implementation for post-handoff review energy dispersion class.
* @author J. Chiang
*
* $Header$
*/

#include <cmath>

#include "astro/SkyDir.h"

#include "latResponse/FitsTable.h"

#include "Edisp.h"

namespace latResponse {

Edisp::Edisp(const std::string & fitsfile, 
             const std::string & extname, size_t nrow) 
   : m_parTables(fitsfile, extname, nrow), m_loge_last(0), m_costh_last(0) {}

auto Edisp::value(double appEnergy,
                    double energy, 
                    const astro::SkyDir & srcDir,
                    const astro::SkyDir & scZAxis,
                    const astro::SkyDir & scXAxis, 
                    double time) const -> double {
   (void)(scXAxis);
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   double phi(0);
   return value(appEnergy, energy, theta, phi, time);
}

auto Edisp::value(double appEnergy, double energy,
                    double theta, double phi, double) const -> double {
   (void)(phi);
   double costh(std::cos(theta*M_PI/180.));
   double xx((appEnergy - energy)/energy);

   /// @todo Find out why this, from handoff_response::RootEval, is
   /// required
   if (xx < -0.9) {
      return 0;
   }
   double * my_pars(pars(energy, costh));
   return old_function(xx, my_pars)/energy;
}

auto Edisp::old_function(double xx, double * pars) -> double {
   double arg(xx/pars[2]);
   double xp1(1. + xx);
   double value(pars[0]*pow(xp1,pars[1])*(arg>40 ? exp(-arg) : 1/(1+exp(arg))));
   value += pars[3]*pow(xp1,pars[4])*exp(-xp1/pars[5]); // for tails
   return value;
}

auto Edisp::pars(double energy, double costh) const -> double * {
   double loge(std::log10(energy));
   if (costh == 1.0) {
      costh = 0.9999;   // restriction from handoff_response::RootEval
   }

   if (loge == m_loge_last && costh == m_costh_last) {
      return m_pars;
   }
   
   m_loge_last = loge;
   m_costh_last = costh;
   
   bool interpolate;
   m_parTables.getPars(loge, costh, m_pars, interpolate=false);
   return m_pars;
}

} // namespace latResponse
