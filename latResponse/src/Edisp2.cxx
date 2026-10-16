/** 
* @file Edisp2.cxx
* @brief Implementation for post-handoff review energy dispersion class.
* @author J. Chiang
*
* $Header$
*/

#include <cmath>
#include <cstdlib>

#include <algorithm>

#include "astro/SkyDir.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/GaussianQuadrature.h"

#include "latResponse/FitsTable.h"
#include "latResponse/IrfLoader.h"

#include "Edisp2.h"

namespace {
   auto binIndex(double x, const std::vector<double> & xx) -> size_t {
      std::vector<double>::const_iterator ix = 
         std::upper_bound(xx.begin(), xx.end(), x);
      return ix - xx.begin();
   }
}

namespace latResponse {

Edisp2::Edisp2(const std::string & fitsfile, 
               const std::string & extname, size_t nrow) 
   : m_parTables(fitsfile, extname, nrow), m_loge_last(0), m_costh_last(0),
     m_renormalized(false), m_fitsfile(fitsfile), m_extname(extname),
     m_nrow(nrow), m_interpolator(0) {
   readScaling(fitsfile);
}

Edisp2::~Edisp2() {
   delete m_interpolator;
}

void Edisp2::renormalize(double logE, double costh, double * params) const {
   double energy(std::pow(10., logE));
   double scale_factor(scaleFactor(logE, costh));
   EdispIntegrand foo(params, energy, scale_factor, *this);
   double err(1e-7);
   int ierr;
   double norm(st_facilities::GaussianQuadrature::dgaus8(foo, energy/10.,
                                                         energy*10, err, ierr));
   params[0] /= norm;
}

void Edisp2::renormalize() {
   const std::vector<double> & logEnergies(m_parTables.logEnergies());
   const std::vector<double> & costhetas(m_parTables.costhetas());
   for (size_t k(0); k < logEnergies.size(); k++) {
      double loge(logEnergies[k]);
      double energy(std::pow(10., loge));
      for (size_t j(0); j < costhetas.size(); j++) {
         double costh(costhetas[j]);
         std::vector<double> my_pars;
         m_parTables.getPars(k, j, my_pars);
         renormalize(loge, costh, &my_pars[0]);
         m_parTables.setPars(k, j, my_pars);
      }
   }
   m_renormalized = true;
}

auto Edisp2::value(double appEnergy,
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

auto Edisp2::evaluate(double emeas, double energy,
                        double theta, double phi, double time,
                        double * pars) const -> double {
   (void)(phi);
   (void)(time);
   double xx((emeas - energy)/energy);
   double costh(std::cos(theta*M_PI/180.));
   costh = std::min(costh, m_parTables.costhetas().back());
   double scale_factor(scaleFactor(std::log10(energy), costh));
   xx /= scale_factor;
   return old_function(xx, pars)/energy/scale_factor;
}

auto Edisp2::value(double appEnergy, double energy,
                     double theta, double phi, double time) const -> double {
   if (::getenv("DISABLE_EDISP_INTERP")) {
      double costh(std::cos(theta*M_PI/180.));
      costh = std::min(costh, m_parTables.costhetas().back());
      double * my_pars(pars(energy, costh));
      return evaluate(appEnergy, energy, theta, phi, time, my_pars);
   }
   if (m_interpolator == 0) {
      m_interpolator = new EdispInterpolator(m_fitsfile, m_extname, m_nrow);
   }
   return m_interpolator->evaluate(*this, appEnergy, energy,
                                   theta, phi, time);
}

auto Edisp2::old_function(double xx, double * pars) const -> double {
// See ::edisp_func in handoff_response/src/gen/Dispersion.cxx
   double tt(std::fabs(xx - pars[3]));
   double s1(pars[1]);
   double s2(pars[4]);
   if (xx > pars[3]) {
      s1 = pars[2];
      s2 = pars[5];
   }
   double g1(std::exp(-0.5*std::pow(tt/s1, m_p1)));
   double g2(std::exp(-0.5*std::pow(tt/s2, m_p2)));
   double nscale(std::exp(0.5*(std::pow(m_t0/s2, m_p2) - 
                               std::pow(m_t0/s1, m_p1))));

   if (tt > m_t0) {
      return pars[0]*nscale*g2;
   }
   return pars[0]*g1;
}

auto Edisp2::scaleFactor(double logE, double costh) const -> double {
   if (!IrfLoader::interpolate_edisp() && m_interpolator==0) {
      // Use midpoint of logE, costh bins in the FITS tabulations
      // instead of the passed values.  This ensures correct
      // normalization via the renormalize() member function. Note
      // that the par values are not interpolated anyways in the
      // pars(...) member function.
      const std::vector<double> & ebounds(m_parTables.ebounds());
      size_t k = binIndex(logE, ebounds) - 1;
      if (k == ebounds.size()-1) {
         k = ebounds.size() - 2;
      }
      const std::vector<double> & tbounds(m_parTables.tbounds());
      size_t j = binIndex(costh, tbounds) - 1;
      if (j == -1) {
         j = 0;
      }

      logE = m_parTables.logEnergies()[k];
      costh = m_parTables.costhetas()[j];
   }

// See handoff_response::Dispersion::scaleFactor
   costh = std::fabs(costh);
   double my_value(m_scalePars.at(0)*logE*logE +
                   m_scalePars.at(1)*costh*costh + 
                   m_scalePars.at(2)*logE + 
                   m_scalePars.at(3)*costh +
                   m_scalePars.at(4)*logE*costh +
                   m_scalePars.at(5));
   return my_value;
}

auto Edisp2::pars(double energy, double costh) const -> double * {
   if (!m_renormalized) {
      const_cast<Edisp2 *>(this)->renormalize();
   }
   double loge(std::log10(energy));
   if (!IrfLoader::interpolate_edisp()) {
      // Ensure use of highest cos(theta) bin.
      costh = std::min(m_parTables.costhetas().back(), costh);
   }

   if (loge == m_loge_last && costh == m_costh_last) {
      return m_pars;
   }
   
   m_loge_last = loge;
   m_costh_last = costh;
   
   bool interpolate;
   m_parTables.getPars(loge, costh, m_pars, interpolate=false);

   if (IrfLoader::interpolate_edisp()) {
      // Ensure proper normalization
      EdispIntegrand foo(m_pars, energy, scaleFactor(loge, costh), *this);
      double err(1e-5);
      int ierr;
      double norm = 
         st_facilities::GaussianQuadrature::dgaus8(foo, energy/10.,
                                                   energy*10., err, ierr);
      m_pars[0] /= norm;
   }

   return m_pars;
}

void Edisp2::readScaling(const std::string & fitsfile, 
                         const std::string & extname) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));

   std::vector<double> values;

   FitsTable::getVectorData(table, "EDISPSCALE", values);

   size_t npars(values.size() - 3);
   m_scalePars.resize(npars);
   std::copy(values.begin(), values.begin() + npars, m_scalePars.begin());

   m_p1 = values.at(npars);
   m_p2 = values.at(npars + 1);
   m_t0 = values.at(npars + 2);

   delete table;
}

} // namespace latResponse
