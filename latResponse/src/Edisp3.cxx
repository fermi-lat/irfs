/** 
* @file Edisp3.cxx
* @brief Implementation for third generation energy dispersion class.
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

#include "irfUtil/IrfHdus.h"

#include "latResponse/FitsTable.h"
#include "latResponse/IrfLoader.h"

#include "latResponse/Edisp3.h"

namespace {
   auto gammln(double x) -> double{
      double tmp, sum;
      static double cof[6] =
         {76.18009173, -86.50532033, 24.01409822,
          -1.231739516, 0.120858003e-2, -0.536382e-5};
      int j;
 
      x -= 1.0;
      tmp = x + 5.5;
      tmp -= (x + 0.5) * log(tmp);
      sum = 1.0;
      for (j = 0; j <= 5; j++) {
         x += 1.0;
         sum += cof[j] / x;
      }
      return (-tmp + log(2.50662827465 * sum));
   }
} // anonymous namespace

namespace latResponse {

Edisp3::Edisp3(const irfUtil::IrfHdus & edisp_hdus, size_t iepoch,
               size_t nrow) 
   : m_fitsfile(edisp_hdus("EDISP").at(iepoch).first), 
     m_extname(edisp_hdus("EDISP").at(iepoch).second), m_nrow(nrow),
     m_parTables(m_fitsfile, m_extname, m_nrow),
     m_loge_last(0), m_costh_last(0), m_interpolator(0) {
   readScaling(edisp_hdus("EDISP_SCALING").at(iepoch).first,
               edisp_hdus("EDISP_SCALING").at(iepoch).second);
}

Edisp3::Edisp3(const std::string & fitsfile, 
               const std::string & extname, 
               const std::string & scaling_extname,
               size_t nrow) 
   : m_parTables(fitsfile, extname, nrow), m_loge_last(0), m_costh_last(0),
     m_fitsfile(fitsfile), m_extname(extname),
     m_nrow(nrow), m_interpolator(0) {
   readScaling(fitsfile, scaling_extname);
}

Edisp3::~Edisp3() {
   delete m_interpolator;
}

auto Edisp3::value(double appEnergy,
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

auto Edisp3::evaluate(double emeas, double energy,
                        double theta, double phi, double time,
                        double * pars) const -> double {
   (void)(phi);
   (void)(time);
   double xx((emeas - energy)/energy);
   double costh(std::cos(theta*M_PI/180.));
   costh = std::min(costh, m_parTables.costhetas().back());
   double scale_factor(scaleFactor(std::log10(energy), costh));
   xx /= scale_factor;
   return thibaut_function(xx, pars)/energy/scale_factor;
}

auto Edisp3::value(double appEnergy, double energy,
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

auto Edisp3::thibaut_function(double xx, double * pars) const -> double {
// See https://confluence.slac.stanford.edu/x/URDlCQ
// Parameter ordering in FITS file: F, S1, K1, BIAS, BIAS2, 
// S2, K2, PINDEX1, PINDEX2
   double F(pars[0]);
   double S1(pars[1]);
   double K1(pars[2]);
   double BIAS(pars[3]);
   double BIAS2(pars[4]);
   double S2(pars[5]);
   double K2(pars[6]);
   double PINDEX1(pars[7]);
   double PINDEX2(pars[8]);
   double value(F*thibaut_base_function(xx, S1, K1, BIAS, PINDEX1)
                + (1 - F)*thibaut_base_function(xx, S2, K2, BIAS2, PINDEX2));
   return value;
}

auto Edisp3::thibaut_base_function(double xx, double sigma, double kk, 
                                     double bb, double pp) const -> double {
   double uu(xx - bb);
   double uup(0);
   double uum(0);
   if (uu >= 0) {
      uup = uu;
   }
   if (uu <= 0) {
      uum = -uu;
   }
   double gamma_one_over_p(std::exp(gammln(1./pp)));
   double exp_pp(std::exp(-std::pow(kk*uup/sigma, pp)));
   double exp_mp(std::exp(-std::pow(uum/kk/sigma, pp)));
   double value(pp/sigma/gamma_one_over_p*kk/(1 + kk*kk)*exp_pp*exp_mp);
   if (value < 0) {
      throw std::runtime_error("negative edisp value");
   }
   return value;
}

auto Edisp3::scaleFactor(double logE, double costh) const -> double {
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

auto Edisp3::pars(double energy, double costh) const -> double * {
   double loge(std::log10(energy));
   // if (!IrfLoader::interpolate_edisp()) {
   //    // Ensure use of highest cos(theta) bin.
   //    costh = std::min(m_parTables.costhetas().back(), costh);
   // }

   if (loge == m_loge_last && costh == m_costh_last) {
      return m_pars;
   }
   
   m_loge_last = loge;
   m_costh_last = costh;
   
   // Do not interpolate on the parameter values!
   bool interpolate;
   m_parTables.getPars(loge, costh, m_pars, interpolate=false);

   // if (IrfLoader::interpolate_edisp()) {
   //    // Ensure proper normalization
   //    EdispIntegrand foo(m_pars, energy, scaleFactor(loge, costh), *this);
   //    double err(1e-5);
   //    int ierr;
   //    double norm = 
   //       st_facilities::GaussianQuadrature::dgaus8(foo, energy/10.,
   //                                                 energy*10., err, ierr);
   //    m_pars[0] /= norm;
   // }

   return m_pars;
}

void Edisp3::readScaling(const std::string & fitsfile, 
                         const std::string & extname) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));
   FitsTable::getVectorData(table, "EDISPSCALE", m_scalePars);
   delete table;
}

void Edisp3::setParams(size_t indx, const std::vector<double>& params) {  
  interpolator().setParams(indx,params);
  std::vector<double> p(m_parTables.logEnergies().size()*
			m_parTables.costhetas().size());
  const std::vector<double>& x = interpolator().energies();
  const std::vector<double>& y = interpolator().thetas();
  for (size_t i(1); i < x.size()-1; i++) {
    for (size_t j(1); j < y.size()-1; j++) {
      size_t idx0 = j*x.size() + i;
      size_t idx1 = (j-1)*(x.size()-2) + i-1;
      p.at(idx1) = params.at(idx0);
    }
  }
  m_parTables.setParams(indx,p);
}


} // namespace latResponse
