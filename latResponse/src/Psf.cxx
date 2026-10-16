/**
 * @file Psf.cxx
 * @brief Implementation for Psf class.
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/GaussianQuadrature.h"
#include "st_facilities/Util.h"

#include "astro/SkyDir.h"

#include "Psf.h"
#include "PsfIntegralCache.h"

namespace {
   auto sqr(double x) -> double {
      return x*x;
   }
}

namespace latResponse {

double Psf::s_ub(10.);

Psf::Psf(const std::string & fitsfile, bool isFront,
         const std::string & extname, size_t nrow)
   : PsfBase(fitsfile, isFront, extname),
     m_parTables(fitsfile, extname, nrow), m_loge_last(0), m_costh_last(0), 
     m_integralCache(0) {
}

Psf::Psf(const Psf & rhs) : PsfBase(rhs), 
                            m_parTables(rhs.m_parTables), 
                            m_par0(rhs.m_par0), m_par1(rhs.m_par1),
                            m_index(rhs.m_index), m_psf_pars(rhs.m_psf_pars),
                            m_loge_last(0), m_costh_last(0), 
                            m_integralCache(0) {}
   
Psf::~Psf() {
   delete m_integralCache;
}

auto Psf::value(const astro::SkyDir & appDir, 
                  double energy, 
                  const astro::SkyDir & srcDir, 
                  const astro::SkyDir & scZAxis,
                  const astro::SkyDir & scXAxis, 
                  double time) const -> double {
   (void)(scXAxis);
   double sep(appDir.difference(srcDir)*180./M_PI);
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   double phi(0);
   return value(sep, energy, theta, phi, time);
}

auto Psf::value(double separation, double energy, double theta,
                  double phi, double time) const -> double {
   (void)(phi);
   (void)(time);
   double * my_pars(pars(energy, std::cos(theta*M_PI/180.)));
   return old_function(separation*M_PI/180., my_pars);
}

auto Psf::angularIntegral(double energy, double theta, 
                            double phi, double radius, double time) const -> double {
   if (energy < 120.) {
      double integral = IPsf::angularIntegral(energy, theta, phi, radius, time);
      return integral;
   }
   double * my_pars(pars(energy, std::cos(theta*M_PI/180.)));
   return old_integral(radius*M_PI/180., my_pars)*(2.*M_PI*::sqr(my_pars[1]));
}

auto Psf::angularIntegral(double energy,
                            const astro::SkyDir & srcDir,
                            const astro::SkyDir & scZAxis,
                            const astro::SkyDir & scXAxis,
                            const AcceptanceConeVector_t & acceptanceCones, 
                            double time) -> double {
   (void)(scXAxis);
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   double phi(0);
   return angularIntegral(energy, srcDir, theta, phi, acceptanceCones, time);
}

auto Psf::angularIntegral(double energy, 
                            const astro::SkyDir & srcDir,
                            double theta, 
                            double phi, 
                            const AcceptanceConeVector_t & acceptanceCones, 
                            double time) -> double {
   (void)(phi);
   (void)(time);
   irfInterface::AcceptanceCone & cone(*acceptanceCones.at(0));
   if (!m_integralCache || cone != m_integralCache->acceptanceCone()) {
      delete m_integralCache;
      m_integralCache = new PsfIntegralCache(*this, cone);
   }
   double psi(srcDir.difference(cone.center()));

   const std::vector<double> & psis(m_integralCache->psis());
   if (psi > psis.back()) {
      std::ostringstream message;
      message << "latResponse::Psf::angularIntegral:\n"
              << "Error evaluating PSF integral.\n"
              << "Requested source location > " 
              << psis.back()*180/M_PI << " degrees from ROI center.";
      throw std::runtime_error(message.str());
   }

   size_t ii(std::upper_bound(psis.begin(), psis.end(), psi) 
             - psis.begin() - 1);

   double * my_pars(pars(energy, std::cos(theta*M_PI/180.)));
   double ncore(my_pars[0]);
   double sigma(my_pars[1]);
   double gcore(my_pars[2]);
   double gtail(my_pars[3]);
   double ntail = ncore*(old_base_function(s_ub, sigma, gcore)
                         /old_base_function(s_ub, sigma, gtail));

   /// Remove sigma**2 scaling imposed by pars(...).  This is put back
   /// in angularIntegral below for each grid value of sigmas.  This
   /// preserves the normalization in the bilinear interpolation by
   /// explicitly putting in the important sigma-dependence.
   ncore *= sigma*sigma;
   ntail *= sigma*sigma;

   double y1(ncore*m_integralCache->angularIntegral(sigma, gcore, ii) + 
             ntail*m_integralCache->angularIntegral(sigma, gtail, ii));
   double y2(ncore*m_integralCache->angularIntegral(sigma, gcore, ii+1) + 
             ntail*m_integralCache->angularIntegral(sigma, gtail, ii+1));

   double y = ((psi - psis.at(ii))/(psis.at(ii+1) - psis.at(ii))
               *(y2 - y1)) + y1;

   return y;
}

auto Psf::old_base_function(double u, double sigma, double gamma) -> double {
   (void)(sigma);
   // ugly kluge because of sloppy programming in handoff_response
   // when setting boundaries of fit parameters for the PSF.
   if (gamma == 1) {
      gamma = 1.001;
   }
   return (1. - 1./gamma)*std::pow(1. + u/gamma, -gamma);
}

auto Psf::old_base_integral(double u, double sigma, double gamma) -> double {
   (void)(sigma);
   return 1. - std::pow(1. + u/gamma, 1. - gamma);
}

auto Psf::old_function(double sep, double * pars) -> double {
   double ncore(pars[0]);
   double sigma(pars[1]);
   double gcore(pars[2]);
   double gtail(pars[3]);
   double ntail = ncore*(old_base_function(s_ub, sigma, gcore)
                         /old_base_function(s_ub, sigma, gtail));
   double r = sep/sigma;
   double u = r*r/2.;
   return (ncore*old_base_function(u, sigma, gcore) +
           ntail*old_base_function(u, sigma, gtail));
}

auto Psf::old_integral(double sep, double * pars) -> double {
   double ncore(pars[0]);
   double sigma(pars[1]);
   double gcore(pars[2]);
   double gtail(pars[3]);
   double ntail = ncore*(old_base_function(s_ub, sigma, gcore)
                         /old_base_function(s_ub, sigma, gtail));
   double r = sep/sigma;
   double u = r*r/2.;
   return (ncore*old_base_integral(u, sigma, gcore) + 
           ntail*old_base_integral(u, sigma, gtail));
}

auto Psf::pars(double energy, double costh) const -> double * {
   double loge(std::log10(energy));
   if (costh == 1.0) {  // Why is this necessary?
      costh = 0.9999;
   }
   
   if (loge == m_loge_last && costh == m_costh_last) {
      return m_pars;
   }
   
   m_loge_last = loge;
   m_costh_last = costh;
   
   m_parTables.getPars(loge, costh, m_pars);
   
   // Rescale the sigma value after interpolation
   m_pars[1] *= scaleFactor(energy);
   
   if (m_pars[1] == 0 || m_pars[2] == 0 || m_pars[3] == 0) {
      std::ostringstream message;
      message << "latResponse::Psf::pars: psf parameters are zero "
              << "when computing solid angle normalization:\n"
              << "\tenergy = " << energy << "\n"
              << "\tm_pars[1] = " << m_pars[1] << "\n"
              << "\tm_pars[2] = " << m_pars[2] << "\n"
              << "\tm_pars[3] = " << m_pars[3] << std::endl;
      std::cerr << message.str() << std::endl;
      throw std::runtime_error(message.str());
   }
   
// Ensure that the Psf integrates to unity.
   double norm;
   static double theta_max(M_PI/2.);
   if (energy < 120.) { // Use the *correct* integral of Psf over solid angle.
      PsfIntegrand foo(m_pars);
      double err(1e-5);
      int ierr;
      norm = st_facilities::GaussianQuadrature::dgaus8(foo, 0, theta_max,
                                                       err, ierr);
      m_pars[0] /= norm*2.*M_PI;
   } else { // Use small angle approximation.
      norm = old_integral(theta_max, m_pars);
      m_pars[0] /= norm*2.*M_PI*m_pars[1]*m_pars[1];
   }

   return m_pars;
}


} // namespace latResponse

