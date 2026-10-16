/**
 * @file Aeff.cxx
 * @brief Implementation for post-handoff review effective area class.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include "astro/SkyDir.h"

#include "tip/TipException.h"

#include "irfUtil/IrfHdus.h"

#include "latResponse/FitsTable.h"
#include "latResponse/ParTables.h"

#include "latResponse/Aeff.h"

namespace latResponse {

Aeff::Aeff(const irfUtil::IrfHdus & irf_hdus, size_t iepoch, size_t nrow) 
   : m_aeffTable(irf_hdus("EFF_AREA").at(iepoch).first,
                 irf_hdus("EFF_AREA").at(iepoch).second,
                 "EFFAREA", nrow),
     m_phiDepPars(0) {
   m_usePhiDependence = false;
   try {
      m_phiDepPars = new ParTables(irf_hdus("PHI_DEP").at(iepoch).first,
                                   irf_hdus("PHI_DEP").at(iepoch).second,
                                   nrow);
      m_usePhiDependence = true;
   } catch (tip::TipException & eobj) {
   }
}

Aeff::Aeff(const std::string & fitsfile, const std::string & extname,
           size_t nrow)
   : m_aeffTable(fitsfile, extname, "EFFAREA", nrow),
     m_phiDepPars(0) {
   try {
      m_phiDepPars = new ParTables(fitsfile, "PHI_DEPENDENCE", nrow);
   } catch (tip::TipException & eObj) {
      m_phiDepPars = 0;
   }
   if (m_phiDepPars) {
      m_usePhiDependence = true;
   } else {
      m_usePhiDependence = false;
   }
}

Aeff::Aeff(const Aeff & other) 
   : irfInterface::IAeff(other),
     m_aeffTable(other.m_aeffTable),
     m_phiDepPars(0) {
   if (other.m_phiDepPars != 0) {
      m_phiDepPars = new ParTables(*other.m_phiDepPars);
   }
}

auto Aeff::operator=(const Aeff & rhs) -> Aeff & {
   if (this != &rhs) {
      irfInterface::IAeff::operator=(rhs);
      m_aeffTable = rhs.m_aeffTable;
      delete m_phiDepPars;
      m_phiDepPars = 0;
      if (rhs.m_phiDepPars != 0) {
         m_phiDepPars = new ParTables(*rhs.m_phiDepPars);
      }
   }
   return *this;
}

Aeff::~Aeff() {
   delete m_phiDepPars;
}

auto Aeff::value(double energy, 
                   const astro::SkyDir & srcDir, 
                   const astro::SkyDir & scZAxis,
                   const astro::SkyDir & scXAxis,
                   double time) const -> double {
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   CLHEP::Hep3Vector yhat(scZAxis().cross(scXAxis()));
   double phi = 180./M_PI*std::atan2(yhat.dot(srcDir()), 
                                     scXAxis().dot(srcDir()));
   return value(energy, theta, phi, time);
}

auto Aeff::value(double energy, double theta, double phi,
                   double time) const -> double {
   (void)(time);
   double costheta(std::cos(theta*M_PI/180.));
   if (costheta < m_aeffTable.minCosTheta()) {
      return 0;
   }
   if (costheta > 0.99999) {
      costheta = 0.99999;  // blech.
   }
   bool interpolate;
   double logE(std::log10(energy));
   double aeff_value(m_aeffTable.value(logE, costheta, interpolate=true)*1e4);
   double phi_mod(phi_modulation(logE, costheta, phi, interpolate=false));
   return aeff_value*phi_mod;
}

auto Aeff::phi_modulation(double par0, double par1, double phi) const -> double {
   if (phi < 0) {
      phi += 360.;
   }
   double norm(1./(1. + par0/(1. + par1)));
   double phi_pv(std::fmod(phi*M_PI/180., M_PI) - M_PI/2.);
   double xx(2.*std::fabs(2./M_PI*std::fabs(phi_pv) - 0.5));
   return norm*(1. + par0*std::pow(xx, par1));
}

auto Aeff::phi_modulation(double logE, double costheta, 
                            double phi, bool interpolate) const -> double {
   double par[2];
   if (!m_phiDepPars || !m_usePhiDependence) {
      return 1.;
   }
   m_phiDepPars->getPars(logE, costheta, par, interpolate);
   return phi_modulation(par[0], par[1], phi);
}

auto Aeff::max_phi_modulation() const -> double {
   if (!m_phiDepPars || !m_usePhiDependence) {
      return 1.;
   }
   const std::vector<std::string> & parNames(m_phiDepPars->parNames());
   std::vector<double> par0;
   m_phiDepPars->getParVector(parNames.at(0), par0);
   std::vector<double> par1;
   m_phiDepPars->getParVector(parNames.at(1), par1);

   double max_mod(0);
   for (size_t i(0); i < par0.size(); i++) {
      double mod0(phi_modulation(par0.at(i), par1.at(i), M_PI/4.));
      if (mod0 > max_mod) {
         max_mod = mod0;
      }
      double mod1(phi_modulation(par0.at(i), par1.at(i), M_PI/2.));
      if (mod1 > max_mod) {
         max_mod = mod1;
      }
   }
   return max_mod;
}

auto 
Aeff::pars(double logE, double costh, bool interpolate) const -> std::pair<double, double> {
   if (!m_phiDepPars || !m_usePhiDependence) {
      return std::make_pair(1., 0);
   }
   double par[2];
   m_phiDepPars->getPars(logE, costh, par, interpolate);
   return std::make_pair(par[0], par[1]);
}

auto Aeff::upperLimit() const -> double {
   return m_aeffTable.maximum()*1e4*max_phi_modulation();
}

} // namespace latResponse
