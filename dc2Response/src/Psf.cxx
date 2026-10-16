/** 
 * @file Psf.cxx
 * @brief Implementation for the DC2 point-spread function class.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Geometry/Vector3D.h"

using CLHEP::RandFlat;
using CLHEP::Hep3Vector;

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/dgaus8.h"
#include "st_facilities/Util.h"

#include "irfInterface/AcceptanceCone.h"

#include "Psf.h"
#include "PsfScaling.h"

namespace {
   auto psfFunc(double x, double gamma) -> double {
      return x*std::pow(1. + x*x/2./gamma, -gamma);
   }
}

namespace dc2Response {

double Psf::s_energy;
double Psf::s_theta;
double Psf::s_phi;
const Psf * Psf::s_self(0);
std::vector<double> Psf::s_gammas;
std::vector<double> Psf::s_psfNorms;
std::vector<double> Psf::s_lowerFractions;
std::vector<double> Psf::s_psi;

Psf::Gint Psf::s_gfunc;

Psf::Psf(const std::string & fitsfile, const std::string & extname)
   : DC2(fitsfile, extname), m_psfScaling(0), m_haveAngularIntegrals(false),
     m_acceptanceCone(0) {
   readData();
   if (s_gammas.empty()) {
      computePsfNorms();
      computeLowerFractions();
   }
}

Psf::~Psf() {
   delete m_psfScaling;
}

Psf::Psf(const Psf & rhs) 
   : IPsf(rhs), DC2(rhs), m_psfScaling(0), m_sigma(rhs.m_sigma),
     m_gamma(rhs.m_gamma), m_logElo(rhs.m_logElo), m_logEhi(rhs.m_logEhi),
     m_logE(rhs.m_logE), m_cosinc(rhs.m_cosinc), m_angScale(rhs.m_angScale),
     m_gamValues(rhs.m_angScale), m_angularIntegral(rhs.m_angularIntegral),
     m_needIntegral(rhs.m_needIntegral),
     m_haveAngularIntegrals(rhs.m_haveAngularIntegrals), 
     m_acceptanceCone(0) {
   if (rhs.m_psfScaling) {
      m_psfScaling = new PsfScaling(*rhs.m_psfScaling);
   }
   if (rhs.m_acceptanceCone) {
      m_acceptanceCone = rhs.m_acceptanceCone->clone();
   }
}

auto Psf::value(const astro::SkyDir & appDir, 
                  double energy, 
                  const astro::SkyDir & srcDir, 
                  const astro::SkyDir & scZAxis,
                  const astro::SkyDir &,
                  double time) const -> double {
// Angle between photon and source directions in radians.
   double separation = appDir.difference(srcDir);
   
// Inclination wrt spacecraft z-axis in radians
   double inc = srcDir.difference(scZAxis);

   return value(separation*180./M_PI, energy, inc*180./M_PI, 0., time);
}

auto Psf::value(double separation, double energy, double theta,
                  double phi, double time) const -> double {
   if (theta < 0) {
      std::ostringstream message;
      message << "dc2Response::Psf::value(...):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }
   (void)(phi);
   (void)(time);
   double logE(std::log(energy));
   double mu(std::cos(theta*M_PI/180.));
   double gam(gamma(logE, mu));
   double angScale(angularScale(energy, mu));
   return value(separation, angScale, gam);
}

auto Psf::value(double separation, double angScale, double gam) const -> double {
   double psfNorm(st_facilities::Util::interpolate(s_gammas, s_psfNorms, gam));
   double x(separation/angScale);
   if (separation == 0) {
      return std::pow(1. + x*x/2./gam, -gam)/angScale*90./M_PI/M_PI
         /(angScale*M_PI/180.)/psfNorm;
   }
   return ::psfFunc(x, gam)/2./M_PI/std::sin(separation*M_PI/180.)
      /(angScale*M_PI/180.)/psfNorm;
}

auto Psf::gamma(double logE, double mu) const -> double {
// Tabulated gammas from AllGamma fits are not reliable for
// extrapolation, so force logE and mu to lie on or within grid
// boundaries.
   mu = std::max(std::min(m_cosinc.front(), mu), m_cosinc.back());
   logE = std::max(std::min(m_logE.back(), logE), m_logE.front());
   double my_gamma = 
      st_facilities::Util::bilinear(m_cosinc, mu, m_logE, logE, m_gamma);
   if (my_gamma < s_gammas.front()) {
      return s_gammas.front();
   }
   if (my_gamma > s_gammas.back()) {
      return s_gammas.back();
   }
   return my_gamma;
}

auto Psf::sigma(double logE, double mu) const -> double {
// Tabulated sigmas from AllGamma fits are not reliable for
// extrapolation, so force logE and mu to lie on or within grid
// boundaries.
   mu = std::max(std::min(m_cosinc.front(), mu), m_cosinc.back());
   logE = std::max(std::min(m_logE.back(), logE), m_logE.front());
   double my_sigma =
      st_facilities::Util::bilinear(m_cosinc, mu, m_logE, logE, m_sigma);
   assert(my_sigma > 0);
   return my_sigma;
}

auto Psf::angularScale(double energy, double mu) const -> double {
   double logE(std::log(energy));
   double scale((*m_psfScaling)(energy, mu));
   double my_value = scale*sigma(logE, mu);
   return my_value;
}

auto Psf::appDir(double energy,
                          const astro::SkyDir & srcDir,
                          const astro::SkyDir & scZAxis,
                          const astro::SkyDir &,
                          double time) const -> astro::SkyDir {
   (void)(time);
   double mu(std::cos(srcDir.difference(scZAxis)));
   double theta(drawOffset(energy, mu)*M_PI/180.);
             
   double xi = RandFlat::shoot();
   double phi = 2.*M_PI*xi;

// These might as well have been passed by value...
   astro::SkyDir appDir = srcDir;
   astro::SkyDir zAxis = scZAxis;

// Create an arbitrary x-direction about which to perform the theta
// rotation.
   Hep3Vector srcVec(appDir());
   Hep3Vector arbitraryVec(srcVec.x() + 1., srcVec.y() + 1., srcVec.z() + 1.);
   Hep3Vector xVec = srcVec.cross(arbitraryVec.unit());

   appDir().rotate(theta, xVec).rotate(phi, srcVec);

   return appDir;
}

auto Psf::drawOffset(double energy, double mu) const -> double {
   double logE(std::log(energy));
   return angularScale(energy, mu)*drawScaledDev(gamma(logE, mu));
}

auto Psf::angularIntegral(double energy,
                            const astro::SkyDir & srcDir,
                            const astro::SkyDir & scZAxis,
                            const astro::SkyDir &,
                            const std::vector<irfInterface::AcceptanceCone *> 
                            & acceptanceCones, 
                            double time) -> double {
   double theta = srcDir.difference(scZAxis)*180./M_PI;
   static double phi(0);
   return angularIntegral(energy, srcDir, theta, phi, acceptanceCones, time);
}

auto Psf::angularIntegral(double energy, 
                            const astro::SkyDir & srcDir,
                            double theta, 
                            double phi, 
                            const std::vector<irfInterface::AcceptanceCone *> 
                            & acceptanceCones, double time) -> double {
   (void)(phi);
   (void)(time);
   if (!m_acceptanceCone || *m_acceptanceCone != *(acceptanceCones[0])) {
      computeAngularIntegrals(acceptanceCones);
      m_haveAngularIntegrals = true;
   }
   double psi = srcDir.difference(m_acceptanceCone->center());
   if (psi >= s_psi.back()) {
      return 0;
   }
   size_t ipsi = std::upper_bound(s_psi.begin(), s_psi.end(), psi) 
      - s_psi.begin();
      
   double mu = std::cos(theta*M_PI/180.);
   double angScale = angularScale(energy, mu);
   if (angScale < m_angScale.front() || angScale >= m_angScale.back()) {
//       std::cout << "angular scale outside bounds: " 
//                 << angScale << "  "
//                 << m_angScale.front() << "  "
//                 << m_angScale.back() << std::endl;
      return 0;
   }
   size_t iang = std::upper_bound(m_angScale.begin(), m_angScale.end(),
                                  angScale) - m_angScale.begin() - 1;

   double gamValue(gamma(std::log(energy), mu));
   size_t igam = std::upper_bound(m_gamValues.begin(), m_gamValues.end(),
                                  gamValue) - m_gamValues.begin() - 1;

   if (iang == m_angScale.size() - 1 || igam == m_gamValues.size() - 1) {
      return psfIntegral(psi, angScale, gamValue);
   }

   size_t is[2] = {iang, iang + 1};
   size_t ig[2] = {igam, igam + 1};

   for (size_t i(0); i < 2; i++) {
      for (size_t j(0); j < 2; j++) {
         size_t indx(is[i]*m_gamValues.size() + ig[j]);
         if (m_needIntegral.at(ipsi).at(indx)) {
            m_angularIntegral.at(ipsi).at(indx) =
               psfIntegral(psi, m_angScale.at(is[i]),
                           m_gamValues.at(ig[j]));
            m_needIntegral.at(ipsi).at(indx) = false;
         }
      }
   }

   return bilinear(angScale, gamValue, ipsi, iang, igam);
}

auto Psf::angularIntegral(double energy, double theta, 
                            double phi, double radius, double time) const -> double {
   (void)(time);
   s_energy = energy;
   s_theta = theta;
   s_phi = phi;
   s_self = this;
   double integral;
   double err(1e-5);
   long ierr(0);
   double zero(0);
   dgaus8_(&coneIntegrand, &zero, &radius, &err, &integral, &ierr);
   return integral;
}

auto Psf::coneIntegrand(double * offset) -> double {
   return s_self->value(*offset, s_energy, s_theta, s_phi)
      *std::sin(*offset*M_PI/180.)*2.*M_PI*M_PI/180.;
}

void Psf::readData() {
   m_psfScaling = new PsfScaling(m_filename);

// Read in the parameters from the FITS file.
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());

   const tip::Table * psf = fileSvc.readTable(m_filename, m_extname);

   tip::Table::ConstIterator it(psf->begin());
   tip::ConstTableRecord & row(*it);

   row["sigma"].get(m_sigma);
   row["gamma"].get(m_gamma);

   std::vector<double> elo, ehi;
   row["energ_lo"].get(elo);
   row["energ_hi"].get(ehi);

   for (size_t k = 0; k < elo.size(); k++) {
      m_logElo.push_back(std::log(elo.at(k)));
      m_logEhi.push_back(std::log(ehi.at(k)));
      m_logE.push_back((m_logElo.back() + m_logEhi.back())/2.);
   }

   std::vector<double> theta;
   row["theta"].get(theta);

   for (size_t i = 0; i < theta.size(); i++) {
      m_cosinc.push_back(std::cos(theta.at(i)));
   }
   delete psf;
}

void Psf::computePsfNorms() {
   size_t ngam(100);
   double gmin(1.001);
   double gmax(10);
   double gstep(std::log(gmax/gmin)/(ngam - 1));
   double lowerLim(0);
   double upperLim(20);
   double err(1e-5);
   long ierr(0);
   double psfNorm(0);
   for (size_t i = 0; i < ngam; i++) {
      s_gammas.push_back(gmin*std::exp(i*gstep));
      dgaus8_(&psfIntegrand, &lowerLim, &upperLim, &err, &psfNorm, &ierr);
      s_psfNorms.push_back(psfNorm);
   }
}

auto Psf::drawScaledDev(double gamma) const -> double {
   double xbreak(std::sqrt(2.*gamma));
   double alpha(std::pow(2.*gamma, gamma));
   double beta(2.*gamma - 1.);

   double lowerFrac(st_facilities::Util::interpolate(s_gammas, 
                                                     s_lowerFractions,
                                                     gamma));
   double x(0);
   while (true) {
      double xi(RandFlat::shoot());
      if (xi < lowerFrac) {
         x = xbreak*std::sqrt(RandFlat::shoot());
         if (RandFlat::shoot() < ::psfFunc(x, gamma)/x) {
            return x;
         }
      } else {
         x = xbreak*std::pow((1. - RandFlat::shoot()), 1./(1. - beta));
         if (RandFlat::shoot() < 
             ::psfFunc(x, gamma)/(alpha*std::pow(x, -beta))) {
            return x;
         }
      }
   }
   return x;
}

void Psf::computeLowerFractions() {
   std::vector<double>::const_iterator it;
   s_lowerFractions.clear();
   for (it = s_gammas.begin(); it != s_gammas.end(); ++it) {
      double gamma = *it;
      double xbreak(std::sqrt(2.*gamma));
      double alpha(std::pow(2.*gamma, gamma));
      double beta(2.*gamma - 1.);
      double lowerIntegral(gamma);
      double upperIntegral(alpha/(beta - 1.)*std::pow(xbreak, 1. - beta));
      s_lowerFractions.push_back(lowerIntegral/(lowerIntegral+upperIntegral));
   }
}

auto Psf::psfIntegrand(double * xx) -> double {
   const double & x(*xx);
   return ::psfFunc(x, s_gammas.back());
}

void Psf::computeAngularIntegrals
(const std::vector<irfInterface::AcceptanceCone *> & cones) {

// Assume the first acceptance cone is the ROI; ignore the rest.
   if (!m_acceptanceCone) {
      m_acceptanceCone = cones.at(0)->clone();
   }

   if (!m_haveAngularIntegrals) {
// Set up the array describing the separation between the center of
// the ROI and the source, s_psi
      if (s_psi.empty()) {
         size_t npsi = 500;
         double psiMin = 0.;
// Sources outside psiMax will be calculated as special cases:
         double psiMax = 60.*M_PI/180.;
         double psiStep = (psiMax - psiMin)/(npsi - 1.);
         for (size_t i = 0; i < npsi; i++) {
            s_psi.push_back(psiStep*i + psiMin);
         }
      }

// The angular scale array should be logrithmically spaced because of
// the strong energy dependence of the psf and will have a range given
// by the extremal values of angularScale(energy, mu).
      size_t nangles = 40;
      m_angScale.clear();
//       double angScaleMin = angularScale(std::exp(m_logE.back()), 
//                                         m_cosinc.front());
//       double angScaleMax = angularScale(std::exp(m_logE.front()), 
//                                         m_cosinc.back());
// There are zero values for sigma in Toby's psf parameter files
// so we cannot rely on sensible values for the angularScale at the
// expected limits of energy and inclination, so we are forced to
// use ad hoc values.
//       double angScaleMin = 0.5*(*m_psfScaling)(std::exp(m_logE.back()),
//                                                m_cosinc.back());
      double angScaleMin(1e-3);
      double angScaleMax = (*m_psfScaling)(std::exp(m_logE.front()),
                                           m_cosinc.front());
      double angScaleStep = std::log(angScaleMax/angScaleMin)/(nangles - 1.);
      for (size_t j = 0; j < nangles; j++) {
         m_angScale.push_back(angScaleMin*std::exp(angScaleStep*j));
      }

// Create an ordered vector of gamma values spanning the range in the 
// parameter vector.
      double gmin(m_gamma.front());
      double gmax(m_gamma.front());
      for (size_t i = 0; i < m_gamma.size(); i++) {
         if (m_gamma.at(i) < gmin) {
            gmin = m_gamma.at(i);
         }
         if (m_gamma.at(i) > gmax) {
            gmax = m_gamma.at(i);
         }
      }
      size_t ngam(40);
      double gstep((gmax - gmin)/(ngam - 1.));
      m_gamValues.clear();
      for (size_t i = 0; i < ngam; i++) {
         m_gamValues.push_back(gmin + i*gstep);
      }
   }

// Fill m_angularIntegrals, i.e., set needIntegral flags to true for
// lazy evaluation.
   m_angularIntegral.clear();
   m_needIntegral.clear();
   size_t npts(m_angScale.size()*m_gamValues.size());
   for (size_t ipsi = 0; ipsi < s_psi.size(); ipsi++) {
      std::vector<double> drow(npts, 0);
      std::vector<bool> brow(npts, true);
      m_angularIntegral.push_back(drow);
      m_needIntegral.push_back(brow);
   }
}

auto Psf::bilinear(double angScale, double gam, size_t ipsi,
                     size_t iang, size_t igam) const -> double {
   double tt = (gam - m_gamValues.at(igam))
      /(m_gamValues.at(igam+1) - m_gamValues.at(igam));
   double uu = (angScale - m_angScale.at(iang))
      /(m_angScale.at(iang+1) - m_angScale.at(iang));
   double y1 = m_angularIntegral.at(ipsi).at(iang*m_gamValues.size()+igam);
   double y2 = m_angularIntegral.at(ipsi).at(iang*m_gamValues.size()+igam+1);
   double y3 = m_angularIntegral.at(ipsi).at((iang+1)*m_gamValues.size()+igam);
   double y4 = 
      m_angularIntegral.at(ipsi).at((iang+1)*m_gamValues.size()+igam+1);
   double value = (1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
      + tt*uu*y3 + (1. - tt)*uu*y4;
   return value;
}

auto Psf::psfIntegral(double psi, double angScale, double gamValue,
                        double roi_radius) -> double {
   if (roi_radius == 0) {
      roi_radius = m_acceptanceCone->radius()*M_PI/180.;
   }
   double one = 1.;
   double mup = cos(roi_radius + psi);
   double mum = cos(roi_radius - psi);
   
   double cp = cos(psi);
   double sp = sin(psi);
   double cr = cos(roi_radius);

   double err = 1e-5;
   long ierr;
   double firstIntegral = 0;

// Check if point is inside or outside the cone            
   if (psi < roi_radius) {
// Integrate the first term...
      s_gfunc = Gint(this, angScale, gamValue);
      dgaus8_(&Psf::gfuncIntegrand, &mum, &one, 
              &err, &firstIntegral, &ierr);
   }

// and the second.
   s_gfunc = Gint(this, angScale, gamValue, cp, sp, cr);
   double secondIntegral;
   dgaus8_(&Psf::gfuncIntegrand, &mup, &mum, 
           &err, &secondIntegral, &ierr);
   
   double my_integral;
   if (psi < roi_radius) {
      my_integral = firstIntegral + secondIntegral;
   } else {
      my_integral = secondIntegral;
   }
   return my_integral;
}

auto Psf::Gint::value(double mu) const -> double {
   double theta = acos(mu)*180./M_PI;
   double my_value = m_psfObj->value(theta, m_angScale, m_gamma);

   if (!m_doFirstTerm) {
      double phimin;
      double arg = (m_cr - mu*m_cp)/sqrt(1. - mu*mu)/m_sp;
      if (arg >= 1.) {
         phimin = 0;
      } else if (arg <= -1.) {
         phimin = M_PI;
      } else {
         phimin = acos(arg);
      }
      return 2.*phimin*my_value;
   }
   return 2.*M_PI*my_value;
}


} // namespace dc2Response
