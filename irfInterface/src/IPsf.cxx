/**
 * @file IPsf.cxx
 * @brief Provide default implementations of factorable member functions.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <algorithm>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "astro/SkyDir.h"

#include "st_facilities/GaussianQuadrature.h"
#include "st_facilities/RootFinder.h"

#include "irfInterface/AcceptanceCone.h"
#include "irfInterface/IPsf.h"

namespace irfInterface {

double IPsf::s_energy(1e3);
double IPsf::s_theta(0);
double IPsf::s_phi(0);
double IPsf::s_time(0);
const IPsf * IPsf::s_self(0);

double IPsf::s_cp(0);
double IPsf::s_sp(0);
double IPsf::s_cr(0);

std::vector<double> IPsf::s_psi_values;

IPsf::IPsf() {
   if (s_psi_values.size() == 0) {
      fill_psi_values();
   }
}

std::vector<double> IPsf::value(const std::vector<double>& separation,
				const std::vector<double>& energy,
				const std::vector<double>& theta,
				double phi, double time) const {
  if(separation.size() != energy.size() || separation.size() != theta.size())
    throw std::runtime_error("Input arrays must have same dimension.");

  std::vector<double> vals;
  vals.reserve(energy.size());
  std::vector<double>::const_iterator itr0 = separation.begin();
  std::vector<double>::const_iterator itr1 = energy.begin();
  std::vector<double>::const_iterator itr2 = theta.begin();  
  for(; (itr0 != separation.end()) && (itr1 != energy.end()) && (itr2 != theta.end()); 
      ++itr0, ++itr1, ++itr2) {
    vals.push_back(value(*itr0,*itr1,*itr2,phi,time));
  }
  return vals;
}

astro::SkyDir IPsf::appDir(double energy,
                           const astro::SkyDir & srcDir,
                           const astro::SkyDir & scZAxis,
                           const astro::SkyDir & scXAxis,
                           double time) const {
   (void)(scXAxis);

// Compute source inclination
   double theta(srcDir.difference(scZAxis)*180./M_PI);

   static double phi(0);
   setStaticVariables(energy, theta, phi, time, this);

   // Draw offset angle.
   double psi;
   if (::getenv("USE_OLD_IPSF_SAMPLER")) {
      // Form cumlative distribution in psi (polar angle from source direction)
      std::vector<double> integrand;
      for (std::vector<double>::iterator psi_it(s_psi_values.begin());
           psi_it != s_psi_values.end(); ++psi_it) {
         integrand.push_back(coneIntegrand(&(*psi_it)));
      }

      std::vector<double> integralDist;
      integralDist.push_back(0);
      for (size_t i = 1; i < s_psi_values.size(); i++) {
         integralDist.push_back(integralDist.at(i-1) + 
                                (integrand.at(i) + integrand.at(i-1))/2.
                                *(s_psi_values.at(i) - s_psi_values.at(i-1)));
      }
   
      double xi(CLHEP::RandFlat::shoot()*integralDist.back());
      size_t indx(std::upper_bound(integralDist.begin(), integralDist.end(), xi)
                  - integralDist.begin() - 1);
      psi = ((xi - integralDist.at(indx))
             /(integralDist.at(indx+1) - integralDist.at(indx))
             *(s_psi_values.at(indx+1) - s_psi_values.at(indx))
             + s_psi_values.at(indx));
   } else {
      const std::vector<double> & xx(s_psi_values);
      std::vector<double> yy;
      std::vector<double> aa;
      std::vector<double> bb;
      yy.push_back(coneIntegrand(&s_psi_values[0]));
      std::vector<double> integralDist;
      integralDist.push_back(0);
      for (size_t i(1); i < s_psi_values.size(); i++) {
         yy.push_back(coneIntegrand(&s_psi_values[i]));
         aa.push_back((yy[i] - yy[i-1])/(xx[i] - xx[i-1]));
         bb.push_back(yy[i-1] - xx[i-1]*aa.back());
         double value(integralDist.back()
                      + aa.back()*(xx[i]*xx[i] - xx[i-1]*xx[i-1])/2. 
                      + bb.back()*(xx[i] - xx[i-1]));
         integralDist.push_back(value);
      }
      double xi(CLHEP::RandFlat::shoot()*integralDist.back());
      size_t indx(std::upper_bound(integralDist.begin(), integralDist.end(), xi)
                  - integralDist.begin() - 1);
      double AA(aa[indx]/2.);
      double BB(bb[indx]);
      double CC(integralDist[indx] - xi 
                - aa[indx]*xx[indx]*xx[indx]/2. - bb[indx]*xx[indx]);
      psi = (std::sqrt(BB*BB - 4.*AA*CC) - BB)/2./AA;
   }
      
   double azimuth(2.*M_PI*CLHEP::RandFlat::shoot());
   
   astro::SkyDir appDir(srcDir);
   astro::SkyDir zAxis(scZAxis);
   
   CLHEP::Hep3Vector srcVec(appDir());
   CLHEP::Hep3Vector arbitraryVec(srcVec.x() + 1., srcVec.y() + 1.,
                                  srcVec.z() + 1.);
   CLHEP::Hep3Vector xVec(srcVec.cross(arbitraryVec.unit()));

   appDir().rotate(psi*M_PI/180., xVec).rotate(azimuth, srcVec);

   return appDir;
}

double IPsf::angularIntegral(double energy, double theta, 
                             double phi, double radius, double time) const {
   setStaticVariables(energy, theta, phi, time, this);
   double integral;
   double err(1e-5);
   long ierr(0);
   double zero(0);
   integral = st_facilities::GaussianQuadrature::integrate(&coneIntegrand,
                                                           zero, radius,
                                                           err, ierr);
   return integral;
}

std::vector<double> IPsf::angularIntegral(const std::vector<double>& energy, 
					  const std::vector<double>& theta, 
					  double phi, double radius, double time) const {
  if(energy.size() != theta.size())
    throw std::runtime_error("Input arrays must have same dimension.");

  std::vector<double> vals;
  vals.reserve(energy.size());
  std::vector<double>::const_iterator itr0 = energy.begin();
  std::vector<double>::const_iterator itr1 = theta.begin();  
  for(; (itr0 != energy.end()) && (itr1 != theta.end()); ++itr0, ++itr1) {
    vals.push_back(angularIntegral(*itr0,*itr1,phi,radius,time));
  }
  return vals;
}

void IPsf::setStaticVariables(double energy, double theta, double phi,
                              double time, const IPsf * self) {
   s_energy = energy;
   s_theta = theta;
   s_phi = phi;
   s_time = time;
   s_self = self;
}

double IPsf::coneIntegrand(double * offset) {
   return s_self->value(*offset, s_energy, s_theta, s_phi, s_time)
      *std::sin(*offset*M_PI/180.)*2.*M_PI*M_PI/180.;
}

void IPsf::fill_psi_values() {
   double psi_min(1e-4);
   double psi_max(90.);
   size_t npsi(100);
   double dpsi(std::log(psi_max/psi_min)/(npsi-1));
   s_psi_values.push_back(0);
   for (size_t i = 0; i < npsi; i++) {
      s_psi_values.push_back(psi_min*std::exp(dpsi*i));
   }
}

double IPsf::angularIntegral(double energy,
                             const astro::SkyDir & srcDir,
                             const astro::SkyDir & scZAxis,
                             const astro::SkyDir & scXAxis,
                             const std::vector<AcceptanceCone *> 
                             & acceptanceCones,
                             double time) {
   return psfIntegral(this, energy, srcDir, scZAxis, scXAxis, acceptanceCones,
                      time);
}

double IPsf::angularIntegral(double energy, 
                             const astro::SkyDir & srcDir,
                             double theta, double phi,
                             const std::vector<AcceptanceCone *> 
                             & acceptanceCones,
                             double time) {
   return psfIntegral(this, energy, srcDir, theta, phi, acceptanceCones, time);
}

double IPsf::angularContainment(double energy, double theta, double phi, 
				double frac, double time, double rtol) const {
  double fmax = angularIntegral(energy,theta,phi,180.);
  IntegralFunctor fn = IntegralFunctor(*this, energy, theta, phi, time);
  return st_facilities::RootFinder::find_root(fn, 0.0, 180.0, frac*fmax, rtol);
}

std::vector<double> IPsf::angularContainment(const std::vector<double>& energy, 
					     const std::vector<double>& theta, 
					     double phi, 
					     double frac, double time, double rtol) const {
  if(energy.size() != theta.size())
    throw std::runtime_error("Input arrays must have same dimension.");

  std::vector<double> vals;
  vals.reserve(energy.size());
  std::vector<double>::const_iterator itr0 = energy.begin();
  std::vector<double>::const_iterator itr1 = theta.begin();  
  for(; (itr0 != energy.end()) && (itr1 != theta.end()); ++itr0, ++itr1) {
    vals.push_back(angularContainment(*itr0,*itr1,phi,frac,time,rtol));
  }
  return vals;
}

double IPsf::psfIntegral(IPsf * self,
                         double energy,
                         const astro::SkyDir & srcDir,
                         const astro::SkyDir & scZAxis,
                         const astro::SkyDir & scXAxis,
                         const std::vector<irfInterface::AcceptanceCone *> 
                         & acceptanceCones,
                         double time) {
   (void)(scXAxis);
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   double phi(0);
   return psfIntegral(self, energy, srcDir, theta, phi, acceptanceCones, time);
}

double IPsf::psfIntegral(IPsf * self,
                         double energy,
                         const astro::SkyDir & srcDir,
                         double theta, 
                         double phi, 
                         const std::vector<irfInterface::AcceptanceCone *> 
                         & acceptanceCones,
                         double time) {
   setStaticVariables(energy, theta, phi, time, self);
   
   const irfInterface::AcceptanceCone & roiCone(*acceptanceCones.front());
   double roi_radius(roiCone.radius()*M_PI/180.);
   double psi(srcDir.difference(roiCone.center()));
   
   double one(1.);
   double mup(std::cos(roi_radius + psi));
   double mum(std::cos(roi_radius - psi));
   
   s_cp = std::cos(psi);
   s_sp = std::sin(psi);
   s_cr = std::cos(roi_radius);

   double err(1e-5);
   long ierr(0);

   double firstIntegral(0);
   if (mum < 0.99) {
      firstIntegral = 
         st_facilities::GaussianQuadrature::integrate(&psfIntegrand1, mum, 
                                                      one, err, ierr);
   }
   
   double secondIntegral(0);
   secondIntegral = 
      st_facilities::GaussianQuadrature::integrate(&psfIntegrand2, mup, 
                                                   mum, err, ierr);

   return firstIntegral + secondIntegral;
}

double IPsf::psfIntegrand1(double * mu) {
   double sep(std::acos(*mu)*180./M_PI);
   return 2.*M_PI*s_self->value(sep, s_energy, s_theta, s_phi, s_time);
}

double IPsf::psfIntegrand2(double * mu) {
   double sep(std::acos(*mu)*180./M_PI);
   double phimin(0);
   double arg((s_cr - *mu*s_cp)/std::sqrt(1. - *mu*(*mu))/s_sp);
   if (arg >= 1.) {
      phimin = 0;
   } else if (arg <= -1.) {
      phimin = M_PI;
   } else {
      phimin = std::acos(arg);
   }
   return 2.*phimin*s_self->value(sep, s_energy, s_theta, s_phi, s_time);
}

double IPsf::IntegralFunctor::operator()(double sep) const {
  return m_psf.angularIntegral(m_energy, m_theta, m_phi, sep, m_time);
}

} // namespace irfInterface
