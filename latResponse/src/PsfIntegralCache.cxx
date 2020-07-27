/**
 * @file PsfIntegralCache.cxx
 * @brief Implementation for a class to cache angular integrals of the
 * Psf over regions of interest as a function of source offset angle,
 * and the psf parameters gamma and sigma.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "st_facilities/GaussianQuadrature.h"

#include "st_stream/StreamFormatter.h"

#include "irfInterface/AcceptanceCone.h"

#include "latResponse/PsfBase.h"
#include "Psf.h"
#include "PsfIntegralCache.h"

namespace latResponse {

PsfIntegralCache::
PsfIntegralCache(const PsfBase & psf, irfInterface::AcceptanceCone & cone) 
   : m_psf(psf), m_acceptanceCone(cone.clone()), 
     m_calls(0), m_interpolations(0), m_cpuTotal(0),
     m_gamma_avg(0), m_sigma_avg(0), m_integralEvals(0),
     m_gamma_max(0), m_gamma_min(100), m_sigma_max(0), m_sigma_min(100) {
   fillParamArrays();
   setupAngularIntegrals();
}

PsfIntegralCache::~PsfIntegralCache() {
   st_stream::StreamFormatter formatter("latResponse", "", 3);
   formatter.info() << "PsfIntegralCache calls: " << m_calls << "\n"
                    << "PsfIntegralCache interpolations: " 
                    << m_calls - m_interpolations << std::endl
                    << "CPU time in psfIntegral (s): " << m_cpuTotal << "\n"
                    << "Avg. CPU time per integral (s): " 
                    << m_cpuTotal/m_interpolations << std::endl
                    << "Avg. gamma: " << m_gamma_avg/m_integralEvals << "\n"
                    << "Min. gamma: " << m_gamma_min << "\n"
                    << "Max. gamma: " << m_gamma_max << "\n"
                    << "Avg. sigma: " << m_sigma_avg/m_integralEvals << "\n"
                    << "Min. sigma: " << m_sigma_min << "\n"
                    << "Max. sigma: " << m_sigma_max << "\n"
                    << std::endl;
   delete m_acceptanceCone;
}

double PsfIntegralCache::
angularIntegral(double sigma, double gamma, size_t ipsi) const {
   m_calls++;
   if (sigma < m_sigmas.front() || sigma > m_sigmas.back() ||
       gamma < m_gammas.front() || gamma > m_gammas.back()) {
      double value = psfIntegral(m_psis.at(ipsi), sigma, gamma)/sigma/sigma;
      return value;
   }

   size_t isig(std::upper_bound(m_sigmas.begin(), m_sigmas.end(), sigma)
               - m_sigmas.begin() - 1);
   size_t igam(std::upper_bound(m_gammas.begin(), m_gammas.end(), gamma)
               - m_gammas.begin() - 1);

   size_t is[2] = {isig, isig + 1};
   size_t ig[2] = {igam, igam + 1};

   for (size_t i(0); i < 2; i++) {
      for (size_t j(0); j < 2; j++) {
         size_t indx(is[i]*m_gammas.size() + ig[j]);
         if (m_needIntegral.at(ipsi).at(indx)) {
            m_interpolations++;
            m_angularIntegral.at(ipsi).at(indx) =
               psfIntegral(m_psis.at(ipsi), m_sigmas.at(is[i]), 
                           m_gammas.at(ig[j]))
               /m_sigmas.at(is[i])/m_sigmas.at(is[i]);
            m_needIntegral.at(ipsi).at(indx) = false;
         }
      }
   }

   double value = bilinear(sigma, gamma, ipsi, isig, igam);
   return value;
}

double PsfIntegralCache::bilinear(double sigma, double gamma, size_t ipsi,
                                  size_t isig, size_t igam) const {
   double tt = ( (gamma - m_gammas.at(igam))
                 /(m_gammas.at(igam+1) - m_gammas.at(igam)) );
   double uu = ( (log(sigma) - log(m_sigmas.at(isig)))
                 /(log(m_sigmas.at(isig+1)) - log(m_sigmas.at(isig))) );
   double y1 = m_angularIntegral.at(ipsi).at(isig*m_gammas.size() + igam);
   double y2 = m_angularIntegral.at(ipsi).at(isig*m_gammas.size() + igam + 1);
   double y3 = m_angularIntegral.at(ipsi).at((isig + 1)*m_gammas.size() + igam);
   double y4 = 
      m_angularIntegral.at(ipsi).at((isig + 1)*m_gammas.size() + igam + 1);
   double value = (1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
      + tt*uu*y3 + (1. - tt)*uu*y4;
   return value;
}

double PsfIntegralCache::
psfIntegral(double psi, double sigma, double gamma) const {
   std::clock_t start_time(std::clock());

   m_gamma_avg += gamma;
   m_sigma_avg += sigma;
   m_integralEvals++;

   if (m_gamma_max < gamma) {
      m_gamma_max = gamma;
   }
   if (m_gamma_min > gamma) {
      m_gamma_min = gamma;
   }

   if (m_sigma_max < sigma) {
      m_sigma_max = sigma;
   }
   if (m_sigma_min > sigma) {
      m_sigma_min = sigma;
   }

   double roi_radius(m_acceptanceCone->radius()*M_PI/180.);
   double one(1.);
   double mup(std::cos(roi_radius + psi));
   double mum(std::cos(roi_radius - psi));

//   double err(1e-5);
   double err(1e-3);
   int ierr(0);

   double firstIntegral(0);
   if ( mum < 0.9999 ) {
      PsfIntegrand1 psfIntegrand1(sigma, gamma);
      firstIntegral = 
         st_facilities::GaussianQuadrature::dgaus8(psfIntegrand1, mum,
                                                   one, err, ierr);
   }
   
   double secondIntegral(0);
   PsfIntegrand2 psfIntegrand2(sigma, gamma, psi, roi_radius);
   secondIntegral =
      st_facilities::GaussianQuadrature::dgaus8(psfIntegrand2, mup, mum, 
                                                err, ierr);
   
   m_cpuTotal += (std::clock() - start_time)/1e6;
   double value = firstIntegral + secondIntegral;
   return value;
}

PsfIntegralCache::PsfIntegrand1::PsfIntegrand1(double sigma, double gamma)
   : m_sigma(sigma), m_gamma(gamma) {}

double PsfIntegralCache::PsfIntegrand1::operator()(double mu) const {
   double r(std::acos(mu)/m_sigma);
   double u(r*r/2.);
   return 2.*M_PI*Psf::old_base_function(u, m_sigma, m_gamma);
}

PsfIntegralCache::PsfIntegrand2::
PsfIntegrand2(double sigma, double gamma, double psi, double roi_radius)
   : m_sigma(sigma), m_gamma(gamma), m_cp(std::cos(psi)), 
     m_sp(std::sin(psi)), m_cr(std::cos(roi_radius)) {}

double PsfIntegralCache::PsfIntegrand2::operator()(double mu) const {
   double r(std::acos(mu)/m_sigma);
   double u(r*r/2.);
   double phimin(0);
   double arg((m_cr - mu*m_cp)/std::sqrt(1. - mu*mu)/m_sp);
   if (arg >= 1.) {
      phimin = 0;
   } else if (arg <= -1.) {
      phimin = M_PI;
   } else {
      phimin = std::acos(arg);
   }
   return 2.*phimin*Psf::old_base_function(u, m_sigma, m_gamma);
}

void PsfIntegralCache::setupAngularIntegrals() {
   m_angularIntegral.clear();
   m_needIntegral.clear();
   size_t npts(m_gammas.size()*m_sigmas.size());
   for (size_t i(0); i < m_psis.size(); i++) {
      std::vector<double> drow(npts, 0);
      std::vector<bool> brow(npts, true);
      m_angularIntegral.push_back(drow);
      m_needIntegral.push_back(brow);
   }
}

void PsfIntegralCache::fillParamArrays() {
   size_t npsi(500);
   double psimin(1E-2*M_PI/180.);
   double psimax(180*M_PI/180.);
   logArray(psimin, psimax, npsi, m_psis);
   m_psis.insert(m_psis.begin(),0.0);

   size_t ngam_fine(100);
   size_t ngam_coarse(25);
// These upper and lower values mirror the parameter fit boundaries in
// PointSpreadFunction.cxx, which are also inside anonymous namespace
// and so are inaccessible outside of that file.
   linearArray(1, 1.2, ngam_fine, m_gammas);
   m_gammas.pop_back();
   linearArray(1.2, 5.1, ngam_fine, m_gammas, false);
//    linearArray(1.2, 5., ngam_coarse, m_gammas, false);
//    linearArray(5., 21., ngam_coarse, m_gammas, false);

   size_t nsig(100);
// Smallest angular scales expected at highest energies.
//    double sigmin(m_psf.scaleFactor(5.62e6)*0.15);
//    double sigmax(m_psf.scaleFactor(30)*2.0);
/// @todo Remove dependence on isFront for range of sigma values.  This is
/// done here in order to check consistency against handoff_response.
   bool isFront;
   double sigmin(m_psf.scaleFactor(5.62e6, isFront=true)*0.15);
   double sigmax(m_psf.scaleFactor(30, isFront=false)*2.0);
   logArray(sigmin, sigmax, nsig, m_sigmas);
}

void PsfIntegralCache::linearArray(double xmin, double xmax, size_t nx, 
                                   std::vector<double> & xx, bool clear) const {
   if (clear) {
      xx.clear();
   }
   double xstep((xmax - xmin)/(nx-1));
   for (size_t i(0); i < nx; i++) {
      xx.push_back(xstep*i + xmin);
   }
}

void PsfIntegralCache::logArray(double xmin, double xmax, size_t nx, 
                                std::vector<double> & xx, bool clear) const {
   if (clear) {
      xx.clear();
   }
   double xstep(std::log(xmax/xmin)/(nx-1));
   for (size_t i(0); i < nx; i++) {
      xx.push_back(xmin*std::exp(xstep*i));
   }
}

} // namespace latResponse
