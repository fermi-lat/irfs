/**
 * @file Psf3.cxx
 * @brief PSF as the sum of two King model functions.
 * See http://confluence.slac.stanford.edu/x/bADIAw.
 * In contrast to Psf2, this class interpolates the distributions 
 * rather than the parameters.
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "irfUtil/IrfHdus.h"

#include "latResponse/Bilinear.h"
#include "latResponse/FitsTable.h"

#include "Psf2.h"
#include "latResponse/Psf3.h"
#include "PsfIntegralCache.h"

namespace {
   double sqr(double x) {
      return x*x;
   }
   class Array {
   public:
      Array(const std::vector<double> & values, size_t nx) 
         : m_values(values), m_nx(nx) {}
      double operator()(size_t iy, size_t ix) const {
         return m_values[iy*m_nx + ix];
      }
   private:
      const std::vector<double> & m_values;
      size_t m_nx;
   };
}

namespace latResponse {

Psf3::Psf3(const irfUtil::IrfHdus & psf_hdus, size_t iepoch, size_t nrow) 
   : PsfBase(psf_hdus("RPSF").at(iepoch).first,
             psf_hdus.convType() == 0,
             psf_hdus("RPSF").at(iepoch).second,
             psf_hdus("PSF_SCALING").at(iepoch).second),
     m_integralCache(0) {
   readFits(psf_hdus("RPSF").at(iepoch).first,
            psf_hdus("RPSF").at(iepoch).second, nrow);
   normalize_pars();
}

Psf3::Psf3(const std::string & fitsfile, bool isFront,
           const std::string & extname, size_t nrow) 
   : PsfBase(fitsfile, isFront, extname), m_integralCache(0) {
   readFits(fitsfile, extname, nrow);
   normalize_pars();
}

Psf3::Psf3(const Psf3 & other) : PsfBase(other),
                                 m_logEs(other.m_logEs),
                                 m_energies(other.m_energies),
                                 m_cosths(other.m_cosths),
                                 m_thetas(other.m_thetas),
                                 m_parVectors(other.m_parVectors),
                                 m_integralCache(0) {}

Psf3 & Psf3::operator=(const Psf3 & rhs) {
   if (this != &rhs) {
      PsfBase::operator=(rhs);
      m_logEs = rhs.m_logEs;
      m_energies = rhs.m_energies;
      m_cosths = rhs.m_cosths;
      m_thetas = rhs.m_thetas;
      m_parVectors = rhs.m_parVectors;
      m_integralCache = 0;
   }
   return *this;
}

Psf3::~Psf3() {
   delete m_integralCache;
}

double Psf3::value(const astro::SkyDir & appDir, 
                    double energy, 
                    const astro::SkyDir & srcDir, 
                    const astro::SkyDir & scZAxis,
                    const astro::SkyDir & scXAxis, 
                    double time) const {
   (void)(scXAxis);
   double sep(appDir.difference(srcDir)*180./M_PI);
   double theta(srcDir.difference(scZAxis)*180./M_PI);
   double phi(0);
   return value(sep, energy, theta, phi, time);
}

double Psf3::value(double separation, double energy, double theta,
                    double phi, double time) const {
   (void)(phi);
   (void)(time);

   double tt, uu;
   std::array<double, 4> cornerEnergies{};
   std::array<size_t, 4> indx{};
   getCornerPars(energy, theta, tt, uu, cornerEnergies, indx);

   double sep(separation*M_PI/180.);
   std::vector<double> yvals(4);
   for (size_t i(0); i < 4; i++) {
      yvals[i] = evaluate(cornerEnergies[i], sep, &m_parVectors[indx[i]][0]);
   }

   double my_value = Bilinear::evaluate(tt, uu, &yvals[0]);
   return my_value;
}

double Psf3::angularIntegral(double energy, double theta, 
                              double phi, double radius, double time) const {
   if (energy < 120.) {
      double value = IPsf::angularIntegral(energy, theta, phi, radius, time);
      return value;
   }

   double tt, uu;
   std::array<double, 4> cornerEnergies{};
   std::array<size_t, 4> indx{};
   getCornerPars(energy, theta, tt, uu, cornerEnergies, indx);
   
   std::vector<double> yvals(4);
   for (size_t i(0); i < 4; i++) {
      yvals[i] = psf_base_integral(cornerEnergies[i], radius,
                                   &m_parVectors[indx[i]][0]);
   }
   double value = Bilinear::evaluate(tt, uu, &yvals[0]);
   return value;
}

double Psf3::psf_base_integral(double energy, double radius, 
                                const double * pars) const {
   double ncore(pars[0]);
   double ntail(pars[1]);
   double score(pars[2]*scaleFactor(energy));
   double stail(pars[3]*scaleFactor(energy));
   double gcore(pars[4]);
   double gtail(pars[5]);

   double sep = radius*M_PI/180.;
   double rc = sep/score;
   double uc = rc*rc/2.;

   double rt = sep/stail;
   double ut = rt*rt/2.;
   if (gcore < 0 || gtail < 0) {
      throw std::runtime_error("gamma < 0");
   }
   return (ncore*Psf2::psf_base_integral(uc, gcore)*2.*M_PI*::sqr(score) + 
           ntail*ncore*Psf2::psf_base_integral(ut, gtail)*2.*M_PI*::sqr(stail));
}

double Psf3::angularIntegral(double energy,
                              const astro::SkyDir & srcDir,
                              double theta, 
                              double phi, 
                              const AcceptanceConeVector_t & acceptanceCones, 
                              double time) {
   (void)(phi);
   (void)(time);

   irfInterface::AcceptanceCone & cone(*acceptanceCones.at(0));
   if (!m_integralCache || cone != m_integralCache->acceptanceCone()) {
      delete m_integralCache;
      m_integralCache = new PsfIntegralCache(*this, cone);
   }

   double tt, uu;
   std::array<double, 4> cornerEnergies{};
   std::array<size_t, 4> indx{};
   getCornerPars(energy, theta, tt, uu, cornerEnergies, indx);

   double psi(srcDir.difference(cone.center()));

   std::array<double, 4> yvals{};
   for (size_t i(0); i < 4; i++) {
      yvals[i] = angularIntegral(cornerEnergies[i], psi,
                                 m_parVectors[indx[i]]);
   }
   double value(Bilinear::evaluate(tt, uu, &yvals[0]));

   return value;
}

double Psf3::angularIntegral(double energy, double psi, 
                              const std::vector<double> & pars) {
   const std::vector<double> & psis(m_integralCache->psis());
   if (psi > psis.back()) {
      std::ostringstream message;
      message << "latResponse::Psf3::angularIntegral:\n"
              << "Error evaluating PSF integral.\n"
              << "Requested source location > " 
              << psis.back()*180/M_PI 
              << " degrees from ROI center.";
      throw std::runtime_error(message.str());
   }
   size_t ii(std::upper_bound(psis.begin(), psis.end(), psi) 
             - psis.begin() - 1);

   double ncore(pars[0]);
   double ntail(pars[1]);
   double score(pars[2]*scaleFactor(energy));
   double stail(pars[3]*scaleFactor(energy));
   double gcore(pars[4]);
   double gtail(pars[5]);

   /// Remove sigma**2 scaling imposed by pars(...).  This is put back
   /// in angularIntegral below for each grid value of sigmas.  This
   /// preserves the normalization in the bilinear interpolation by
   /// explicitly putting in the important sigma-dependence.
   double norm_core(ncore*score*score);
   double norm_tail(ntail*stail*stail);

   double y1 = 
      norm_core*m_integralCache->angularIntegral(score, gcore, ii) + 
      norm_tail*ncore*m_integralCache->angularIntegral(stail, gtail, ii);
   double y2 =
      norm_core*m_integralCache->angularIntegral(score, gcore, ii+1) + 
      norm_tail*ncore*m_integralCache->angularIntegral(stail, gtail, ii+1);

   double y = ((psi - psis[ii])/(psis[ii+1] - psis[ii])*(y2 - y1)) + y1;
   return y;
}

double Psf3::evaluate(double energy, double sep, const double * pars) const {
   double ncore(pars[0]);
   double ntail(pars[1]);
   double score(pars[2]*scaleFactor(energy));
   double stail(pars[3]*scaleFactor(energy));
   double gcore(pars[4]);
   double gtail(pars[5]);

   double rc = sep/score;
   double uc = rc*rc/2.;

   double rt = sep/stail;
   double ut = rt*rt/2.;
   return (ncore*Psf2::psf_base_function(uc, gcore) +
           ntail*ncore*Psf2::psf_base_function(ut, gtail));
}

void Psf3::readFits(const std::string & fitsfile,
                     const std::string & extname, 
                     size_t nrow) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));
   const std::vector<std::string> & validFields(table->getValidFields());

   // The first four columns *must* be "ENERG_LO", "ENERG_HI", "CTHETA_LO",
   // "CTHETA_HI", in that order.
   const char * boundsName[] = {"energ_lo", "energ_hi", 
                                "ctheta_lo", "ctheta_hi"};
   for (size_t i(0); i < 4; i++) {
      if (validFields.at(i) != boundsName[i]) {
         std::ostringstream message;
         message << "latResponse::ParTables::ParTables: "
                 << "invalid header in " << fitsfile << "  "
                 << validFields.at(i) << "  " << i;
         throw std::runtime_error(message.str());
      }
   }

   // Push boundary values onto energy and theta arrays, replicating
   // parameter values along outer boundary.

   std::vector<double> elo, ehi;
   FitsTable::getVectorData(table, "ENERG_LO", elo, nrow);
   FitsTable::getVectorData(table, "ENERG_HI", ehi, nrow);
   std::vector<double> logEs;
   for (size_t k(0); k < elo.size(); k++) {
      logEs.push_back(std::log10(std::sqrt(elo[k]*ehi[k])));
   }

   std::vector<double> mulo, muhi;
   FitsTable::getVectorData(table, "CTHETA_LO", mulo, nrow);
   FitsTable::getVectorData(table, "CTHETA_HI", muhi, nrow);
   std::vector<double> cosths;
   for (size_t i(0); i < muhi.size(); i++) {
      cosths.push_back((mulo[i] + muhi[i])/2.);
   }

   size_t par_size(elo.size()*mulo.size());

   std::vector<double> values;
   for (size_t i(4); i < validFields.size(); i++) {
      const std::string & tablename(validFields[i]);
      FitsTable::getVectorData(table, tablename, values, nrow);
      if (values.size() != par_size) {
         std::ostringstream message;
         message << "Parameter array size does not match "
                 << "expected size based on energy and costheta arrays "
                 << "for table " << tablename
                 << " in  " << fitsfile;
         throw std::runtime_error(message.str());
      }

      std::vector<double> my_values;
      generateBoundaries(logEs, cosths, values, 
                         m_logEs, m_cosths, my_values);

      if (i == 4) {
         m_parVectors.resize(m_logEs.size()*m_cosths.size(),
                             std::vector<double>());
      }
      for (size_t j(0); j < my_values.size(); j++) {
         m_parVectors[j].push_back(my_values[j]);
      }
   }
   for (size_t k(0); k < m_logEs.size(); k++) {
      m_energies.push_back(std::pow(10., m_logEs[k]));
   }
   for (size_t j(0); j < m_cosths.size(); j++) {
      m_thetas.push_back(std::acos(m_cosths[j])*180./M_PI);
   }
   if (m_parVectors[0].size() != 6) {
      std::ostringstream message;
      message << "Number of PSF parameters in "
              << fitsfile
              << " does no match the expected number of 6.";
      throw std::runtime_error(message.str());
   }
   delete table;
}

void Psf3::generateBoundaries(const std::vector<double> & x,
                               const std::vector<double> & y,
                               const std::vector<double> & values,
                               std::vector<double> & xout,
                               std::vector<double> & yout,
                               std::vector<double> & values_out, 
                               double xlo, double xhi, double ylo, double yhi) {
   xout.resize(x.size() + 2);
   std::copy(x.begin(), x.end(), xout.begin() + 1);
   xout.front() = xlo;
   xout.back() = xhi;

   yout.resize(y.size() + 2);
   std::copy(y.begin(), y.end(), yout.begin() + 1);
   yout.front() = ylo;
   yout.back() = yhi;

   Array array(values, x.size());
   values_out.push_back(array(0, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(0, i));
   }
   values_out.push_back(array(0, x.size()-1));
   for (size_t j(0); j < y.size(); j++) {
      values_out.push_back(array(j, 0));
      for (size_t i(0); i < x.size(); i++) {
         values_out.push_back(array(j, i));
      }
      values_out.push_back(array(j, x.size()-1));
   }
   values_out.push_back(array(y.size()-1, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(y.size()-1, i));
   }
   values_out.push_back(array(y.size()-1, x.size()-1));
}
                               

void Psf3::normalize_pars(double radius) {
   double phi(0);
   double time(0);
   size_t indx(0);
   for (size_t j(0); j < m_thetas.size(); j++) {
      for (size_t k(0); k < m_energies.size(); k++, indx++) {
         double energy(m_energies[k]);
         double norm;
         if (energy < 120.) {
            norm = IPsf::angularIntegral(energy, m_thetas[j], phi, 
                                         radius, time);
         } else {
            norm = psf_base_integral(energy, radius, &m_parVectors[indx][0]);
         }
         m_parVectors[indx][0] /= norm;
      }
   }
}

void Psf3::getCornerPars(double energy, double theta,
                          double & tt, double & uu,
                          std::array<double, 4> & cornerEnergies,
                          std::array<size_t, 4> & indx) const {
   double logE(std::log10(energy));
   double costh(std::cos(theta*M_PI/180.));
   int i(findIndex(m_logEs, logE));
   int j(findIndex(m_cosths, costh));

   tt = (logE - m_logEs[i-1])/(m_logEs[i] - m_logEs[i-1]);
   uu = (costh - m_cosths[j-1])/(m_cosths[j] - m_cosths[j-1]);
   cornerEnergies[0] = m_energies[i-1];
   cornerEnergies[1] = m_energies[i];
   cornerEnergies[2] = m_energies[i];
   cornerEnergies[3] = m_energies[i-1];

   size_t xsize(m_energies.size());
   indx[0] = xsize*(j-1) + (i-1);
   indx[1] = xsize*(j-1) + (i);
   indx[2] = xsize*(j) + (i);
   indx[3] = xsize*(j) + (i-1);
}

int Psf3::findIndex(const std::vector<double> & xx, double x) {
   typedef std::vector<double>::const_iterator const_iterator_t;

   const_iterator_t ix(std::upper_bound(xx.begin(), xx.end(), x));
   if (ix == xx.end() && x != xx.back()) {
      std::cout << xx.front() << "  "
                << x << "  "
                << xx.back() << std::endl;
      throw std::invalid_argument("Psf3::findIndex: x out of range");
   }
   if (x == xx.back()) {
      ix = xx.end() - 1;
   } else if (x <= xx.front()) {
      ix = xx.begin() + 1;
   }
   int i(ix - xx.begin());
   return i;
}

std::vector<double> Psf3::params(size_t indx) const { 

  if(indx >= m_parVectors[0].size())
    throw std::runtime_error("Parameter index out of range.");

  std::vector<double> vals(m_parVectors.size(),0.0);
  std::vector< std::vector<double> >::const_iterator itr0 = m_parVectors.begin();
  std::vector<double>::iterator itr1 = vals.begin();
  for( ;itr0 != m_parVectors.end() && itr1 != vals.end();
       ++itr0, ++itr1)
    *itr1 = (*itr0)[indx];

  return vals;
}

void Psf3::setParams(size_t indx, const std::vector<double>& params) {
  
  if(indx >= m_parVectors[0].size())
    throw std::runtime_error("Parameter index out of range.");
  else if(params.size() != m_parVectors.size())
    throw std::runtime_error("Wrong size for parameter array.");

  delete m_integralCache;
  m_integralCache = 0;

  std::vector<double>::const_iterator itr0 = params.begin();
  std::vector< std::vector<double> >::iterator itr1 = m_parVectors.begin();

  for( ;itr0 != params.end() && itr1 != m_parVectors.end();
       ++itr0, ++itr1)
    (*itr1)[indx] = *itr0;

  normalize_pars();
}

} // namespace latResponse
