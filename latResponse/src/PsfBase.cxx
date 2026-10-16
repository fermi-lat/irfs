/**
 * @file PsfBase.cxx
 * @brief Base class for latResponse Psfs
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <algorithm>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "latResponse/FitsTable.h"

#include "latResponse/PsfBase.h"

namespace {
   auto sqr(double x) -> double { 
      return x*x;
   }
}

namespace latResponse {

PsfBase::PsfBase(const std::string & fitsfile, bool isFront,
                 const std::string & extname,
                 const std::string & scaling_extname) {
   readScaling(fitsfile, isFront, scaling_extname);
}

PsfBase::PsfBase(const PsfBase & other) : irfInterface::IPsf(other),
                                          m_par0(other.m_par0),
                                          m_par1(other.m_par1),
                                          m_index(other.m_index),
                                          m_psf_pars(other.m_psf_pars) {}

auto PsfBase::operator=(const PsfBase & rhs) -> PsfBase & {
   if (this != &rhs) {
      irfInterface::IPsf::operator=(rhs);
      m_par0 = rhs.m_par0;
      m_par1 = rhs.m_par1;
      m_index = rhs.m_index;
      m_psf_pars = rhs.m_psf_pars;
   }
   return *this;
}

auto PsfBase::scaleFactor(double energy) const -> double {
   double tt(std::pow(energy/100., m_index));
   return std::sqrt(::sqr(m_par0*tt) + ::sqr(m_par1));
}

auto PsfBase::scaleFactor(double energy, bool isFront) const -> double {
   /// For Pass 8 event_types, only one set of scaling parameters is
   /// passed in array of size 3, so the isFront parameter is not
   /// used.
   (void)(isFront);
   double par0, par1;
   par0 = m_psf_pars.at(0);
   par1 = m_psf_pars.at(1);
   double tt(std::pow(energy/100., m_index));
   return std::sqrt(::sqr(par0*tt) + ::sqr(par1));
}

void PsfBase::readScaling(const std::string & fitsfile, bool isFront, 
                          const std::string & extname) {
   /// For Pass 8 event_types, only one set of scaling parameters is
   /// passed in array of size 3, so the isFront parameter is not
   /// used.
   (void)(isFront);

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));

   std::vector<double> values;

   FitsTable::getVectorData(table, "PSFSCALE", values);
   
   m_par0 = values.at(0);
   m_par1 = values.at(1);
   m_index = values.at(2);

   m_psf_pars.resize(values.size());
   std::copy(values.begin(), values.end(), m_psf_pars.begin());

   delete table;
}

} // namespace latResponse
