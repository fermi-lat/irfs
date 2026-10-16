/**
 * @file PsfScaling.cxx
 * @brief Energy scaling of angular deviation for fitting Psf.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <stdexcept>

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "PsfScaling.h"

namespace {
   auto sqr(double x) -> double {
      return x*x;
   }

   void fitsReportError(int status, std::string routine) {
      if (status == 0) {
         return;
      }
      fits_report_error(stderr, status);
      std::ostringstream message;
      message << routine << ": CFITSIO error " << status;
      throw std::runtime_error(message.str());
   }

   void readRowVector(fitsfile * fptr, const std::string & colname, 
                      int row, int nelements, std::vector<double> & data) {
      std::string routineName("::readRowVector");
      int status(0);
      int colnum(0);
      fits_get_colnum(fptr, CASEINSEN, const_cast<char *>(colname.c_str()),
                      &colnum, &status);
      fitsReportError(status, routineName);
      
      long nrows(0);
      fits_get_num_rows(fptr, &nrows, &status);
      fitsReportError(status, routineName);
      if (row >= nrows) {
         std::ostringstream message;
         message << routineName << "\n"
                 << "Request for row " << row << " from a table that "
                 << "has only " << nrows << "rows.";
         throw std::runtime_error(message.str());
      }

      int anynul(0), nulval(0);
      data.resize(nelements);
      fits_read_col(fptr, TDOUBLE, colnum, row+1, 1, nelements, &nulval,
                    &data[0], &anynul, &status);
      fitsReportError(status, routineName);
   }
}

namespace dc2Response {

PsfScaling::PsfScaling(const std::vector<double> & pars, bool useFront) 
   : m_pars(pars), m_useFront(useFront) {}

PsfScaling::PsfScaling(const std::string & psfFile) {
   const tip::Table * table =
      tip::IFileSvc::instance().readTable(psfFile, "PSF_SCALING_PARAMS");
//    tip::ConstTableRecord & row(*table->begin());
//    row["PSFSCALE"].get(m_pars);
   std::string latClass;
   table->getHeader()["LATCLASS"].get(latClass);
   if (latClass.find("FRONT") != std::string::npos) {
      m_useFront = true;
   } else {
      m_useFront = false;
   }
   delete table;

// Read in the data from the single row containing the scaling parameters.
   int status(0);
   fitsfile * fptr(0);
   fits_open_file(&fptr, psfFile.c_str(), READONLY, &status);
   fitsReportError(status, "PsfScaling::PsfScaling");
   int hdutype(0);
   fits_movabs_hdu(fptr, 3, &hdutype, &status);
   fitsReportError(status, "PsfScaling::PsfScaling");
   readRowVector(fptr, "PSFSCALE", 0, 13, m_pars);
   fits_close_file(fptr, &status);
   fitsReportError(status, "PsfScaling::PsfScaling");
}

auto PsfScaling::operator()(double McEnergy, double McZDir) const -> double {
   double t(powerLawScaling(McEnergy));
   double x(trendline(McEnergy));
   double scaleFactor(1.);
   if (m_useFront) {
      scaleFactor *= (m_pars[4]*x*x - m_pars[5]*x + m_pars[6]);
      scaleFactor *= zfactor(McZDir)*std::sqrt(::sqr(m_pars[7]*t) +
                                               ::sqr(m_pars[8]));
   } else {
      scaleFactor *= (m_pars[9]*x + m_pars[10]);
      scaleFactor *= zfactor(McZDir)*std::sqrt(::sqr(m_pars[11]*t) +
                                               ::sqr(m_pars[12]));
   }
   return scaleFactor*180./M_PI;
}

auto PsfScaling::powerLawScaling(double McEnergy) const -> double {
   return std::pow(McEnergy/100., m_pars[0]);
}

auto PsfScaling::zfactor(double McZDir) const -> double {
   return 1. + m_pars[1]*(1. - std::fabs(McZDir));
}

auto PsfScaling::trendline(double McEnergy) const -> double {
   return m_pars[2]*std::log10(McEnergy/100.) + m_pars[3];
}

} //namespace dc2Response
