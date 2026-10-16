/**
 * @file HdCaldb.cxx
 * @brief Wrapper class for HEASARC CALDB routines.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include <stdexcept>

#include "st_facilities/Environment.h"
#include "st_facilities/Util.h"

#include "Hdcal.h"

#include "irfUtil/HdCaldb.h"

namespace irfUtil {

HdCaldb::HdCaldb(const std::string & telescope, const std::string & instrument)
   : m_telescope(telescope), m_instrument(instrument), m_filenamesize(1024) {

// Check that CALDB files exist.
   st_facilities::Util::file_ok(
      st_facilities::Environment::getEnv("CALDBALIAS"));
   st_facilities::Util::file_ok(
      st_facilities::Environment::getEnv("CALDBCONFIG"));

   for (int i = 0; i < s_maxret; i++) {
      m_filenames[i] = new char[m_filenamesize];
      m_online[i] = new char[m_filenamesize];
   }
}

HdCaldb::~HdCaldb() {
   for (int i = 0; i < s_maxret; i++) {
      delete [] m_filenames[i];
      delete [] m_online[i];
   }
}

auto 
HdCaldb::operator()(const std::string & detName, 
                    const std::string & respName,
                    const std::string & expression, 
                    const std::string & filter,
                    const std::string & date, 
                    const std::string & time) -> std::pair<std::string, int> {
   int nret;
   int nfound(0);
   int status(0);

   std::string startdate = date;
   std::string starttime = time;
// Leave stop date/time disabled.
   std::string stopdate = "-";
   std::string stoptime = "-";

   HDgtcalf(m_telescope.c_str(), m_instrument.c_str(),
            detName.c_str(), filter.c_str(), respName.c_str(),
            startdate.c_str(), starttime.c_str(), 
            stopdate.c_str(), stoptime.c_str(), expression.c_str(),
            s_maxret, m_filenamesize, m_filenames, m_extnums, 
            m_online, &nret, &nfound, &status);

   if (status != 0) {
      throw std::runtime_error("HdCaldb::operator(): HDgtcalf error.");
   }
   if (nfound == 0) {
      throw std::runtime_error("HdCaldb::operator(): no files found for "
                               + expression);
   }

   std::string filename = m_filenames[0];
   int extnum = m_extnums[0];

   return std::make_pair(filename, extnum);
}

void HdCaldb::getFiles(std::vector<std::string> & files,
                       std::vector<int> & extnums,
                       const std::string & detName,
                       const std::string & respName,
                       const std::string & irfName) {
   int nret;
   int nfound(0);
   int status(0);

// Disable boundary dates by setting start/stop strings to be "-".
// This will return all IRF files regardless of validity date.
   std::string startdate = "-";
   std::string starttime = "-";
   std::string stopdate = "-";
   std::string stoptime = "-";

   std::string filter("NONE");

   std::string expression("VERSION.eq." + irfName);
   HDgtcalf(m_telescope.c_str(), m_instrument.c_str(),
            detName.c_str(), filter.c_str(), respName.c_str(),
            startdate.c_str(), starttime.c_str(), 
            stopdate.c_str(), stoptime.c_str(), expression.c_str(),
            s_maxret, m_filenamesize, m_filenames, m_extnums, 
            m_online, &nret, &nfound, &status);
   if (status != 0) {
      throw std::runtime_error("HdCaldb::getFiles: HDgtcalf error.");
   }
   if (nfound == 0) {
      throw std::runtime_error("HdCaldb::getFiles: no files found for "
                               + irfName);
   }

   for (int i(0); i < nfound; i++) {
      files.push_back(m_filenames[i]);
      extnums.push_back(m_extnums[i]);
   }
}

} // namespace irfUtil
