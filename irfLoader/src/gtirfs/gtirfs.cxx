/**
 * @file gtirfs.cxx
 * @brief Apply event by event energy correction
 * @author J. Chiang
 *
 * $Header$
 */

#include <cstdlib>

#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "irfInterface/Irfs.h"

#include "irfLoader/Loader.h"

class PrintIrfs : public st_app::StApp {

public:

   PrintIrfs() : st_app::StApp(),
                     m_pars(st_app::StApp::getParGroup("gtirfs")) {
      try {
         setVersion(s_cvs_id);
      } catch (std::exception & eObj) {
         std::cerr << eObj.what() << std::endl;
         std::exit(1);
      } catch (...) {
         std::cerr << "Caught unknown exception in PrintIrfs constructor." 
                   << std::endl;
         std::exit(1);
      }
   }

   virtual ~PrintIrfs() throw() {
      try {
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }  
   }

   virtual void run();
   virtual void banner() const;

private:

   st_app::AppParGroup & m_pars;

   auto par(const std::string & key) const -> std::string;

   std::vector<irfInterface::Irfs *> m_irfs;

   static std::string s_cvs_id;

};

std::string PrintIrfs::s_cvs_id("$Name$");

st_app::StAppFactory<PrintIrfs> myAppFactory("gtirfs");

void PrintIrfs::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

auto PrintIrfs::par(const std::string & key) const -> std::string {
   std::string value = m_pars[key.c_str()];
   return value;
}

void PrintIrfs::run() {
   m_pars.Prompt();
   m_pars.Save();

   st_stream::StreamFormatter formatter("gtirfs", "", 2);

   irfLoader::Loader::go();

   typedef std::vector<std::string> IrfVector_t;
   typedef std::map<std::string, std::vector<std::string> > IrfMap_t;

   IrfMap_t respIds(irfLoader::Loader::respIds());

   for (IrfMap_t::const_iterator irf(respIds.begin());
        irf != respIds.end(); ++irf) {
      int level(2);
      if (irf->first.find("TEST") != std::string::npos ||
          irf->first.find("G25") != std::string::npos ||
          irf->first.find("DC2") != std::string::npos ||
          irf->first.find("DC1") != std::string::npos ||
          irf->first.find("HANDOFF") != std::string::npos ||
          irf->first.find("P5_v0") != std::string::npos ||
          irf->first.find("Pass4") != std::string::npos) {
         level = 3;
      }
      formatter.info(level) << irf->first;
      if (irf->second.size() > 1 || irf->second.front() != irf->first) {
         formatter.info(level) << " ( = ";
         for (IrfVector_t::const_iterator component(irf->second.begin());
              component != irf->second.end(); ++component) {
            formatter.info(level) << *component;
            if (component != irf->second.end()-1) {
               formatter.info(level) << " + ";
            }
         }
         formatter.info(level) << " )";
      }
      formatter.info(level) << std::endl;
   }
}
