/**
 * @file IrfsFactory.cxx
 * @brief Provide pointers to a standard set of Irfs.
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#define ST_DLL_EXPORTS
#include "irfInterface/IrfsFactory.h"
#undef ST_DLL_EXPORTS

namespace irfInterface {

IrfsFactory * IrfsFactory::s_instance(0);

auto IrfsFactory::instance() -> IrfsFactory * {
   if (s_instance == 0) {
      s_instance = new IrfsFactory();
   }
   return s_instance;
}

void IrfsFactory::delete_instance() {
   delete s_instance;
   s_instance = 0;
}

IrfsFactory::~IrfsFactory() {
   std::map<std::string, Irfs *>::iterator itor = m_prototypes.begin();
   for ( ; itor != m_prototypes.end(); ++itor) {
      delete itor->second;
   }
   m_irfNames.clear();
}

auto IrfsFactory::create(const std::string & name) const -> Irfs * {
   std::map<std::string, Irfs *>::const_iterator itor 
      = m_prototypes.find(name);
   if (itor == m_prototypes.end()) {
      std::string message("irfInterface::IrfsFactory::create: ");
      message += "Cannot create Irfs object named " + name + ".\n";
      message += "Valid names are\n";
      for (itor = m_prototypes.begin(); itor != m_prototypes.end(); ++itor) {
         message += itor->first + "\n";
      }
      throw std::invalid_argument(message);
   } 
   return itor->second->clone();
}

void IrfsFactory::addIrfs(const std::string & name, Irfs * irfs, 
                          bool verbose) {
   if (m_prototypes.count(name) && verbose) {
      std::cerr << "irfInterface::IrfsFactory::addIrfs: "
                << "An Irfs object named " + name 
                << " already exists and is being replaced."
                << std::endl;
      delete m_prototypes[name];
   }
   m_prototypes[name] = irfs;
   if (!std::count(m_irfNames.begin(), m_irfNames.end(), name)) {
      m_irfNames.push_back(name);
   }
}

void IrfsFactory::getIrfsNames(std::vector<std::string> & names) const {
   names = m_irfNames;
//    names.clear();
//    std::map<std::string, Irfs *>::const_iterator itor = m_prototypes.begin();
//    for (; itor != m_prototypes.end(); ++itor) {
//       names.push_back(itor->first);
//    }
}

} // namespace irfInterface
