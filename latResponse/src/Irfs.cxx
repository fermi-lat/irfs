/**
 * @file Irfs.cxx
 * @brief Implementation of Irfs container class for the latResponse
 * package.  This class defers reading of the IRF files until the
 * corresponding IRF objects (aeff, psf, edisp) are accessed.
 * @author J. Chiang
 *
 * $Header$
 */

#include "irfUtil/EventTypeMapper.h"
#include "irfUtil/IrfHdus.h"

#include "latResponse/IrfLoader.h"

#include "Irfs.h"

namespace latResponse {

Irfs::Irfs(const std::string & irfName, const std::string & eventType)
   : irfInterface::Irfs(0, 0, 0, 0), m_irfName(irfName),
     m_eventType(eventType) {
   setIrfID(irfUtil::EventTypeMapper::instance().bitPos(irfName,
                                                        eventType).first);
}

Irfs::Irfs(const Irfs & other) 
   : irfInterface::Irfs(other), m_irfName(other.m_irfName),
     m_eventType(other.m_eventType) {
}

auto Irfs::aeff() -> irfInterface::IAeff * {
   if (irfInterface::Irfs::aeff() == 0) {
      irfUtil::IrfHdus hdus(irfUtil::IrfHdus::aeff(m_irfName, m_eventType));
      setAeff(IrfLoader::aeff(hdus));
   }
   return irfInterface::Irfs::aeff();
}

auto Irfs::efficiencyFactor() const -> const irfInterface::IEfficiencyFactor * {
   if (irfInterface::Irfs::efficiencyFactor() == 0) {
      irfUtil::IrfHdus hdus(irfUtil::IrfHdus::aeff(m_irfName, m_eventType));
      irfInterface::IEfficiencyFactor * eff(IrfLoader::efficiency_factor(hdus));
      if (eff) {
         const_cast<latResponse::Irfs *>(this)->setEfficiencyFactor(eff);
         delete eff;
      }
   }
   return irfInterface::Irfs::efficiencyFactor();
}

auto Irfs::psf() -> irfInterface::IPsf * {
   if (irfInterface::Irfs::psf() == 0) {
      irfUtil::IrfHdus hdus(irfUtil::IrfHdus::psf(m_irfName, m_eventType));
      setPsf(IrfLoader::psf(hdus));
   }
   return irfInterface::Irfs::psf();
}

auto Irfs::edisp() -> irfInterface::IEdisp * {
   if (irfInterface::Irfs::edisp() == 0) {
      irfUtil::IrfHdus hdus(irfUtil::IrfHdus::edisp(m_irfName, m_eventType));
      setEdisp(IrfLoader::edisp(hdus));
   }
   return irfInterface::Irfs::edisp();
}

} // namespace latResponse
