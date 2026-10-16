/**
 * @file test.cxx
 * @brief Test program for irfInterface.
 * @author J. Chiang
 *
 * $Header$
 */

#if defined(TRAP_FPE) || defined(HEADAS)
#include <fenv.h>
#endif

#include <algorithm>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "facilities/commonUtilities.h"

#include "astro/SkyDir.h"

#include "irfInterface/AcceptanceCone.h"
#include "irfInterface/IrfsFactory.h"

#include "Aeff.h"
#include "Psf.h"
#include "Edisp.h"
#include "MyIrfLoader.h"

using namespace irfInterface;

class irfInterfaceTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(irfInterfaceTests);

   CPPUNIT_TEST(test_addIrfs);
   CPPUNIT_TEST(test_create);
   CPPUNIT_TEST_EXCEPTION(test_creation_failure, std::invalid_argument);
   CPPUNIT_TEST(test_getIrfsNames);
   CPPUNIT_TEST(psf_normalization);
   CPPUNIT_TEST(psf_integral);
   CPPUNIT_TEST(edisp_normalization);
   CPPUNIT_TEST(test_IrfRegistry);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void test_addIrfs();
   void test_create();
   void test_creation_failure();
   void test_getIrfsNames();
   void psf_normalization();
   void psf_integral();
   void edisp_normalization();
   void test_IrfRegistry();

private:

   std::map<std::string, Irfs *> m_irfs;
   std::vector<std::string> m_irfNames;
   std::vector<std::string>::iterator m_name;
   IrfsFactory * m_irfsFactory;

};

void irfInterfaceTests::setUp() {
   m_irfNames.push_back("Moe");
   m_irfNames.push_back("Larry");
   m_irfNames.push_back("Curly");
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      m_irfs[m_irfNames[i]] = new Irfs(new Aeff(), new Psf(), new Edisp(), i);
   }
   m_irfsFactory = IrfsFactory::instance();
}

void irfInterfaceTests::tearDown() {
   IrfsFactory::delete_instance();
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      delete m_irfs[m_irfNames[i]];
   }
   m_irfNames.clear();
}

void irfInterfaceTests::test_addIrfs() {
   for (m_name = m_irfNames.begin(); m_name != m_irfNames.end(); ++m_name) {
      m_irfsFactory->addIrfs(*m_name, m_irfs[*m_name]->clone());
   }
}

void irfInterfaceTests::test_create() {
   test_addIrfs();
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      Irfs * my_irfs = m_irfsFactory->create(m_irfNames[i]);
      CPPUNIT_ASSERT(my_irfs->irfID() == static_cast<int>(i));
   }
}

void irfInterfaceTests::test_creation_failure() {
   test_addIrfs();
   Irfs * my_irfs = m_irfsFactory->create("Shemp");
   (void)(my_irfs);
}

void irfInterfaceTests::test_getIrfsNames() {
   test_addIrfs();
   std::vector<std::string> names;
   m_irfsFactory->getIrfsNames(names);
   for (unsigned int i = 0; i < names.size(); i++) {
      m_name = std::find(m_irfNames.begin(), m_irfNames.end(), names[i]);
      CPPUNIT_ASSERT(m_name != m_irfNames.end());
   }
}

void irfInterfaceTests::psf_normalization() {
   double energy(100);
   double theta(0);
   double phi(0);
   double maxSep(10);

   Psf psf(maxSep);
   double tol(1e-4);

   double integral(psf.angularIntegral(energy, theta, phi, maxSep));
   CPPUNIT_ASSERT(std::fabs(integral - 1.) < tol);

   double radius(maxSep/2);
   integral = psf.angularIntegral(100., 0., 0., radius);
   double value(2.*M_PI*(1. - std::cos(radius*M_PI/180.))
                *psf.value(radius, energy, theta, phi));
   CPPUNIT_ASSERT(std::fabs(integral - value) < tol);
}

void irfInterfaceTests::psf_integral() {
   double energy(1000);
   double theta(0);
   double phi(0);

   double tol(1e-4);

   astro::SkyDir srcDir(0, 0);
   astro::SkyDir roiCenter(10, 0);
   double theta_roi(15);
   AcceptanceCone roiCone(roiCenter, theta_roi);
   std::vector<AcceptanceCone *> cones;
   cones.push_back(&roiCone);
   double maxSep(theta_roi - roiCenter.difference(srcDir)*180./M_PI);

   Psf psf(maxSep);
   double integral(psf.angularIntegral(energy, srcDir, theta, phi, cones));

// @todo replace this with a stringent test
   Psf psf2(2.*maxSep);
   integral = psf2.angularIntegral(energy, srcDir, theta, phi, cones);
   
   CPPUNIT_ASSERT(std::fabs(integral - 1.) != tol);
}

void irfInterfaceTests::edisp_normalization() {
   double theta(0);
   double phi(0);
   double time(0);

   Edisp edisp;

   double e0(100.);
   double emax(1000);
   double emin(0);

   double tol(1e-4);
   double integral(edisp.integral(emin, e0, e0, theta, phi, time));
   CPPUNIT_ASSERT(std::fabs(integral - 0.5) < tol);

   integral = edisp.integral(emin, emax, e0, theta, phi, time);
   CPPUNIT_ASSERT(std::fabs(integral - 1) < tol);
}
   
void irfInterfaceTests::test_IrfRegistry() {
   IrfRegistry & registry(IrfRegistry::instance());
   registry.registerLoader(new MyIrfLoader());

   registry.loadIrfs("my_classes");

   char * class_names[] = {"FrontA", "BackA", "FrontB", "BackB"};
   const std::vector<std::string> & classes(registry["my_classes"]);
   for (size_t i(0); i < 4; i++) {
      CPPUNIT_ASSERT(classes.at(i) == class_names[i]);
   }
}

auto main() -> int {
#if defined(TRAP_FPE) || defined(HEADAS)
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#else
      // throw std::runtime_error("Floating point exception trapping "
      //                          "cannot be enabled for this build.");
#endif
   CppUnit::TextTestRunner runner;
   
   runner.addTest(irfInterfaceTests::suite());
    
   bool result = runner.run();
   if (result) {
      return 0;
   } else {
      return 1;
   }
}
