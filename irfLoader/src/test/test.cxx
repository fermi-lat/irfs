/**
 * @file test.cxx
 * @brief Test program for irfLoader.
 * @author J. Chiang
 *
 * $Header$
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <algorithm>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "irfInterface/IrfsFactory.h"
#include "irfInterface/IrfRegistry.h"

#include "irfLoader/Loader.h"
#include "MyLoader.h"

using irfInterface::IrfsFactory;
using namespace irfLoader;

class irfLoaderTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(irfLoaderTests);
   CPPUNIT_TEST(load_single_irfs);
   CPPUNIT_TEST_EXCEPTION(access_missing_irfs, std::invalid_argument);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void load_single_irfs();
   void access_missing_irfs();

   void test_IrfRegistry();

private:

};

void irfLoaderTests::setUp() {
// Tabula rasa
   IrfsFactory::delete_instance();
   irfInterface::IrfRegistry::instance().registerLoader(new MyLoader());
}

void irfLoaderTests::tearDown() {
}

void irfLoaderTests::load_single_irfs() {
   Loader::go();
   IrfsFactory * myFactory = IrfsFactory::instance();
   std::vector<std::string> names;
   myFactory->getIrfsNames(names);
   CPPUNIT_ASSERT(std::find(names.begin(), names.end(), 
                            std::string("my_classes::FrontA")) != names.end());
   CPPUNIT_ASSERT(std::find(names.begin(), names.end(), 
                            std::string("my_classes::BackA")) != names.end());
   CPPUNIT_ASSERT(std::find(names.begin(), names.end(), 
                            std::string("my_classes::FrontB")) != names.end());
   CPPUNIT_ASSERT(std::find(names.begin(), names.end(), 
                            std::string("my_classes::BackB")) != names.end());
}

void irfLoaderTests::access_missing_irfs() {
   Loader::go("DEV");
}
   
auto main() -> int {
   CppUnit::TextTestRunner runner;
   
   runner.addTest(irfLoaderTests::suite());
    
   bool result = runner.run();
   if (result) {
      return 0;
   } else {
      return 1;
   }
}
