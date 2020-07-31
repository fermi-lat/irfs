#
# Script to add standard set of IRF files to caldb.indx using the udcif FTOOL
#
# $Header$
#
import os
import os

#udcif = os.system('udcif')

#suffixes = ("P6_v4_trans", "P6_v4_source",
#            "P6_v5_trans", "P6_v5_source", "P6_v5_diff")
#suffixes = ("P7_v1_diff", "P7_v1_source", "P7_v1_trans", "P7_v1_S3")
#suffixes = ("P6_v8_diff", )
#suffixes = ("P6_v7_trans", "P6_v7_source", "P6_v7_diff")
#suffixes = ("p6_v3_dataclean",)
#suffixes = ("P7TRANSIENT_V6", "P7SOURCE_V6", "P7CLEAN_V6", "P7ULTRACLEAN_V6", 
#            "P7SOURCE_V6MC", 
#            "P7REP_TRANSIENT_V10", "P7REP_SOURCE_V10", 
#            "P7REP_CLEAN_V10", "P7REP_ULTRACLEAN_V10")
suffixes = ("P8R3_TRANSIENT100E_V2",
            "P8R3_TRANSIENT100_V2",
            "P8R3_TRANSIENT020E_V2",
            "P8R3_TRANSIENT020_V2",
            "P8R3_TRANSIENT010E_V2",
            "P8R3_TRANSIENT010_V2",
            "P8R3_SOURCE_V2",
            "P8R3_CLEAN_V2",
            "P8R3_ULTRACLEAN_V2",
            "P8R3_ULTRACLEANVETO_V2",
            "P8R3_SOURCEVETO_V2",
            "P8R3_TRANSIENT100S_V2",
            "P8R3_TRANSIENT015S_V2")

for item in suffixes:
    for prefix, loc in zip(("aeff", "psf", "edisp"), ("ea", "psf", "edisp")):
        os.chdir(loc)
        for section in ("front", "back", "EDISP", "FB", "PSF"):
            os.system("printenv CALDB")
            infile = "%s_%s_%s.fits" % (prefix, item, section)
            cmd = "udcif "+infile+" ../../caldb.indx editc=no quality=0"
            os.system(cmd)
            #udcif.run(infile="%s_%s_%s.fits" % (prefix, item, section),
            #          cif='../../caldb.indx', editc="no",
            #          quality=0, newentries=1)
        os.chdir('..')
