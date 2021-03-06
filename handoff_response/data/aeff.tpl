# $Id$
#
SIMPLE  =                    T / Dummy Created by MWRFITS v1.4a                
BITPIX  =                    8 / Dummy primary header created by MWRFITS       
NAXIS   =                    0 / No data is associated with this header        
EXTEND  =                    T / Extensions may (will!) be present             
ORIGIN  = 'LISOC'              / name of organization making this file
END

XTENSION= 'BINTABLE'           / binary table extension                        
BITPIX  =                    8 / 8-bit bytes                                   
NAXIS   =                    2 / 2-dimensional binary table                    
NAXIS1  =                15224 / width of table in bytes                       
NAXIS2  =                    1 / number of rows in table                       
PCOUNT  =                    0 / size of special data area                     
GCOUNT  =                    1 / one data group (required keyword)             
TFIELDS =                    5 / number of fields in each row                  
ORIGIN  = 'LISOC'              / name of organization making this file
DATE    =                      / file creation date (YYYY-MM-DDThh:mm:ss UT)
EXTNAME = 'EFFECTIVE AREA'     / name of this binary table extension           
TTYPE1  = 'ENERG_LO'           /                                               
TFORM1  = '54E     '           /                                               
TUNIT1  = 'MeV     '           /                                               
TTYPE2  = 'ENERG_HI'           /                                               
TFORM2  = '54E     '           /                                               
TUNIT2  = 'MeV     '           /                                               
TTYPE3  = 'CTHETA_LO'          /                                               
TFORM3  = '32E     '           /                                               
TUNIT3  = '        '           /
TTYPE4  = 'CTHETA_HI'          /                                               
TFORM4  = '32E     '           /                                               
TUNIT4  = '        '           /
TTYPE5  = 'EFFAREA'            /                                               
TFORM5  = '1728E   '           /                                               
TUNIT5  = 'm2      '           /                                               
TDIM5   = '( 54, 32)'
#TTYPE6  = 'STAT_MIN'           /                                              
#TFORM6  = '1728E   '           /                                              
#TUNIT6  = 'm2      '           /                                              
#TDIM6   = '( 54, 32)'
TELESCOP= 'GLAST   '           /                                               
INSTRUME= 'LAT     '           /                                               
DETNAM  = 'FRONT   '
HDUCLASS= 'OGIP    '           /                                               
HDUDOC  = 'CAL/GEN/92-019'     /                                               
HDUCLAS1= 'RESPONSE'           /                                               
HDUCLAS2= 'EFF_AREA'           /                                               
HDUVERS = '1.0.0   '           /                                               
EARVERSN= '1992a   '           /                                               
1CTYP5  = 'ENERGY  '           / Always use log(ENERGY) for interpolation      
2CTYP5  = 'COSTHETA'           / Off-axis angle cosine                         
CREF5   = '(ENERG_LO:ENERG_HI,CTHETA_LO:CTHETA_HI)' /                          
CSYSNAME= 'XMA_POL '           /                                               
CCLS0001= 'BCF     '           /                                               
CDTP0001= 'DATA    '           /                                               
CCNM0001= 'EFF_AREA'           /                                               
CBD10001= 'VERSION(HANDOFF)'   /
CBD20001= 'CLASS(FRONTA)'      /
CBD30001= 'ENERG(18-560000)MeV' /
CBD40001= 'CTHETA(0.2-1)'      /
CBD50001= 'PHI(0-360)deg'      /
CBD60001= 'NONE'               /
CBD70001= 'NONE'               /
CBD80001= 'NONE'               /
CBD90001= 'NONE'               /
CVSD0001= '2007-01-17'         / Dataset validity start date (UTC)
CVST0001= '00:00:00'           /                                               
CDES0001= 'GLAST LAT effective area' /                    
END

XTENSION= 'BINTABLE'           / binary table extension                        
BITPIX  =                    8 / 8-bit bytes                                   
NAXIS   =                    2 / 2-dimensional binary table                    
NAXIS1  =                15224 / width of table in bytes                       
NAXIS2  =                    1 / number of rows in table                       
PCOUNT  =                    0 / size of special data area                     
GCOUNT  =                    1 / one data group (required keyword)             
TFIELDS =                    6 / number of fields in each row                  
EXTNAME = 'PHI_DEPENDENCE'     / name of this binary table extension
ORIGIN  = 'LISOC'              / name of organization making this file
DATE    =                      / file creation date (YYYY-MM-DDThh:mm:ss UT)
TTYPE1  = 'ENERG_LO'           /                                               
TFORM1  = '60E     '           /                                               
TUNIT1  = 'MeV     '           /                                               
TTYPE2  = 'ENERG_HI'           /                                               
TFORM2  = '60E     '           /                                               
TUNIT2  = 'MeV     '           /                                               
TTYPE3  = 'CTHETA_LO'          /                                               
TFORM3  = '32E     '           /                                               
TUNIT3  = '        '           /
TTYPE4  = 'CTHETA_HI'          /                                               
TFORM4  = '32E     '           /                                               
TUNIT4  = '        '           /
TTYPE5  = 'PHIDEP0 '           /                                               
TFORM5  = '1920E   '           /                                               
TUNIT5  = '        '           /                                               
TDIM5   = '( 60, 32)'
TTYPE6  = 'PHIDEP1 '           /                                               
TFORM6  = '1920E   '           /                                               
TUNIT6  = '        '           /                                               
TDIM6   = '( 60, 32)'
TELESCOP= 'GLAST   '           /                                               
INSTRUME= 'LAT     '           /                                               
DETNAM  = 'FRONT   '
HDUCLASS= 'OGIP    '           /                                               
HDUDOC  = 'CAL/GEN/92-019'     /                                               
HDUCLAS1= 'RESPONSE'           /                                               
HDUCLAS2= 'RPSF    '           /                                               
HDUVERS = '1.0.0   '           /                                               
EARVERSN= '1992a   '           /                                               
1CTYP5  = 'ENERGY  '           / Always use log(ENERGY) for interpolation      
2CTYP5  = 'COSTHETA'           / Off-axis angle cosine                         
CREF5   = '(ENERG_LO:ENERG_HI,CTHETA_LO:CTHETA_HI)' /                          
CCLS0001= 'BCF     '           /                                               
CDTP0001= 'DATA    '           /                                               
CCNM0001= 'EFF_AREA'           /                                               
CBD10001= 'VERSION(HANDOFF)'   /
CBD20001= 'CLASS(FRONTA)'      /
CBD30001= 'ENERG(18-560000)MeV' /
CBD40001= 'CTHETA(0.2-1)'      /
CBD50001= 'PHI(0-360)deg'      /
CBD60001= 'NONE'               /
CBD70001= 'NONE'               /
CBD80001= 'NONE'               /
CBD90001= 'NONE'               /
CVSD0001= '2007-01-17'         / Dataset validity start date (UTC)
CVST0001= '00:00:00'           /                                               
CDES0001= 'GLAST LAT effective area' /                    
END

XTENSION =             BINTABLE / Binary table extension
BITPIX   =                    8 / Bits per pixel
NAXIS    =                    2 / Required value
NAXIS1   =                   48 / Number of bytes per row
NAXIS2   =                    6 / Number of rows
PCOUNT   =                    0 / Normally 0 (no varying arrays)
GCOUNT   =                    1 / Required value
TFIELDS  =                    1 / Number of columns in table
EXTNAME  =    EFFICIENCY_PARAMS / Extension name
TELESCOP =                GLAST / Telescope or mission name
INSTRUME =                  LAT / Instrument name
FILTER   =                 NULL / Instrument filter in use (if any)
HDUCLASS =                 OGIP / Format confirms to OGIP standard
HDUCLAS1 =             RESPONSE / Extension contains response data
HDUCLAS2 =            PSFPARAMS / Extension contains response matrix
HDUVERS  =                1.2.0 / Version number of the format
TTYPE1   = 'EFFICIENCY_PARS'    / Trigger rate dependent efficiency parameters
TFORM1   = 6E                   / Data format of this field
CCLS0001= 'BCF     '           /                                               
CDTP0001= 'DATA    '           /                                               
CCNM0001= 'EFF_AREA'           /                                               
CBD10001= 'VERSION(HANDOFF)'   /
CBD20001= 'CLASS(FRONTA)'      /
CBD30001= 'ENERG(18-560000)MeV' /
CBD40001= 'CTHETA(0.2-1)'      /
CBD50001= 'PHI(0-360)deg'      /
CBD60001= 'NONE'               /
CBD70001= 'NONE'               /
CBD80001= 'NONE'               /
CBD90001= 'NONE'               /
CVSD0001= '2007-01-17'         / Dataset validity start date (UTC)
CVST0001= '00:00:00'           /                                               
CDES0001= 'GLAST LAT effective area' /                    
END
