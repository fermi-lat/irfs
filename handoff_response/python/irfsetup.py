import yaml
import os,sys
from math import *
from xml.dom import minidom
import argparse

usage = "%(prog)s [options] "
description = """This generates directory structure and configuration
files for IRF generation."""
parser = argparse.ArgumentParser(usage=usage, description=description)

parser.add_argument('--meritdir',
                    default = '/afs/slac/g/glast/groups/canda/irfs/p8_merit/P8V6/allGamma', 
                    help = 'ROOT directory for input merit files.')

parser.add_argument('--livetime_bin',default=None, 
                    help='Peform setup for livetime correction.')

args = parser.parse_args()

###### CONFIGURATION, EDIT HERE #############

#common parameters
#yaml file with class definitions
yamfile='pass8_event_selections_noalias_ssd.yaml'
#IRF definitions
xmlfile="EvtClassDefs_P8R2.xml"

if not args.livetime_bin:
    #selection of files to use
    globfile='*.root'
    #location of root files: END WITH A '/'
    meritdir=args.meritdir
else:
    meritdir=os.path.join(args.meritdir,'livetimebins')
    globfile='*%s*.root'%args.livetime_bin

#where to put the path structure
target_dir='.'
#update existing setup? if not fails if directory exists
update=True

#items: yaml file class identifier, irfs class name, variant
#variant is either None or a yaml file identifier (e.g. "P8R1_EDISP0")
def get_class_types_from_xml(xmlfile):
    xmldoc = minidom.parse(xmlfile)
    class_version =xmldoc.getElementsByTagName('EventClass')[0].attributes['version'].value
    itemlist = xmldoc.getElementsByTagName('EventMap')
    classnames=[l.attributes['name'].value for l in itemlist[0].getElementsByTagName("EventCategory") if l.attributes['name'].value!="LLE"]
    typenames=[l.attributes['name'].value for l in itemlist[1].getElementsByTagName("EventCategory")]
    return class_version,classnames,typenames

class_version,classnames,typenames = get_class_types_from_xml(xmlfile)

params=yaml.load(open("params_P8R2_V6.yaml",'r'))
generated=params['MC']
irf_version=params['irf_version']
################ END CONFIG ########################################
####################################################################

def makesetup(class_version,classname,variant,irfname,livetime_bin=None):
    c = yaml.load(open(yamfile,'r'))

    setup_string = """
from IRFdefault import *
import glob
from math import *\n
"""
    if livetime_bin is None:
        setup_string+="Prune.fileName = 'skim_%s.root'\n" % (variant)
    else:
        setup_string+="Prune.fileName = 'skim_%s_%s.root'\n" % (variant,livetime_bin)
    p8class = '%s_%s'%(class_version,classname)
    if variant == "FRONT":
        variant_cut = "Tkr1FirstLayer>5.5"
    elif variant == "BACK":
        variant_cut = "Tkr1FirstLayer<5.5"
    else:
        variant_cut = c['%s_%s'%(class_version,variant)]
    setup_string+="Prune.cuts = '(%s)&&(%s)'\n" % (variant_cut,c[p8class])
    
    setup_string+="""
Prune.branchNames = '''McEnergy  McLogEnergy
McXDir  McYDir  McZDir
Tkr1FirstLayer WP8Best* WP8CT* Evt*
FswGamState TkrNumTracks CalTrackAngle Tkr1SSDVeto'''.split()

meritDir = '%s'
meritFiles = '%s'
Data.files = sorted(glob.glob(meritDir + meritFiles))\n
""" % (meritdir,globfile)

    if livetime_bin is None:
        setup_string+="Data.generated = %s\n" % str([eval(generated[key][0]) for key in ['allE','lowE','highE','vlowE']])
        setup_string+="Data.logemin = %s\n" % str([eval(str(generated[key][1])) for key in ['allE','lowE','highE','vlowE']])
        setup_string+="Data.logemax = %s\n" % str([generated[key][2] for key in ['allE','lowE','highE','vlowE']])
    else:
        n_generated=float(params['Livetime']['ngen'][livetime_bin])
        setup_string+="Data.generated = [%g]\n"%n_generated
        setup_string+="Data.logemin = [1.25]\n" 
        setup_string+="Data.logemax = [5.75]\n" 

    setup_string+="""
Data.friends = {}
Data.tree_name = 'MeritTuple'\n
"""

    setup_string+="className='%s'\n" % (classname) 
    setup_string+="root_dir='.'\n"
#    setup_string+="selectionName='%s'" % (frontstr)

    setup_string+="""
Data.var_xdir = 'WP8BestXDir'
Data.var_ydir = 'WP8BestYDir'
Data.var_zdir = 'WP8BestZDir'
Data.var_energy = 'WP8BestEnergy'

Bins.logemin= 0.75
Bins.logemax = 6.5

Bins.edisp_energy_overlap = 0
Bins.edisp_angle_overlap = 0
    
Bins.psf_energy_overlap = 0
Bins.psf_angle_overlap = 0

"""
    setup_string+="PSF.fit_pars={'ncore':[%g,1e-6,1.],'ntail':[%g,1e-6,1.0], 'score':[%g,0.1,5.0], 'stail':[%g,0.1,5], 'gcore':[%g,1.001,20.], 'gtail':[%g,1.001,20.]}\n\n" % tuple(params['PSF']['fit_pars'][variant])

    if livetime_bin is None:
        setup_string+="parameterFile = 'par_%s.root'\n\n" % (irfname)
    else:
        setup_string+="parameterFile = 'par_%s_%s.root'\n\n" % (irfname,livetime_bin)
    setup_string+=('Edisp.scaling_pars=%s\n'%str(params['Edisp']['scaling_pars'][variant]))
    setup_string+=('PSF.scaling_pars=%s\n'%str(params['PSF']['scaling_pars'][variant]))
    setup_string+="selectionName=\'%s\'"%variant
    return setup_string

if not args.livetime_bin:
    for classname in classnames:
        for typename in typenames:
            irfname='%s_%s_%s_%s'%(class_version,classname,irf_version,typename)
            print 'Doing',irfname
            dirname=target_dir+'/'+irfname+'/'
            try:
                os.mkdir(dirname)
            except OSError as err:
                if update:
                    pass #do not care if direcotry exists
                else:
                    raise err
            f = open(dirname+'setup.py','w')
            f.write(makesetup(class_version,classname,typename,irfname))
            f.close()
else:
    livetime_bin=args.livetime_bin
    for classname in ["CLEAN"]:
        for typename in typenames:
            irfname='%s_%s_%s_%s'%(class_version,classname,irf_version,typename)
            print 'Doing',irfname
            dirname=target_dir+'/'+irfname+'/'
            try:
                os.mkdir(dirname)
            except OSError as err:
                if update:
                    pass #do not care if direcotry exists
                else:
                    raise err
            f = open(dirname+'setup_%s.py'%livetime_bin,'w')
            f.write(makesetup(class_version,classname,typename,irfname,livetime_bin))
            f.close()
    
