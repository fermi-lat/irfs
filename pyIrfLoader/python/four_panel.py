"""
@brief Example "four panel" performance plots using pyIrfLoader.

@author J. Chiang
"""
#
# $Header$
#

import bisect
import numpy as num
import pylab
import pyIrfLoader

class FunctionWrapper(object):
    def __init__(self, func):
        self.func = func
    def __call__(self, xx, **kwds):
        try:
            y = []
            for x in xx:
                y.append(self.func(x, **kwds))
            if isinstance(xx, num.ndarray):
                y = num.array(y)
            return y
        except TypeError:
            return self.func(xx, **kwds)
    def __getattr__(self, attrname):
        return getattr(self.func, attrname)

_win_id = 0

class Window(object):
    def __init__(self, id=None):
        global _win_id
        if id is None:
            id = _win_id
            _win_id += 1
        self.fig = pylab.figure(id)
        self.axes = self.fig.add_subplot(111)
        self.id = id
    def set_title(self, title):
        self.axes.set_title(title)

def setAxis(xrange=None, yrange=None):
    axisrange = list(pylab.axis())
    if xrange is not None:
        axisrange[:2] = xrange
    if yrange is not None:
        axisrange[2:] = yrange
    pylab.axis(axisrange)

def plot_curve(x, y, xlog=0, ylog=0, xname='x', yname='y', 
               oplot=0, color='k', lineStyle='-', linewidth=1,
               xrange=None, yrange=None):
    if oplot == 0:
        win = Window()
    else:
        win = None
    marker = '%s%s' % (color, lineStyle)
    if xlog and ylog:
        pylab.loglog(x, y, marker, markersize=3, linewidth=linewidth)
    elif xlog:
        pylab.semilogx(x, y, marker, markersize=3, linewidth=linewidth)
    elif ylog:
        pylab.semilogy(x, y, marker, markersize=3, linewidth=linewidth)
    else:
        pylab.plot(x, y, marker, markersize=3, linewidth=linewidth)
    if not oplot:
        pylab.xlabel(xname)
        pylab.ylabel(yname)
    setAxis(xrange, yrange)
    return win

logspace = lambda xmin, xmax, nx : num.logspace(num.log10(xmin),
                                                num.log10(xmax), nx)

pyIrfLoader.Loader_go()

factory = pyIrfLoader.IrfsFactory_instance()

irfName = "P7SOURCE_V6MC"

front = factory.create(irfName + "::FRONT")
back = factory.create(irfName + "::BACK")

psf_f = front.psf()
psf_b = back.psf()

radii = logspace(1e-2, 30., 30)

@FunctionWrapper
def theta_68(energy, psf=None, inc=0, phi=0, frac=0.68):
    f = FunctionWrapper(lambda x : psf.angularIntegral(energy, inc, phi, x))
    y = f(radii)
    indx = bisect.bisect(y, frac) - 1
    return ((frac - y[indx])/(y[indx+1] - y[indx])
            *(radii[indx+1] - radii[indx]) + radii[indx])

energies = logspace(20., 3e5, 40)
plot1 = plot_curve(energies, theta_68(energies, psf=psf_f),
                   xlog=1, ylog=1, xname='Energy (MeV)',
                   yname='theta_68 (deg)')
plot1.set_title('normal incidence')
plot_curve(energies, theta_68(energies, psf=psf_b), oplot=1, lineStyle=':')

aeff_f = front.aeff()
aeff_b = back.aeff()

@FunctionWrapper
def aeff(energy, aeffObj=None, inc=0, phi=0):
    return aeffObj.value(energy, inc, phi)

plot2 = plot_curve(energies, aeff(energies, aeffObj=aeff_f), xlog=1,
                   xname='Energy (MeV)', yname='eff. area (cm^2)')
plot2.set_title('normal incidence')
plot_curve(energies, aeff(energies, aeffObj=aeff_b), oplot=1, lineStyle=':')

@FunctionWrapper
def aeff_profile(inc, aeffObj=None, energy=1e3, phi=0):
    return aeffObj.value(energy, inc, phi)

thetas = num.arange(70, dtype=float)

plot3 = plot_curve(thetas, aeff_profile(thetas, aeffObj=aeff_f),
                   xname='inclination (deg)', yname='eff. area (cm^2)')
plot3.set_title('E = 1 GeV')
plot_curve(thetas, aeff_profile(thetas, aeffObj=aeff_b), oplot=1,
           lineStyle=':')

@FunctionWrapper
def th68_profile(inc, psf=None, energy=1e3, phi=0, frac=0.68):
    f = FunctionWrapper(lambda x : psf.angularIntegral(energy, inc, phi, x))
    y = f(radii)
    indx = bisect.bisect(y, frac) - 1
    return ((frac - y[indx])/(y[indx+1] - y[indx])
            *(radii[indx+1] - radii[indx]) + radii[indx])

plot4 = plot_curve(thetas, th68_profile(thetas, psf=psf_f),
                   xname='inclination (deg)', yname='theta_68 (deg)')
plot4.set_title('E = 1 GeV')
plot_curve(thetas, th68_profile(thetas, psf=psf_b), oplot=1,
           lineStyle=':')
