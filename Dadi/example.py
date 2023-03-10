#! /usr/bin/env python

import numpy
import dadi
import Selection

#the demographic function we will be using for our analysis. This describes
#a two epoch model with one historical size change.
def two_epoch(params, ns, pts):
    nu, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

#set demographic parameters and theta. this is usually inferred from
#synonymous sites
demog_params = [2, 0.05]
theta_ns = 4000.
ns = numpy.array([250])

#the demography + selection function. single size change and selection.
def two_epoch_sel(params, ns, pts):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

#integrate over a range of gammas without defining breaks
pts_l = [600, 800, 1000]
spectra = Selection.spectra(demog_params, ns, two_epoch, pts_l=pts_l, 
                            int_bounds=(1e-5,500), Npts=300, echo=True, 
                           mp=True)

#load sample data
data = dadi.Spectrum.from_file('example.sfs')

#fit a DFE to the data
#initial guess and bounds
sel_params = [0.2, 1000.]
lower_bound = [1e-3, 1e-2]
upper_bound = [1, 50000.]
p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                              upper_bound=upper_bound)
popt = Selection.optimize_log(p0, data, spectra.integrate, Selection.gamma_dist, 
                              theta_ns, lower_bound=lower_bound, 
                              upper_bound=upper_bound, verbose=len(sel_params), 
                              maxiter=30)
print(popt)
#get expected SFS for MLE
model_sfs = spectra.integrate(popt[1], Selection.gamma_dist, theta_ns)
print(model_sfs)
