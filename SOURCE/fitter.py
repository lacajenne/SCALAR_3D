from __future__ import print_function   # makes this work for python2 and 3

import sys
import collections
import gvar as gv
import numpy as np
import corrfitter as cf

T  = int(sys.argv[1])
tm = int(sys.argv[2])
t0 = int(sys.argv[3])

def main():
    data = make_data(filename='mydata')
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in range(2, 6):                                   
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0)  
        p0 = fit.pmean    
        print_results(fit)
    fastfit = cf.fastfit(G=data['cdata'], ampl='1(1)', dE='0.5(5)', tmin=tm, tp=T)
    print(fastfit)

def make_data(filename):
    """ Read data, compute averages/covariance matrix for G(t). """
    return gv.dataset.avg_data(cf.read_dataset(filename))

def make_models():
    """ Create corrfitter model for G(t). """
    corr = cf.Corr2(
        datatag='cdata', tp=T, tdata=range(T), tfit=range(t0, T-t0), 
        a='a', b='a', dE='dE'
        )
    return [corr]

def make_prior(N):
    """ Create prior for N-state fit. """
    prior = collections.OrderedDict()
    prior['a'] = gv.gvar(N * ['0(1)'])
    prior['logdE'] = gv.log(gv.gvar(N * ['0.5(5)']))                       
    return prior

def print_results(fit):
    p = fit.p                               
    E = np.cumsum(p['dE'])                              
    a = p['a']                                          
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1])) 
    print('{:2}  {:15}  {:15}\n'.format('a', a[0], a[1]))

if __name__ == '__main__':
    main()

