# cython: embedsignature=True
from astropy import constants as const
import numpy as np

#python setup.py build_ext --inplace

#List of available C funcs in math.h:
#https://en.wikipedia.org/wiki/C_mathematical_functions
"""
Docstring at top of cytools.pyx
"""
def mean_corr_len(x):
    m = mean_corr_len_cy(x)
    return m

cdef double mean_corr_len_cy(double [:] x):
    """
    Calculates the average correlation length
    of the input parameter's data in a boxcar style.
    The correlation length follows that described in Tegmark2004.
    ie. The step in the chain where the variance is half the total chain's variance.
    We perform this repeatedly across the chain to give a more accurate average.
    """
    cdef double var_all, half_var_all, n_cl, s_x, s_xx, v
    cdef int i_last, npts, cl_tot

    npts = x.shape[0]
    var_all = var_calc(x)
    half_var_all = var_all/2.0
    n_cl = 0
    cl_tot = 0
    i_last = 0
    s_x = 0
    s_xx = 0
    for i in range(npts):
        s_x += x[i]
        s_xx += x[i]*x[i]
        v = (s_xx/float(i-i_last+1))-(s_x/float(i-i_last+1))*(s_x/float(i-i_last+1))
        if v>half_var_all:
            cl_tot += i-i_last+1
            i_last = i+1
            n_cl += 1.0
            s_x = 0
            s_xx = 0

    return float(cl_tot)/n_cl

cdef double var_calc(double [:] x):
    """
    This will calculate the "bias-corrected sample variance"
    and uses an advanced "corrected two-pass algorithm" to reduce roundoff error,
    a modified version of function on pg 728 of Numerical Recipes 3rd.
    """
    cdef double var,ep,s,ave
    cdef int npts
    var = 0
    ep = 0
    npts = x.shape[0]
    ave = mean_calc(x)
    if npts>1:
        for i in range(npts):
            s = x[i]-ave
            ep+=s
            var+=s*s
        var = (var-(ep*ep)/npts)/(npts-1)
    return var

cdef double mean_calc(double [:] x):
    """
    Calculate the mean of a 1-D array of doubles.
    """
    cdef double sm
    cdef int i, npts
    ave = 0
    i = 0
    sm = 0
    npts = x.shape[0]
    while i<npts:
        sm+=x[i]
        i+=1
    return sm/float(npts)

cdef double ecc_anomaly(double p, double tc, double to, double ecc,
               double epoch):
    """
    Calculate the Eccentric Anomaly (E).
    Remember that in RV, there is occasionally a phase shift due to
    the Tc!=To, that doesn't exist in DI.
    So, for RV pass in both values, for DI just set both to To.

    Following this, you need to call ta_anomaly to get T.

    E returned in units of radians.
    """
    cdef double ma, multiples, e_prime, ea, pi, days_per_year
    cdef int newton_count, maxiter
    cdef bint warnings_on

    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double fabs(double _x)
        double floor(double _x)

    pi = np.pi
    days_per_year = 365.2422
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #warnings_on = True  #$$$$ for debugging, kill this after finished testing?
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ## Calc Mean Anomaly
    ma = (2.0*pi*(epoch-2.0*to+tc))/(p*days_per_year)
    #Find if over one full orbit and shift into [-2pi,2pi]
    multiples = floor(ma/(2.0*pi))
    ma -= multiples*2.0*pi
    # shift into [0,2pi]
    if ma<0:
        ma += 2.0*pi

    ## if not circular, calculate their specific values
    if ecc<=0 or ma==0 or ma==2.0*pi:
        # set E to M by default for circular orbits
        return ma
    else:
        ## Perform Newton's Method to determine E
        maxiter = 50
        ea_prime = ma + ecc*sin(ma) + ((ecc*ecc)/2.0*ma)*sin(2.0*ma)
        newton_count = 0
        for newton_count in range(maxiter):
            ea = ea_prime
            ea_prime = ea - ( (ea-ecc*sin(ea)-ma))/(1.0-ecc*cos(ea))
            if fabs(ea-ea_prime)<1e-10:
                ## Convergence
                break
            elif newton_count== maxiter -1:
                ## Failure!
                raise
        ## double check E found satisfies original equation
        #if (fabs((ea-ecc*sin(ea))-ma)>1e-5):
        #    if warnings_on:
        #        print "PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"
        #        if True:
        #            print "M = "+str(ma)
        #            print "E = "+str(ea)
        #            print "T0 = "+str(to)
        #            print "Tc = "+str(tc)
        #            print "P = "+str(p)
        #            print "Eprime = "+str(ea_prime)
        #            print "NewtonCount = "+str(newton_count)
    return ea

cdef double ta_anomaly(double ea, double ecc):
    """
    Calculate the True Anomaly from the Eccentric Anomaly.
    Thus, call ea_anomaly func before this one to get E.

    TA returned in units of radians.
    """
    cdef double ta, pi

    cdef extern from "math.h":
        double cos(double _x)
        double acos(double _x)

    pi = np.pi

    # calculate TA from E
    ta = acos( (cos(ea)-ecc)/(1.0-ecc*cos(ea)) )
    # both increase properly 0->180, but as E properly continues to
    # increase from 180->360, TA decreases.  So correcting for that here.
    if ea>pi:
        ta = (2.0*pi) - ta

    return ta

cdef void orbit_rv(double [:] epochs, double ecc, double to, double tc, double p,
              double inc, double arg_peri_rv, double k, double [:] offsets,
               int [:] rv_inst_num, double [:] rv_model):
    """
    Calculates the predicted rv for a each epoch in epochs array.

    model value = calculated rv + the instrument specific offset.
    """
    cdef double top, ea, ta, rv, arg_peri_rad, pi
    cdef int npts
    cdef extern from "math.h":
        double cos(double _x)

    pi = np.pi

    npts = epochs.shape[0]
    arg_peri_rad = arg_peri_rv*(pi/180.0)

    ## calculate rv for each epoch in the data
    for i in range(npts):
        ## Call anomalies
        ea = ecc_anomaly(p, tc, to, ecc, epochs[i])
        ta = ta_anomaly(ea, ecc)
        ## Calc predicted RV + the instrument specific offset
        rv = k*( cos(ta+arg_peri_rad) + ecc*cos(arg_peri_rad) )
        rv_model[i] = rv + offsets[rv_inst_num[i]]

cdef void orbit_di(double [:] epochs, double long_an, double ecc, double to,
              double tc, double p, double inc, double arg_peri_di,
              double a_tot_au, double parallax, double [:] rapa_model,
              double [:] decsa_model, bint pasa):
    """
    Calculates the predited x/RA/PA and y/Dec/SA for a single epoch.
    """
    cdef double ea, a_thi, b_thi, f_thi, g_thi, x_thi, y_thi, ra, dec, pa, sa, pi
    cdef double long_an_rad, arg_peri_rad, inc_rad, sep_arcsec
    cdef int npts
    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double sqrt(double _x)
        double atan2(double _x,double _y)

    pi = np.pi
    npts = epochs.shape[0]
    long_an_rad = long_an*(pi/180.0)
    arg_peri_rad = arg_peri_di*(pi/180.0)
    inc_rad = inc*(pi/180.0)
    sep_arcsec = a_tot_au*(parallax/1000.0)

    ## Calculate predicted x/RA/PA and y/Dec/SA values for each epoch in the data
    ## Use Thiele-Innes Method to calculate predicted x/RA/PA and y/Dec/SA
    a_thi = sep_arcsec * ( cos(long_an_rad)*cos(arg_peri_rad) - sin(long_an_rad)*sin(arg_peri_rad)*cos(inc_rad))
    b_thi = sep_arcsec * ( sin(long_an_rad)*cos(arg_peri_rad) + cos(long_an_rad)*sin(arg_peri_rad)*cos(inc_rad))
    f_thi = sep_arcsec * (-cos(long_an_rad)*sin(arg_peri_rad) - sin(long_an_rad)*cos(arg_peri_rad)*cos(inc_rad))
    g_thi = sep_arcsec * (-sin(long_an_rad)*sin(arg_peri_rad) + cos(long_an_rad)*cos(arg_peri_rad)*cos(inc_rad))
    for i in range(npts):
        ## Call anomalies
        ea = ecc_anomaly(p, tc, to, ecc, epochs[i])
        # The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
        x_thi = cos(ea)-ecc
        y_thi = sqrt(1.0-ecc*ecc)*sin(ea)
        ## Calculate the predicted x&y in ["], or PA[deg], SA["]
        #  KEY NOTE: x_TH-I = y_plot = North = Dec = A*X +F*Y
        #            y_TH-I = x_plot = East  = RA  = B*X +G*Y
        dec = a_thi*x_thi + f_thi*y_thi
        ra  = b_thi*x_thi + g_thi*y_thi
        ## check which and then store RA/Dec or PA/SA accordingly
        if pasa:
            # convert RA and Dec to PA and SA
            pa = atan2(ra,dec)
            if (pa<0.0):
                pa+=2.0*pi
            sa = sqrt(ra*ra+dec*dec)
            rapa_model[i] = pa*(180.0/pi)
            decsa_model[i] = sa
        else:
            rapa_model[i] = ra
            decsa_model[i] = dec

def model_input_pars(double [:] pars, bint low_ecc, bint tc_equal_to,
                     bint vary_tc, str data_mode, double omega_offset_di,
                     double omega_offset_rv, double [:] model_in_pars):
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    # model_in_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]

    cdef double m1, m2, parallax, long_an, ecc, p, inc, arg_peri, tc, pi,sec_per_year, days_per_year
    cdef double sqrte_sinomega, sqrte_cosomega
    cdef double a_tot, top, arg_peri_di, arg_peri_rv
    cdef double ta_temp, half_ea, m_t_tc, delta_t

    cdef extern from "math.h":
        double sin(double _x)
        double cos(double _x)
        double atan2(double _x, double _y)
        double sqrt(double _x)
        double pow(double _x, double _y)

    pi = np.pi
    days_per_year = 365.2422
    sec_per_year = 60*60*24*days_per_year

    ## push pars in array into meaningful variables.
    ########### special conversions for specialty parametrizations ############
    ## Convert sqrt(e)sin(omega)&sqrt(e)cos(omega) => e and omega if required.
    ecc, arg_peri = pars[4], pars[8]
    if low_ecc:
        sqrte_sinomega, sqrte_cosomega = pars[4],pars[8]
        if 0.0 not in [sqrte_sinomega, sqrte_cosomega]:
            ecc = sqrte_sinomega**2  +  sqrte_cosomega**2
            arg_peri = (180.0/pi)*atan2(sqrte_sinomega, sqrte_cosomega)
    ###########################################################################
    [m1, m2, parallax, long_an] = pars[0:4]
    [to, p, inc] = pars[5:8]

    ## update varied omega values to omegaDI and omegaRV here
    # Get the model version of each omega.
    # This is required as the RV values are for the primary and thus omega+pi
    # compared to the value used for the secondary in the DI data.
    ## Another way to think of this is arg_peri_di = arg_peri_companion
    ##                                 arg_peri_rv = arg_peri_primary
    arg_peri_di = arg_peri + omega_offset_di
    arg_peri_rv = arg_peri + omega_offset_rv +180.0

    ## Calculate a_tot and a_tot_au
    a_tot = 0.0
    if m1!=0:
        top = p*p*sec_per_year*sec_per_year *const.G.value*const.M_sun.value*(m1+m2)
        a_tot =pow( top/(4.0*pi*pi) , (1.0/3.0))
    a_tot_au = a_tot/const.au.value

    ## set K, Tc/To to defaults, then check if they need to be calculated
    #  properly for use in the RV model.
    k = 0.0
    tc = to
    if data_mode!="DI":
        ## Calculate K
        #  NOTE: both of the below versions produce identical values of K.
        #        Using 'semi-major version' by default for simplicity/speed.
        # semi-major axis version
        top =  2.0*pi * (a_tot/(1.0+(m1/m2))) * sin(inc*(pi/180.0))
        k =  top / (p*sec_per_year*sqrt(1.0-ecc*ecc))
        # masses version
        #cdef double part1, part2, part3
        #part1 = pow( (2.0*pi*const.G.value)/(p*sec_per_year) , (1.0/3.0) )
        #part2 = pow( (m2*const.M_sun.value)/((m1+m2)*const.M_sun.value) , (2.0/3.0) )
        #part3 = sin(inc*(pi/180.0)) / sqrt(1.0-ecc*ecc)
        #k = part1 * part2 * part3

        ## Calculate Tc <-> T if needed.  Note, only relevant to RV data.
        if tc_equal_to==False:
            ta_temp = (pi/2.0)-arg_peri_rv*(pi/180.0)
            half_ea = atan2( sqrt(1.0-ecc)*sin(ta_temp/2.0) , sqrt(1.0+ecc)*cos(ta_temp/2.0) )
            m_t_tc = 2.0*half_ea-ecc*sin(2.0*half_ea);
            delta_t = (m_t_tc*p*days_per_year)/(2.0*pi);
            if vary_tc:
                # If directly varying Tc, rather than To, then
                # value in the pars array is Tc, so exchange vars
                to = tc - delta_t
            else:
                tc = to + delta_t

    #model_in_pars[0:9] = m1, m2, parallax, long_an, ecc, to, tc, p, inc
    #model_in_pars[9:14] = arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k
    #### Below is a hacky way to do it, but I kept getting the error  ###$$$$$$$$$$$$
    #### due to indexing being wrong.  Fix this !!! #$$$$$
    model_in_pars[0],model_in_pars[1],model_in_pars[2] = m1, m2, parallax
    model_in_pars[3],model_in_pars[4],model_in_pars[5] = long_an, ecc, to
    model_in_pars[6],model_in_pars[7],model_in_pars[8] = tc, p, inc
    model_in_pars[9],model_in_pars[10] = arg_peri, arg_peri_di
    model_in_pars[11],model_in_pars[12],model_in_pars[13] = arg_peri_rv, a_tot_au, k

def orbit(double [:] model_in_pars, double [:] offsets, bint pasa, str data_model,
          double [:] epochs_di, double [:] epochs_rv, int [:] rv_inst_num,
          double [:] rapa_model, double [:] decsa_model, double [:] rv_model):
    """
    The version of orbit that is called from Python.
    Pass in the parameters and output data array to have loaded up in place.
    This will call fast cdef functions to do the actual calculations.
    """
    cdef double m1, m2, parallax, long_an, ecc, to, tc, p, inc
    cdef double arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k
    cdef int npts_di, npts_rv

    npts_di = epochs_di.shape[0]
    npts_rv = epochs_rv.shape[0]

    ## push model input pars in array into meaningful variables.
    [m1, m2, parallax, long_an, ecc, to, tc, p, inc] = model_in_pars[0:9]
    [arg_peri, arg_peri_di, arg_peri_rv, a_tot_au, k] = model_in_pars[9:14]

    ## run through each epoch and calc predicted x,y,rv as requested
    if (npts_rv>0) and (data_model!='DI'):
        ## calc predicted rv if necessary
        orbit_rv(epochs_rv, ecc, to, tc, p, inc, arg_peri_rv, k, offsets,
                  rv_inst_num,rv_model)

    if (npts_di>0) and (data_model!='RV'):
        ## calc predicted DI if necessary
        orbit_di(epochs_di, long_an, ecc, to, tc, p, inc, arg_peri_di,
                 a_tot_au, parallax, rapa_model, decsa_model, pasa)

#EOF
