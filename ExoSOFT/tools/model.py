#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import copy
import KMlogger
from astropy import constants as const
from six.moves import range

## import from modules in ExoSOFT ##
#from . import constants as const
from .cytools import orbit, model_input_pars
from .utils import  load_di_data, load_rv_data

log = KMlogger.getLogger('main.model',lvl=100,addFH=False)

class ExoSOFTmodel(object):
    """
    """
    def __init__(self,sd):
        
        ####################
        ## member variables
        ####################    
        #resulting fit values   
        self.chi_squared_3d = 0
        self.chi_squared_di = 0
        self.chi_squared_rv = 0
        self.prior = 0
        self.sd = sd
        ## TRACK BEST CHI SQUAREDS FOUND SO FAR IN HERE?
        ## makes more sense to change this to 'ExoSOFTresults' and name the object 'Results'??!!
        ## load in the RV and Astrometry (DI) data
        (epochs_di, rapa, rapa_err, decsa, decsa_err) = load_di_data(self.sd['di_dataFile'])
        (epochs_rv, rv, rv_err, rv_inst_num) = load_rv_data(self.sd['rv_dataFile'])
        
        ## prior functions??
        self.Params = ExoSOFTparams(self.sd['omega_offset_di'], 
             self.sd['omega_offset_rv'], self.sd['vary_tc'], self.sd['tc_equal_to'], 
             self.sd['data_mode'], self.sd['low_ecc'], self.sd['range_maxs'], self.sd['range_mins'], 
             self.sd['num_offsets'])
    
        self.Data = ExoSOFTdata(epochs_di, epochs_rv, rapa, rapa_err, decsa, decsa_err,
                 rv, rv_err, rv_inst_num,self.sd['data_mode'], self.sd['pasa'])
        
        ExoSOFTpriors = self.sd['ExoSOFTpriors']
        
        self.Priors = ExoSOFTpriors(const=const, ecc_prior=self.sd['ecc_prior'], 
             p_prior=self.sd['p_prior'], inc_prior=self.sd['inc_prior'], 
             m1_prior=self.sd['m1_prior'], m2_prior=self.sd['m2_prior'], 
             para_prior=self.sd['para_prior'], inc_min=self.sd['inc_min'],
             inc_max=self.sd['inc_max'], p_min=self.sd['p_min'], p_max=self.sd['p_max'],
             para_est=self.sd['para_est'], para_err=self.sd['para_err'], 
             m1_est=self.sd['m1_est'], m1_err=self.sd['m1_err'], m2_est=self.sd['m2_est'], 
             m2_err=self.sd['m2_err'],ecc_min=self.sd['ecc_min'],ecc_max=self.sd['ecc_max'])
    
class ExoSOFTparams(object):
    """
    
    
    
    +---+--------------------+---------------+-------------------+-------+
    |   |  Directly Varried  | Model Inputs  | Stored Parameters |       |
    +---+--------------------+---------------+-------------------+-------+
    |   |    direct_pars     | model_in_pars |   stored_pars     |       |
    +---+--------------------+---------------+-------------------+-------+
    | i |     Parameter      |   Parameter   |     Parameter     | units |
    +===+====================+===============+===================+=======+
    | 0 |Mass of Primary (m1)|      m1       |        m1         |  Msun |
    +---+--------------------+---------------+-------------------+-------+
    .
    .
    .$$ FILL THIS OUT!!!!
    
    """
    def __init__(self, omega_offset_di, omega_offset_rv, vary_tc, tc_equal_to, 
                 di_only, low_ecc, range_maxs, range_mins, num_offsets):
        # params that effect calculating the full list of params from the directly varied one
        self.omega_offset_di = omega_offset_di
        self.omega_offset_rv = omega_offset_rv
        self.vary_tc = vary_tc
        self.tc_equal_to = tc_equal_to
        self.di_only = di_only
        self.low_ecc = low_ecc
        ## max/min ranges
        self.maxs = range_maxs 
        self.mins = range_mins
        ## prep versions of all param arrays
        self.num_offsets = num_offsets
        # direct_pars: [m1,m2,parallax,long_an,e OR sqrt(e)*sin(arg_peri),to/tc,p,inc,arg_peri OR sqrt(e)*cos(arg_peri),v1,v2...]
        self.direct_pars = np.zeros((9+num_offsets),dtype=np.dtype('d'))
        # model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        self.model_in_pars = np.zeros((14),dtype=np.dtype('d'))
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        self.stored_pars = np.zeros((13+num_offsets),dtype=np.dtype('d'))
        self.offsets = np.zeros((num_offsets),dtype=np.dtype('d'))
        #check_pars: [m1, m2, parallax, long_an, e, to/tc, p, inc, arg_peri]
        self.check_pars = np.zeros((9+num_offsets),dtype=np.dtype('d'))
        
    def make_model_in(self):
        """
        Convert directly varied parameters into a comprehensive list
        of those used ans inputs to during model calculations.
        
        model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        """        
        model_input_pars(self.direct_pars, self.low_ecc, self.tc_equal_to, 
                   self.vary_tc, self.di_only, self.omega_offset_di, 
                   self.omega_offset_rv, self.model_in_pars)
        self.offsets = self.direct_pars[9:]
        #print('self.offsets = '+repr(self.offsets))
        ## Wrap periodic params into allowed ranges.  ie. long_an and arg_peri
        m_par_ints = [3,9]
        min_max_ints = [3,8]
        for i in [0,1]:
            if self.mins[min_max_ints[i]] > self.model_in_pars[m_par_ints[i]]:
                #print('par was '+str(model_input_pars[m_par_ints[i]]))
                self.model_in_pars[m_par_ints[i]]+=360.0
                #print('now '+str(model_input_pars[m_par_ints[i]]))
            elif self.model_in_pars[m_par_ints[i]] > self.maxs[min_max_ints[i]]:
                #print('par was '+str(model_input_pars[m_par_ints[i]]))
                self.model_in_pars[m_par_ints[i]]-=360.0
                #print('now '+str(model_input_pars[m_par_ints[i]]))
        #print(repr(self.model_in_pars))
    
    def stored_to_direct(self,pars):
        """ take a set of parameters matching 'stored_pars' and make the 
        directly varied versions matching 'direct_pars'.
        Note:
        direct_pars: [m1,m2,parallax,long_an,e OR sqrt(e)*sin(arg_peri),to/tc,p,inc,arg_peri OR sqrt(e)*cos(arg_peri),v1,v2...]
        stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        """
        direct_ary = np.zeros((9+self.num_offsets),dtype=np.dtype('d'))
        direct_ary[0:4] = pars[0:4]
        if self.low_ecc:
            direct_ary[4] = np.sqrt(pars[4])*np.sin(np.radians(pars[9]))
            direct_ary[8] = np.sqrt(pars[4])*np.cos(np.radians(pars[9]))
        else:
            direct_ary[4] = pars[4]
            direct_ary[8] = pars[9]
        if self.vary_tc:
            direct_ary[5] = pars[6]
        else:
            direct_ary[5] = pars[5]
        direct_ary[6:8] = pars[7:9]
        direct_ary[9:] = pars[13:]
        return direct_ary
    
    def direct_to_stored(self,pars):
        """ Take a single set of parameters in 'direct' format and return 
        the matching set in 'stored' format.
        
        direct_pars: [m1,m2,parallax,long_an,e OR sqrt(e)*sin(arg_peri),to/tc,p,inc,arg_peri OR sqrt(e)*cos(arg_peri),v1,v2...]
        stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        """
        self.direct_pars = pars
        self.make_model_in()
        self.make_stored(1.0e6)
        return copy.deepcopy(self.stored_pars)
        
    def make_stored(self,chi_squared):
        """ 
        Push values in model_in_params, offsets and the resulting 
        chi_squared_3d into an array to be stored on disk during ExoSOFT.  
        Not sure how to make this work with emcee or other tools...
        """
        # model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        self.stored_pars[0:10] = self.model_in_pars[0:10]
        self.stored_pars[10] = self.model_in_pars[12] #a_tot_au
        self.stored_pars[11] = chi_squared
        self.stored_pars[12] = self.model_in_pars[13] #K
        self.stored_pars[13:] = self.offsets[:]
         
    def check_range(self):
        """Determine if all parameters in the full list are within their 
        allowed ranges.
        Range arrays corrispond to parameters in:
        [m1, m2, parallax, long_an, e, to/tc, p, inc, arg_peri, v1,v2,...]
        """        
        debugging = False
        self.check_pars[0:5] = self.model_in_pars[0:5]
        if self.vary_tc:
            self.check_pars[5] = self.model_in_pars[6]
        else:
            self.check_pars[5] = self.model_in_pars[5]
        self.check_pars[6:9] = self.model_in_pars[7:10]
        self.check_pars[9:] = self.offsets[:]
        
        if len(self.check_pars)!=len(self.maxs)!=len(self.mins):
            print("LENGTH OF CHECK_PARAMS IS NOT EQUAL TO LENGTH OF MINS OR MAXS!!!")
        in_range = True
        for i in range(len(self.check_pars)):
            if (self.check_pars[i]>self.maxs[i]) or (self.check_pars[i]<self.mins[i]):
                in_range = False
                if debugging:
                    print("Param # "+str(i)+" out of range")
                    print(str(self.mins[i])+"!> "+str(self.check_pars[i])+" OR !< "+str(self.maxs[i]))
        return in_range

class ExoSOFTdata(object):
    """
    An object to contain all the necessary data arrays and parameters to 
    calculate matching predicted data with the model.  All member variables 
    will remain constant throughout.
    
    Notes:
    -Except for rv_inst_num array, all other arrays must be ndarrays of double 
     precision floating point numbers (dtype=np.dtype('d')).
    -Arrays, epochs_di, rapa, rapa_err, decsa, and decsa_err must all have same length.
    -Arrays, epochs_rv, rv, rv_err and rv_inst_num must all have same length.
    
    Inputs:
    rv_inst_num = ndarray of positive signed or unsigned integers, of same length
                  as epochs_rv, rv, and rv_err.
    """
    def __init__(self, epochs_di, epochs_rv, rapa, rapa_err, decsa, decsa_err,
                 rv, rv_err, rv_inst_num, data_mode, pasa=False):
        
        self.epochs_di = epochs_di
        self.epochs_rv = epochs_rv
        # x/RA/PA
        self.rapa = rapa
        self.rapa_err = rapa_err
        self.rapa_model = np.zeros((len(epochs_di)),dtype=np.dtype('d'))
        # y/Dec/SA
        self.decsa = decsa
        self.decsa_err = decsa_err
        self.decsa_model = np.zeros((len(epochs_di)),dtype=np.dtype('d'))
        # RV
        self.rv = rv
        self.rv_err = rv_err
        self.rv_model = np.zeros((len(epochs_rv)),dtype=np.dtype('d'))
        # dataset/instrument number
        self.rv_inst_num = rv_inst_num
        self.data_mode = data_mode
        self.pasa = pasa

def ln_posterior(pars, Model):
    """
    Calculates the likelihood for a given set of inputs.
    Then calculate the natural logarithm of the posterior probability.
    
    -Model is of type ExoSOFTmodel.  Currently just holds resulting fit values.
    -Data is of type ExoSOFTdata, containing all input data and params to 
    produce predicted values of matching units, and arrays for predicted values.
    -Params is of type ExoSOFTparams, an class containing functions for 
    calculating versions of the 'pars' used as model inputs, and a version 
    that would be for storing to disk when ran in ExoSOFT.
    -Priors is of type ExoSOFTpriors, containing funtions for each parameter's
    prior, a function calculate combined prior given list of params, and any 
    variables necessary for those calculations.
    
    """    
    speed_test = False#$$$$$$$$$$$$$$$$$$$$
    ## convert params from raw values
    Model.Params.direct_pars = pars
    Model.Params.make_model_in()
        
    ## Range check on proposed params, set ln_post=zero if outside ranges.
    ln_post = -np.inf
    if speed_test: #$$$$$$$$$$$$$$$$$$$$
        in_range=True#$$$$$$$$$$$$$$$$$$$$
    else:#$$$$$$$$$$$$$$$$$$$$
        in_range = Model.Params.check_range()
    if in_range:         
        ## Call Cython func to calculate orbit. ie. -> predicted x,y,rv values.
        orbit(Model.Params.model_in_pars, Model.Params.offsets, Model.Data.pasa, 
              Model.Data.data_mode, Model.Data.epochs_di, Model.Data.epochs_rv, 
              Model.Data.rv_inst_num, Model.Data.rapa_model, 
              Model.Data.decsa_model, Model.Data.rv_model)
        if speed_test==False:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            chi_sqr_rv, chi_sqr_rapa, chi_sqr_decsa = 0, 0, 0
            if (len(Model.Data.epochs_rv)>0) and (Model.Data.data_mode!='DI'):
                #print('rv diffs\n'+repr(np.sort(Model.Data.rv-Model.Data.rv_model)))
                chi_sqr_rv = np.sum((Model.Data.rv-Model.Data.rv_model)**2 / Model.Data.rv_err**2)
            if (len(Model.Data.epochs_di)>0) and (Model.Data.data_mode!='RV'):
                chi_sqr_rapa = np.sum((Model.Data.rapa-Model.Data.rapa_model)**2 / Model.Data.rapa_err**2)
                chi_sqr_decsa = np.sum((Model.Data.decsa-Model.Data.decsa_model)**2 / Model.Data.decsa_err**2)
            chi_sqr_3d = chi_sqr_rv + chi_sqr_rapa + chi_sqr_decsa
            # Remember that chisqr = -2*log(Likelihood).  OR,
            ln_lik = -0.5*chi_sqr_3d
            #print('ln_lik',ln_lik)
            ## Make version of params with chi_sqr_3d for storing during ExoSOFT
            Model.Params.make_stored(chi_sqr_3d)
            #print('stored_pars',Model.stored_pars)
            ## store the chi sqr values in model object for printing in ExoSOFT.
            #print('chi_sqr_3d',chi_sqr_3d)
            Model.chi_squared_3d = chi_sqr_3d
            Model.chi_squared_di = chi_sqr_rapa + chi_sqr_decsa
            Model.chi_squared_rv = chi_sqr_rv
            
            ## Calculate priors
            prior = Model.Priors.priors(Model.Params.model_in_pars)
            Model.prior = prior
            #print('np.log(prior)',np.log(prior))
            #print('prior ',prior)
            ## calculate lnpost
            ln_post = np.log(prior) + ln_lik
            #print('ln_post ',ln_post)
        
    return ln_post

#EOF