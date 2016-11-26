from __future__ import absolute_import
import numpy as np
import copy
import os
import timeit
from six.moves import range
import tools



# get settings and instantiate model object
const = tools.constants
s2 = "/mnt/HOME/MEGA/Dropbox/EclipseWorkspaceDB/ExoSOFT2/examples/settings.py"
sd2 = tools.load_settings_dict(s2)
Model = tools.ExoSOFTmodel(sd2)

#2452511.427 -0.2728288668 0.01020893362 0.01904839453 0.01020893362
#2452511.427 9.096162642 0.3562644929

if False:
		Model.Data.epochs_di = np.array([2452511.427],dtype=np.dtype('d'))
		Model.Data.rapa = np.array([-0.2728288668],dtype=np.dtype('d'))
		Model.Data.rapa_err = np.array([0.01020893362],dtype=np.dtype('d'))
		Model.Data.decsa = np.array([0.01904839453],dtype=np.dtype('d'))
		Model.Data.decsa_err = np.array([0.01020893362],dtype=np.dtype('d'))

		Model.Data.epochs_rv = np.array([2452511.427],dtype=np.dtype('d'))
		Model.Data.rv = np.array([9.096162642],dtype=np.dtype('d'))
		Model.Data.rv_err = np.array([0.3562644929],dtype=np.dtype('d'))

		Model.Data.rapa_model = np.ones((1),dtype=np.dtype('d'))
		Model.Data.decsa_model = np.ones((1),dtype=np.dtype('d'))
		Model.Data.rv_model = np.ones((1),dtype=np.dtype('d'))
		Model.Data.rv_inst_num = np.zeros((1),dtype=np.dtype('i'))


##############################################################################################
## define a set of starting parameters
m2 = const.kg_per_mjup/const.kg_per_msun
orbParams = [1.0,m2,50,100.6,0.048,2450639.0,2450639.0,11.9,45.0,14.8,0,0.0,0,0.0]
#ExoSOFT2 format8.81
params_direct = copy.deepcopy(Model.Params.stored_to_direct(orbParams))
#print(repr(params_direct))
#Model.Params.direct_pars = params_direct
#Model.Params.make_model_in()
#print(repr(Model.Params.model_in_pars))
params_direct = np.array(params_direct,dtype=np.dtype('d'))
#Run in ExoSOFT2
ln_post = tools.ln_posterior(params_direct, Model)

t_tot = 0.0
ntrials = 100000
for _ in range(ntrials):
	tic=timeit.default_timer()
	ln_post = tools.ln_posterior(params_direct, Model)
	t_tot+=timeit.default_timer()-tic
	#print("this round took:",tools.timeStrMaker(timeit.default_timer()-tic))
print(str(ntrials)+" rounds took: "+str(t_tot))
print("Thus, an average of "+str(float(t_tot)/float(ntrials))+" seconds per round")
#print(str(ntrials)+" rounds had an average time of:",tools.timeStrMaker(t_tot/ntrials))

#print('stored_pars outside:\n'+repr(copy.deepcopy(Model.Params.stored_pars)))


if False:
	print('rapa_model:\n'+repr(np.sort(copy.deepcopy(Model.Data.rapa_model))))
	print('decsa_model:\n'+repr(np.sort(copy.deepcopy(Model.Data.decsa_model))))
	print('rv_model:\n'+repr(np.sort(copy.deepcopy(Model.Data.rv_model))))

if True:
		## print a few basic results of the fit
    print('chi_squared_3d '+str(Model.chi_squared_3d))
    print('reduced chi_squared_3d '+str(Model.chi_squared_3d/35))
    print('reduced chi_squared_3d '+str(Model.chi_squared_di/14))
    print('reduced chi_squared_3d '+str(Model.chi_squared_rv/19))
    print('prior '+str(Model.prior))
    print('ln_post '+str(ln_post))
    
