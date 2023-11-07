import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Tkagg')

import sys
def get_best(res, **kwargs):
    """Get the maximum a posteriori parameters.                                                                                                               
    From prospect.utils.plotting                                                                                                                              
    """
    imax = np.argmax(res['lnprobability'])
    theta_best = res['chain'][imax, :].copy()

    return theta_best

galaxy = int(sys.argv[1])

infile = f'galaxy{galaxy}.h5'
z=6.014

def plot_sed(infile, infile_label, sedfile_label):
	from astropy import units as u

	print('importing model and sps')
	sys.path.append('./')
	from run_prosp import build_model, build_sps	
	import prospect.io.read_results as reader
	filename = infile
	results, obs, model = reader.results_from(filename)

#	mod_sim = model
	mod_sim = build_model(z, 0.7)
	res = results
	try:
		sps_sim = reader.get_sps(res)
	except:
		sps_sim = build_sps()
	thetas_sim = get_best(res)



	spec_sim, _, _ = mod_sim.mean_model(thetas_sim, res['obs'], sps_sim)
	sim_wav = sps_sim.wavelengths*(1+z)

	obs = res['obs']
	maggies = obs['pd_sed']
	redshifted_wav = obs['pd_wav']



	#plotting
	plt.loglog(sim_wav, spec_sim, label=infile_label)
	plt.loglog(redshifted_wav, maggies, label=sedfile_label)

fig, ax = plt.subplots(figsize=(10,5))
plot_sed(infile, 'Prospector Fitted SED', 'Powderday Mock SED')
ax.set_xlabel(r'$\lambda$ [$\AA$]')
ax.set_ylabel('Flux (mJy)')
plt.xlim([10**3.7, 10**7])
plt.ylim([10**(-12), 10**(-8)])
ax.grid()
plt.legend()

plt.title(f'snap 74 galaxy {galaxy}')
#plt.show()
plt.savefig(f'snap74_galaxy{galaxy}_sed.png')
