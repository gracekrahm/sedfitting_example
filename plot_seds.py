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
sedfile = f'snap74.galaxy{galaxy}.rtout.sed'
z=6.014

def plot_sed(infile, sedfile, infile_label, sedfile_label):
	from prospect.models import priors, sedmodel
	from astropy.cosmology import FlatLambdaCDM
	from astropy import units as u
	cosmo = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)

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
	sim_wav = sps_sim.wavelengths
	from hyperion.model import ModelOutput
	from astropy import units as u
	from astropy import constants
	m = ModelOutput(sedfile)
	wav,flux = m.get_sed(inclination=0,aperture=-1)

	#converting units and redshifting galaxy
	wav  = np.asarray(wav)*u.micron 
	wav = wav.to(u.AA)
	flux = np.asarray(flux)*u.erg/u.s
	dl = cosmo.luminosity_distance(z).to('cm')
	flux /= (4.*3.14*dl**2.)
	nu = constants.c.cgs/(wav.to(u.cm))
	nu = nu.to(u.Hz)
	flux /= nu
	flux = flux.to(u.Jy)
	flux /= 3631
	sim_wav = (sim_wav*u.AA).to(u.micron)
	wav = wav.to(u.micron)

	#plotting
	plt.loglog(10000*wav, flux, label=sedfile_label)
	plt.loglog(10000*sim_wav, spec_sim, label=infile_label)

fig, ax = plt.subplots(figsize=(10,5))
plot_sed(infile, sedfile, 'Prospector Fitted SED', 'Powderday Mock SED')
ax.set_xlabel(r'$\lambda$ [$\AA$]')
ax.set_ylabel('Flux (mJy)')
plt.xlim([10**3, 10**7])
plt.ylim([10**(-12), 10**(-8)])
ax.grid()
plt.legend()

plt.title(f'snap 74 galaxy {galaxy}')
plt.savefig(f'snap74_galaxy{galaxy}_sed_plot.png')
plt.show()
