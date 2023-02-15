import matplotlib.pyplot as plt
import matplotfont.parameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib
from QEINTERFACE.QEMODULE import ABIIO
from QEINTERFACE.QEMODULE import PHIO
import numpy as np
matplotlib.use("TkAgg")
plt.rcParams.update(matplotfont.parameters.params)
fig, ax = plt.subplots(1,3);
for i in range(3):
    ax[i].xaxis.set_minor_locator(AutoMinorLocator())
    ax[i].yaxis.set_minor_locator(AutoMinorLocator())
workingfolder = "/workspace/jiahaoz/FIXP/QECAL/";
natom = 5;
startconfig = ABIIO(natom, "{0:s}btoE.in".format(workingfolder), "{0:s}bto.in".format(workingfolder));
ITERMAXTIMES = 7;
forcenorm = [];
efieldnorm = [];
Polarization = [];
Energy = [];
startconfig.obtaindipolescf("{0:s}btoE.out".format(workingfolder));
Polarization.append(startconfig.polarization);
for i in range(ITERMAXTIMES):
    interconfig = ABIIO(natom,  "{0:s}ITER{1:d}".format(workingfolder, i), "{0:s}ITERNOE{1:d}".format(workingfolder, i))
    guidingPH = PHIO(natom);
    guidingPH.obtainph("{0:s}PHOUT{1:d}".format(workingfolder, i));
    guidingPH.obtaindyn("{0:s}DYNOUT{1:d}".format(workingfolder, i));
    interconfig.obtainforce("{0:s}ITER.out{1:d}".format(workingfolder, i));
    forcenorm.append(np.linalg.norm(interconfig.force));
    interconfig.obtainefield("{0:s}ITER{1:d}".format(workingfolder, i))
    efieldnorm.append(interconfig.efield);
    interconfig.obtaindipolescf("{0:s}ITER.out{1:d}".format(workingfolder, i));
    interconfig.offsetbyperiod(startconfig, guidingPH)
    Polarization.append(interconfig.polarization);
ax[0].plot(np.arange(ITERMAXTIMES), forcenorm, '-d');
ax[0].set_xlabel('Iteration(#)')
ax[0].set_ylabel('Force Norm')
ax[1].plot(np.arange(ITERMAXTIMES), efieldnorm, '-d');
ax[1].set_xlabel('Iteration(#)')
ax[1].set_ylabel('Efield Norm')
Polarization = np.array(Polarization);
for i in range(3):
    ax[2].plot(np.arange(ITERMAXTIMES + 1), Polarization[:, i] - np.ones(len(Polarization)) * Polarization[0, i],'-d', label= 'Dir = {0:d}'.format(i));
ax[2].set_xlabel('Iteration(#)')
ax[2].set_ylabel(r'Polarization (C/m$^2$)')
ax[2].legend();
fig.tight_layout(w_pad = 0.000001);
plt.show()
