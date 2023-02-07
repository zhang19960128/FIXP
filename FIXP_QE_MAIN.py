import numpy as np
from iterate import iteration
from QEINTERFACE.QEMODULE import ABIIO, PHIO
import os
if __name__=='__main__':
    natom = 5;
    workingfolder = './QECAL/';
    pwcommand = "/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 48 {0:s}pw.x".format(workingfolder);
    phcommand = "/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 48 {0:s}ph.x".format(workingfolder);
    dynmatfilename ="bto.dyn"
    startconfig = ABIIO(natom, workingfolder + 'btoE.in', workingfolder + 'bto.in');
    startconfig.obtaindipolescf(workingfolder + "btoE.out");
    startconfig.obtainforce(workingfolder + "btoE.out");
    startconfig.obtainefield(workingfolder + 'btoE.in');
    startPH = PHIO(natom);
    startPH.obtainph(workingfolder + 'ph.out');
    startPH.obtaindyn(workingfolder + dynmatfilename);
    DeltaP = np.array([0.0, 0.0, 0.05]);
    Ptarget = DeltaP + startconfig.polarization;
    iter = iteration(natom, Ptarget, 0.00001, 0.00001);
    ITERMAXTIMES = 8;
    for i in range(ITERMAXTIMES):
        iter.loading(startconfig.efield, startconfig.atomp, startconfig.axis, startconfig.polarization);
        iter.nextsolution(startconfig.force, startPH.dynmatrix, startPH.atomcharge, startPH.edie, startconfig.polarization);
        startconfig.writescfE(iter.Echain[-1], iter.Pchain[-1], workingfolder + "ITER{0:d}".format(i));
        startconfig.writenewscfnoe(iter.Pchain[-1], workingfolder + "ITERNOE{0:d}".format(i));
        os.system("{0:s} < {2:s}ITER{1:d} > {2:s}ITER.out{1:d}".format(pwcommand, i, workingfolder));
        os.system("{0:s} < {2:s}ITERNOE{1:d} > {2:s}ITERNOE.out{1:d}".format(pwcommand, i, workingfolder));
        os.system("{0:s} < {2:s}ph.in > {2:s}PHOUT{1:d}".format(phcommand, i, workingfolder));
        os.system("cp {0:s}{2:s} DYNOUT{1:d}".format(workingfolder, i, dynmatfilename));
        nextconfig = ABIIO(natom, workingfolder + 'ITER{0:d}'.format(i), workingfolder + "ITERNOE{0:d}".format(i));
        nextconfig.obtaindipolescf(workingfolder + "ITER{0:d}.out".format(i));
        nextconfig.obtainforce(workingfolder + "ITER{0:d}.out".format(i));
        nextconfig.obtainefield(workingfolder + 'ITER{0:d}'.format(i));
        nextPH = PHIO(natom);
        nextPH.obtainph("{0:s}PHOUT{1:d}".format(workingfolder, i));
        nextPH.obtaindyn("{0:s}DYNOUT{1:d}".format(workingfolder, i));
        print(startconfig.polarization, nextconfig.polarization);
        startconfig = nextconfig;
        startPH = nextPH;
