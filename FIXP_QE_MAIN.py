import numpy as np
from iterate import iteration
from QEINTERFACE.QEMODULE import ABIIO, PHIO
import os
if __name__=='__main__':
    natom = 5;
    workingfolder = "/workspace/jiahaoz/FIXP/QECAL/";
    pwcommand = "/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 24 {0:s}pw.x".format(workingfolder);
    phcommand = "/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 48 {0:s}ph.x -nk 12 -nd 4".format(workingfolder);
    dynmatfilename ="bto.dyn"
    startconfig = ABIIO(natom, "{0:s}btoE.in".format(workingfolder), "{0:s}bto.in".format(workingfolder));
    startconfig.obtaindipolescf("{0:s}btoE.out".format(workingfolder));
    startconfig.obtainforce("{0:s}btoE.out".format(workingfolder));
    startconfig.obtainefield("{0:s}btoE.in".format(workingfolder));
    startPH = PHIO(natom);
    startPH.obtainph("{0:s}ph.out".format(workingfolder));
    startPH.obtaindyn("{0:s}{1:s}".format(workingfolder, dynmatfilename));
    DeltaP = np.array([0.0, 0.0, 0.05]);
    Ptarget = DeltaP + startconfig.polarization;
    iter = iteration(natom, Ptarget, 0.00001, 0.00001);
    ITERMAXTIMES = 10;
    os.chdir(workingfolder);
    iter.loading(startconfig.efield, startconfig.atomp, startconfig.axis, startconfig.polarization);
    for i in range(ITERMAXTIMES):
        iter.nextsolution(startconfig.force, startPH.dynmatrix, startPH.atomcharge, startPH.edie, startconfig.polarization, i);
        startconfig.writescfE(iter.Echain[-1], iter.Pchain[-1], "{0:s}ITER{1:d}".format(workingfolder, i));
        startconfig.writenewscfnoe(iter.Pchain[-1], "{0:s}ITERNOE{1:d}".format(workingfolder, i));
        os.system("{0:s} < {2:s}ITER{1:d} > {2:s}ITER.out{1:d}".format(pwcommand, i, workingfolder));
        os.system("{0:s} < {2:s}ITERNOE{1:d} > {2:s}ITERNOE.out{1:d}".format(pwcommand, i, workingfolder));
        os.system("{0:s} < {2:s}ph.in > {2:s}PHOUT{1:d}".format(phcommand, i, workingfolder));
        os.system("cp {0:s}{2:s} {0:s}DYNOUT{1:d}".format(workingfolder, i, dynmatfilename));
        nextconfig = ABIIO(natom, "{0:s}ITER{1:d}".format(workingfolder, i), "{0:s}ITERNOE{1:d}".format(workingfolder, i));
        nextconfig.obtaindipolescf("{0:s}ITER.out{1:d}".format(workingfolder, i));
        nextconfig.obtainforce("{0:s}ITER.out{1:d}".format(workingfolder, i));
        nextconfig.obtainefield("{0:s}ITER{1:d}".format(workingfolder, i));
        nextPH = PHIO(natom);
        nextPH.obtainph("{0:s}PHOUT{1:d}".format(workingfolder, i));
        nextPH.obtaindyn("{0:s}DYNOUT{1:d}".format(workingfolder, i));
        nextconfig.offsetbyperiod(startconfig, nextPH)
        print(startconfig.polarization, nextconfig.polarization);
        startconfig = nextconfig;
        startPH = nextPH;
