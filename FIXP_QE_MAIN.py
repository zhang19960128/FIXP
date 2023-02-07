import numpy as np
from iterate import iteration
from QEINTERFACE.QEMODULE import ABIIO, PHIO
if __name__=='__main__':
    natom = 5;
    workingfolder = './QECAL/'
    startconfig = ABIIO(natom, workingfolder+'bto.in', workingfolder+'btoE.in');
    iter = iteration(natom, 0.04, 0.00001, 0.00001);
    Spolar = startconfig.obtaindipolescf(workingfolder+"btoE.out");
    startconfig.obtainefield(workingfolder+'btoE.in');
    DeltaP = np.array([0.0, 0,0, 0.02]);
    Ptarget = DeltaP + Spolar;
    iter.loading(natom, startconfig.efield, startconfig.atomp, startconfig.axis, Spolar);