from QEINTERFACE.QEMODULE import ABIIO
from QEINTERFACE.QEMODULE import PHIO
import numpy as np
natom = 5;
startconfig = ABIIO(natom, './QECAL/bto.in', './QECAL/btoE3.in');
startconfig.readscf();
startconfig.obtainforce('./QECAL/btoE.out3')
startconfig.obtaindipolescf('./QECAL/btoE.out3')
startconfig.obtaindipolescf('./QECAL/btoE.out4')
startconfig.obtaindipolescf('./QECAL/btoE.out6')
PH = PHIO(natom);
PH.obtainph('QECAL/ph.out')
PH.obtaindyn('QECAL/bto.dyn')
print(PH.edie)