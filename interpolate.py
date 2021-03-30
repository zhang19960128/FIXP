import update
import numpy as np
from CalPatE import scanE
filebefore="DFTAFTER00.in";
fileafter="DFTAFTER01.in";
natoms=20;
beforeatoms=update.readposition(filebefore,natoms);
afteratoms=update.readposition(fileafter,natoms);
diffp=afteratoms-beforeatoms;
inter=5;
for i in range(6):
  update.writenewscf(natoms,filebefore,beforeatoms+diffp/inter*i,'INTERPOLATE'+str(i),'Angstrom')
#  print(np.linalg.norm(diffp/inter*i))
pwEscf=scanE('./PWOUT','ph.out0','dyn.out0','ite.out0','ite0');
scfoutlist=[];
scfzerolist=[];
scfinlist=[];
pathlist=[];
path="/workspace/jiahaoz/BiFeO3/EfieldEnergy/05_phase/INTERPOLATE/";
for i in range(0,6):
  ptemp=path+"{0:01d}".format(i)+'/bto.out';
  scfoutlist.append(ptemp);
  ptemp=path+"{0:01d}".format(i)+'/INTERPOLATE{0:01d}'.format(i);
  scfinlist.append(ptemp);
  ptemp=path+"{0:01d}".format(i)+'/';
  pathlist.append(ptemp);
  ptemp=path+"{0:01d}".format(i)+'/bto.out';
  scfzerolist.append(ptemp);
pwEscf.scanP(scfinlist,scfzerolist,scfzerolist);
