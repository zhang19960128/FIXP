import update
filebefore="DFTAFTER00.in";
fileafter="DFTAFTER01.in";
natoms=20;
beforeatoms=update.readposition(filebefore,natoms);
afteratoms=update.readposition(fileafter,natoms);
diffp=afteratoms-beforeatoms;
inter=5;
for i in range(5):
  update.writenewscf(natoms,filebefore,beforeatoms+diffp/inter*i,'INTERPOLATE'+str(i),'Angstrom')
