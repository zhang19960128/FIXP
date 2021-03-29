import numpy as np
import sys
def obtainnatoms(filedftinput):
  f=open(filedftinput,'r');
  lines=f.readlines();
  for i in lines:
    if i.find("nat")!=-1:
      arra=i.split('=');
      arra=arra[1].split(',');
      natom=int(arra[0]);
  f.close();
  return natom;
def readposition(scfin,natoms):
  scfinput=open(scfin,'r');
  lines=scfinput.readlines();
  length=len(lines);
  atomp=np.zeros([natoms,3]);
  for i in range(length):
    if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("(angstrom)")!=-1:
      for j in range(natoms):
        for k in range(3):
          atomp[j,k]=float(lines[i+j+1].split()[k+1]);
    if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("crystal")!=-1:
      for j in range(natoms):
        for k in range(3):
          atom[j,k]=float(lines[i+j+1].split()[k+1]);
  return atomp;
def writenewscf(natoms,scfin,atomposition,filename,unit):
    scffiles=open(scfin,'r');
    newfilename=open(filename,'w');
    lines=scffiles.readlines();
    tick1=0;
    tick2=0;
    skiplines=0;
    for i in range(len(lines)):
      if "ATOMIC_POSITIONS".lower() in lines[i].lower():
        newfilename.write("ATOMIC_POSITIONS ("+unit+")\n");
        tick2=i;
        for j in range(natoms):
          temp="";
          for t in range(3):
            temp=temp+" "+'{:12.8f}'.format(atomposition[j,t]);
          newfilename.write(lines[j+i+1].split()[0]+" "+temp+"\n");
        skiplines=natoms+skiplines+1;
      if skiplines > 0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
    newfilename.close();
    scffiles.close();
def updatedft(outputpath,dftinput,dftoutput,unit):
    dftainput=outputpath+'DFTAFTER.in';
    natoms=obtainnatoms(dftinput);
    filere=open(dftoutput,'r');
    lines=filere.readlines();
    Endsignal=0;
    atomsp=np.zeros((natoms,3))
    length=len(lines);
    for i in range(length):
      if lines[i].find("End final")!=-1:
        Endsignal=1;
      if lines[i].find("ATOMIC_POSITIONS")!=-1:
        for j in range(natoms):
          line=lines[i+j+1].split();
          for k in range(3):
            atomsp[j][k]=float(line[k+1]);
    writenewscf(natoms,dftinput,atomsp,dftainput,unit)
def switcryandcar(dftinput):
    dftainput='DFTAFTER.in'
    natoms=obtainnatoms(dftinput);
    filescfin=open(dftinput,'r');
    lines=filescfin.readlines();
    length=len(lines);
    axis=np.zeros((3,3));
    atomsp=np.zeros((natoms,3));
    atomtrans=np.zeros((natoms,3));
    for i in range(length):
      if "CELL_PARAMETERS".lower() in lines[i].lower():
        for j in range(3):
          line=lines[i+j+1].split();
          for k in range(3):
            axis[j][k]=float(line[k]);
      if "ATOMIC_POSITIONS".lower() in lines[i].lower() and "angstrom" in lines[i].lower():
        for j in range(natoms):
          line=lines[i+j+1].split();
          for k in range(3):
            atomsp[j][k]=float(line[k+1]);
        transmatrix=np.matrix.transpose(np.linalg.inv(axis));
        for j in range(natoms):
          atomtrans[j][:]=np.matmul(transmatrix,atomsp[j][:]);
        writenewscf(natoms,dftinput,atomtrans,dftainput,'crystal');
      if "ATOMIC_POSITIONS".lower() in lines[i].lower() and "crystal" in lines[i].lower():
        for j in range(natoms):
          line=lines[i+j+1].split();
          for k in range(3):
            atomsp[j][k]=float(line[k+1]);
        transmatrix=np.matrix.transpose(axis);
        for j in range(natoms):
          atomtrans[j][:]=np.matmul(transmatrix,atomsp[j][:]);
        writenewscf(natoms,dftinput,atomtrans,dftainput,'angstrom');
if __name__=='main':
  updatedft('/workspace/jiahaoz/BiFeO3/EfieldEnergy/05_phase/','bto.in','bto.out','angstrom');
