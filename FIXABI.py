import numpy as np
import re
import math
import sys
import os
class abiout:
  def __init__(self,scfin,zeroout,phin,phout):
    self.phfile=phout;
    self.zerofile=zeroout;
    self.scfin=scfin;
    self.phin=phin;
    natoms=self.obtainnatoms(self.scfin);
    self.obtain(natoms);
    self.atomp=self.readposition(self.scfin,natoms);
  def obtainaxis(self,filedftinput):
    f=open(filedftinput,'r');# the units of axis is Angstrom.
    lines=f.readlines();
    self.axis=np.zeros([3,3]);
    acell=[0,0,0];
    autoA=0.52917721067121;
    for i in range(len(lines)):
      if lines[i].find("acell")!=-1:
        arra=lines[i].split();
        if lines[i].lower().find('angstrom')!=-1:
          autoA=1.0;
        elif lines[i].lower().find('bohr')!=-1:
          autoA=0.52917721067121;
        else:
          autoA=0.52917721067121;
        for j in range(3):
          acell[j]=float(arra[j+1])*autoA;
      if lines[i].find("rprim")!=-1:
        for j in range(3):
          arra=lines[i+j].split();
          for k in range(-3,0,1):
            self.axis[j][k]=float(arra[k]);
    for j in range(3):
      for k in range(3):
        self.axis[j][k]=acell[j]*self.axis[j][k];
    f.close();
  def obtainnatoms(self,filedftinput):
    f=open(filedftinput,'r');
    lines=f.readlines();
    for i in lines:
      if i.find("natom")!=-1:
        arra=i.split();
        natom=int(arra[1]);
    f.close();
    return natom;
  def obtain(self,natom):
    au=0.52917721067121;
    self.natoms=natom;
    self.atomp=np.zeros([self.natoms,3]);
    self.atommass=np.zeros([self.natoms,1]);
    self.atomcharge=np.zeros([self.natoms,3,3]);
    self.atomdis=np.arange(self.natoms*3);
    self.dynmatrix=np.zeros([self.natoms,self.natoms,3,3]);
    self.edie=np.zeros([3,3]);
    self.force=np.zeros([self.natoms,3]);
    self.polar=np.zeros(3);
    self.obtainborn(self.phfile);
    self.obtaindyn(self.phfile)
    self.obtainedie(self.phfile);
    self.obtainforce(self.zerofile);
    self.obtainaxis(self.scfin)
    self.vol=np.linalg.det(self.axis)/au/au/au;
  def obtaindyn(self,dynfile):
    dyn=open(dynfile,'r');
    lines=dyn.readlines(); # note that the dynamical matrix obtain from abinit calculation is acutally force constant matrix Units Ha/bohr/bohr
    for t in range(len(lines)):
      if lines[t].find("Dynamical matrix, in cartesian coordinates")!=-1:
        tick=0;
        for i in range(self.natoms):
          for m in range(3):
            for j in range(self.natoms):
              for n in range(3):
                line=lines[t+tick+5].split();
                for s in range(4):
                  line[s]=int(line[s]);
                line[4]=float(line[4]);
                self.dynmatrix[line[1]-1][line[3]-1][line[0]-1][line[2]-1]=line[4];
                tick=tick+1; 
            tick=tick+1; # remove one blank
    dyn.close();
  def obtainborn(self,dynfile):
    dyn=open(dynfile,'r');
    lines=dyn.readlines();
    for t in range(len(lines)):
      if lines[t].find("from electric field response")!=-1 and lines[t-1].find("Effective charges, in cartesian coordinates")!=-1:
        tick=0;
        for m in range(3):
          for i in range(self.natoms):
            for n in range(3):
              line=lines[t+tick+5].split();
              for s in range(4):
                line[s]=int(line[s]);
              line[4]=float(line[4]);
              self.atomcharge[line[1]-1][line[0]-1][line[2]-1]=line[4];
              tick=tick+1;
          tick=tick+1;# remove one blank
    dyn.close();
  def obtainedie(self,dynfile):
    dyn=open(dynfile,'r');
    lines=dyn.readlines();
    for t in range(len(lines)):
      if lines[t].find("Dielectric tensor, in cartesian coordinates")!=-1:
        tick=0;
        for m in range(3):
          for n in range(3+1):
            line=lines[t+tick+4].split();
            tick=tick+1;
            if len(line) < 4:
              continue;
            for s in range(4):
              line[s]=int(line[s]);
            line[4]=float(line[4]);
            self.edie[line[0]-1][line[2]-1]=line[4];
    dyn.close();
  def obtainforce(self,scfoutfile):
    scfout=open(scfoutfile,'r');
    lines=scfout.readlines();
    for t in range(len(lines)):
      if lines[t].find("cartesian forces (hartree/bohr) at end")!=-1:
        for i in range(self.natoms):
          line=lines[t+i+1].split();
          for j in range(3):
            self.force[i][j]=float(line[j+1]);
    scfout.close();
  def obtainpolarization(self,scfoutfile):
    scfout=open(scfoutfile,'r');
    lines=scfout.readlines();
    polar=np.zeros(3);
    for t in range(len(lines)):
      if lines[t].find("Polarization in cartesian coordinates (a.u.)")!=-1:
        line=lines[t+5].split();
        for j in range(-3,0,1):
          polar[j]=float(line[j]);
    scfout.close();
    return polar; # the units of Polarization is e/bohr^2;
  def obtainefield(self,scfin):
    f=open(scfin,'r');
    lines=f.readlines();
    efield=np.zeros(3);
    for t in range(len(lines)):
      if lines[t].find("efield")!=-1:
        line=lines[t].split();
        for j in range(-3,0,1):
          efield[j]=float(line[j]);
    f.close();
    return efield;
  def readposition(self,scfin,natoms):
    scfinput=open(scfin,'r');
    lines=scfinput.readlines();
    length=len(lines);
    atomp=np.zeros([natoms,3]);
    atomreturn=np.zeros([natoms,3]);
    for t in range(length):
      if lines[t].find("xred")!=-1:
        for i in range(natoms):
          line=lines[t+i].split();
          for j in range(-3,0,1):
            atomp[i][j]=float(line[j]);
        for i in range(natoms):
          for j in range(3):
            atomreturn[i]=atomreturn[i]+self.axis[j]*atomp[i][j];
      if lines[t].find("xcart")!=-1:
        for i in range(natoms):
          line=lines[t+i].split();
          if len(lines[t+natoms-1].split())==3:
            autoA=0.52917721067121;
            for i in range(natoms):
              line=lines[t+i].split();
              for j in range(-3,0,1):
                atomreturn[i][j]=float(line[j])*autoA;
          elif len(lines[t+natoms-1].split())==4:
            if lines[t+natoms-1].lower().find("bohr"):
              autoA=0.52917721067121;
            elif lines[t+natoms-1].lower().find("angstrom"):
              autoA=1.0;
            line=lines[t+0].split();
            for j in range(-3,0,1):
              atomreturn[0][j]=float(line[j])*autoA;
            for i in range(1,natoms):
              line=lines[t+i].split();
              for j in range(3):
                atomreturn[i][j]=float(line[j])*autoA;
    return atomreturn;
  def obtaindiff(self,bscfout,ascfout):
    pb=self.obtainpolarization(bscfout);
    pa=self.obtainpolarization(ascfout);
    return pa-pb
  def obtaindiffp_at_zeroE(self,bscfin,bscfzeroout,ascfin,ascfzeroout):
    print('-----------------------------------------');
    print('-- Polarization Difference by Position --')
    atombefore=self.readposition(bscfin,self.natoms);
    atomafter=self.readposition(ascfin,self.natoms);
    period=np.zeros(3);
    autoA=0.52917721067121;
    for i in range(3):
      period[i]=self.axis[i][i]/autoA;
    diffp=atomafter-atombefore;
    for i in range(self.natoms):
      for j in range(3):
        diffp[i][j]=diffp[i][j]-round(diffp[i][j]/period[j])*period[j];
    dipoleion=np.zeros(3);
    for i in range(self.natoms):
      dipoleion=dipoleion+np.matmul(self.atomcharge[i,0:3,0:3],diffp[i,0:3])/autoA;
    dp=self.obtaindiff(bscfzeroout,ascfzeroout);
    dp=dp*np.linalg.det(self.axis/autoA);
    print("Before Correction: ",dp,'File: ',bscfzeroout,ascfzeroout);
    for i in range(3):
      dp[i]=dp[i]-round(dp[i]/period[i])*period[i];
      temp=period[i]*round(dipoleion[i]/period[i]);
      dp[i]=dp[i]+temp;
      dp[i]=dp[i]-dipoleion[i]-round((dp[i]-dipoleion[i])/period[i])*period[i]+dipoleion[i];
    print("After Correction: ",dp,'Reference: ',dipoleion);
    print('-----------------------------------------');
    return dp
  def obtaindiffp_at_zeroP(self,scfin,scfout,scfzeroout):
    print('-----------------------------------------');
    print('-- Polarization Difference by Efield ----');
    dpe=self.obtaindiff(scfzeroout,scfout);
    autoA=0.529;
    dpe=dpe*np.linalg.det(self.axis/autoA);
    efield=self.obtainefield(scfin);
    epsilzero=8.8541878*10**-12;
    angstrom=10**-10;
    efieldabi=514220624373.482;
    echarge=1.60217662*10**-19;
    au=0.529*10**-10;
    dipoleqe=epsilzero*angstrom**3*efieldabi/echarge/au;
    polarestE=np.matmul(self.edie,efield)*np.linalg.det(self.axis)*dipoleqe;# e*au;
    period=np.zeros(3);
    for i in range(3):
      period[i]=self.axis[i][i]/autoA;
    print('----- Before Correction -----',dpe);
    for i in range(3):
      dpe[i]=dpe[i]-round(dpe[i]/period[i])*period[i];
      temp=period[i]*round(polarestE[i]/period[i]);
      dpe[i]=dpe[i]+temp;
      dpe[i]=dpe[i]-polarestE[i]-round((dpe[i]-polarestE[i])/period[i])*period[i]+polarestE[i];
    print('After Correction: ',dpe,'Reference: ',polarestE,'Efield: ',efield);
    print('-----------------------------------------');
    return dpe;
  def obtaindipolediffperiodtwo(self,bscfin,bscfout,bscfzeroout,ascfin,ascfout,ascfzeroout):
    autoA=0.529;
    total=self.obtaindiffp_at_zeroE(bscfin,bscfzeroout,ascfin,ascfzeroout);
    totalE2=self.obtaindiffp_at_zeroP(ascfin,ascfout,ascfzeroout);
    totalE1=self.obtaindiffp_at_zeroP(bscfin,bscfout,bscfzeroout);
    return (totalE2+total-totalE1)/np.linalg.det(self.axis/autoA);
  def writenewscf(self,Efield,atomposition,filename):
    efieldabi=514220624373.482;# Efield is Mv/cm, atmposition is $\AA$
    changeunits=Efield*10**6/10**(-2)/efieldabi;
    scffiles=open(self.scfin,'r');
    newfilename=open(filename,'w');
    lines=scffiles.readlines();
    skiplines=0;
    for i in range(len(lines)):
      if lines[i].find("efield")!=-1:
        newfilename.write("efield "+str(changeunits[0])+" "+str(changeunits[1])+" "+str(changeunits[2])+"\n");
        skiplines=skiplines+1;
      elif lines[i].find("xred")!=-1:
        newfilename.write("xred ");
        for j in range(self.natoms):
          temp="";
          atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(self.axis.transpose()));
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      elif lines[i].find("xcart")!=-1:
        if lines[i+self.natoms-1].lower().find('bohr')!=-1:
          autoA=0.529;
        elif lines[i+self.natoms-1].lower().find('angstrom')!=-1:
          autoA=1.0;
        newfilename.write("xcart ");
        for j in range(self.natoms):
          temp="";
          atomp=atomposition[j]/autoA;
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      if skiplines >0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
    newfilename.close();
    scffiles.close();
  def writenewscfnoe(self,atomposition,filename):
    scffiles=open(self.scfin,'r');
    newfilename=open(filename,'w');
    lines=scffiles.readlines();
    skiplines=0;
    for i in range(len(lines)):
      if lines[i].find('efield')!=-1:
        skiplines=skiplines+1;
      elif lines[i].find('berryopt')!=-1:
        skiplines=skiplines+1;
      elif lines[i].find("xred")!=-1:
        newfilename.write("xred ");
        for j in range(self.natoms):
          temp="";
          atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(self.axis.transpose()));
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      elif lines[i].find("xcart")!=-1:
        if lines[i+self.natoms-1].lower().find('bohr')!=-1:
          autoA=0.529;
        elif lines[i+self.natoms-1].lower().find('angstrom')!=-1:
          autoA=1.0;
        newfilename.write("xcart ");
        for j in range(self.natoms):
          temp="";
          atomp=atomposition[j]/autoA;
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      if skiplines >0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
    newfilename.close();
    scffiles.close();
  def writenewscfnoedipole(self,atomposition,filename):
    scffiles=open(self.scfin,'r');
    newfilename=open(filename,'w');
    lines=scffiles.readlines();
    skiplines=0;
    for i in range(len(lines)):
      if lines[i].find('efield')!=-1:
        newfilename.write('efield 0.0 0.0 0.0\n');
        skiplines=skiplines+1;
      elif lines[i].find("xred")!=-1:
        newfilename.write("xred ");
        for j in range(self.natoms):
          temp="";
          atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(self.axis.transpose()));
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      elif lines[i].find("xcart")!=-1:
        if lines[i+self.natoms-1].lower().find('bohr')!=-1:
          autoA=0.529;
        elif lines[i+self.natoms-1].lower().find('angstrom')!=-1:
          autoA=1.0;
        newfilename.write("xcart ");
        for j in range(self.natoms):
          temp="";
          atomp=atomposition[j]/autoA;
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      if skiplines >0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
    newfilename.close();
    scffiles.close();
  def writenewdfpt(self,atomposition,filename):
    phfiles=open(self.phin,'r');
    lines=phfiles.readlines();
    newfilename=open(filename,'w');
    skiplines=0;
    for i in range(len(lines)):
      if lines[i].find("xred")!=-1:
        newfilename.write("xred ");
        for j in range(self.natoms):
          temp="";
          atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(self.axis.transpose()));
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      elif lines[i].find("xcart")!=-1:
        if lines[i+self.natoms-1].lower().find('bohr')!=-1:
          autoA=0.529;
        elif lines[i+self.natoms-1].lower().find('angstrom')!=-1:
          autoA=1.0;
        newfilename.write("xcart ");
        for j in range(self.natoms):
          temp="";
          atomp=atomposition[j]/autoA;
          for k in range(3):
            temp=temp+' {:12.8f}'.format(atomp[k]);
          if j==0:
            newfilename.write(temp+'\n');
          else:
            newfilename.write("    "+temp+'\n');
          skiplines=skiplines+1;
      if skiplines >0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
  def solve(self,force,deltaP):
    #solve the linear equation Ax=y;
    self.A=np.zeros([3*self.natoms+3,3*self.natoms+3]);
    self.x=np.zeros(3*self.natoms+3);
    self.y=np.zeros(3*self.natoms+3);
    for i in range(self.natoms):
        for j in range(self.natoms):
            self.A[3*i:3*(i+1),3*j:3*(j+1)]=np.copy(self.dynmatrix[i,j,0:3,0:3]);
    for i in range(self.natoms):
        self.A[3*i:3*(i+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]*0.0003884098/2.0);
    for i in range(self.natoms):
        self.A[3*(self.natoms):3*(self.natoms+1),3*i:3*(i+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]);
    self.A[3*(self.natoms):3*(self.natoms+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.edie[0:3,0:3]*0.00001545*self.vol);
    for i in range(self.natoms):
        for j in range(3):
            self.y[i*3+j]=force[i][j];
    self.y[3*self.natoms+0]=self.vol*deltaP[0];
    self.y[3*self.natoms+1]=self.vol*deltaP[1];
    self.y[3*self.natoms+2]=self.vol*deltaP[2];
    self.x=np.matmul(np.linalg.inv(self.A),self.y);
    return np.matmul(np.linalg.inv(self.A),self.y);
  def iterateintermediate(self,startingpoint,times,step):
    startP=self.obtainpolarization(self.zerofile); #units e/bohr^2;
    endP=startP+np.array([0,0,1])*step/57.137;#units e/bohr^2;step units is C/m^2
    deltaP=startP-endP;
    print('the aim is: ',endP*57.137,'Please also be notifying that forces should also be zero');
    self.obtainedie(self.phfile);
    self.obtaindyn(self.phfile);
    # starting the first guess: th  ink the polarization increase partly come from the electrical contribution
    efieldabi=514220624373.482;
    efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.00001545; # the units of efield is Mv/cm
    startingefield=self.obtainefield('scf'+str(startingpoint))*efieldabi*(10**-8);
    efieldzero=efield+startingefield;#only count as 10% from the electrical contribution, the first step.
    au=0.52917721067121;
    self.writenewscf(efieldzero,self.atomp,"scf"+str(startingpoint+1));
    self.writenewscfnoedipole(self.atomp,'itenoe'+str(startingpoint+1));
    self.writenewdfpt(self.atomp,'dfpt'+str(startingpoint+1));
    accup=np.zeros([self.natoms,3])
    accue=np.zeros(3);
    accupolar=np.zeros(3);
    for i in range(startingpoint+1,startingpoint+times,1):
        self.obtainforce('scf'+str(i)+'.abo');
        scfbefore='scf'+str(i-1);
        scfafter='scf'+str(i);
        scfbeforeout='scf'+str(i-1)+'.abo';
        scfafterout='scf'+str(i)+'.abo';
        scfzerobeforeout="itenoe"+str(i-1)+'.abo';
        scfzeroafterout="itenoe"+str(i)+'.abo';
        diffp=self.obtaindipolediffperiodtwo(scfbefore,scfbeforeout,scfzerobeforeout,scfafter,scfafterout,scfzeroafterout);
        self.obtainedie('dfpt'+str(i)+'.abo');
        self.obtaindyn('dfpt'+str(i)+'.abo');
        accupolar=accupolar+diffp;
        print("The polarization difference is:(C/m^2) ",accupolar*57.137)
        deltax=self.solve(self.force,accupolar+startP-endP);
        accup=accup+au*np.reshape(deltax[0:self.natoms*3],[self.natoms,3]);
        accue=accue+deltax[self.natoms*3:(self.natoms+1)*3];
        posit=accup+self.atomp;
        efield=accue+efieldzero;
        self.writenewscf(efield,posit,'ite'+str(i+1));
        self.writenewscfnoedipole(posit,'itenoe'+str(i+1));
        self.writenewdfpt(posit,'dfpt'+str(i+1));
abi=abiout("./scf0","./scf0.abo","./dfpt0","./dfpt0.abo");
startingpoint=0;
endpoint=int(sys.argv[1]);
abi.iterateintermediate(startingpoint,endpoint,float(sys.argv[2]));
