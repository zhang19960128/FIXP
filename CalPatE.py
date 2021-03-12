from FIXP import pwout
import os
from update import updatedft
from update import obtainnatoms
class scanE(pwout):
  def updatescflist(self,pathlist,scfinlist,scfoutlist):
      length=len(scfinlist);
      for i in range(length):
        updatedft(pathlist[i],scfinlist[i],scfoutlist[i],'angstrom')
  def scanP(self,scfinlist,scfoutlist,scfzerooutlist):
      length=len(scfinlist);
      natoms=obtainnatoms(scfinlist[0]);
      print(natoms)
      self.obtain(natoms);
      dpdiff=[];
      for i in range(length):
        dp=obtaindipolediffperiodtwo(self.scfin,self.zerofile,self.zerofile,scfinlist[i],scfoutlist[i],scfzerooutlist[i]);
        dpdiff.append(dp);
      return dpdiff;
  def absolutep(self,highsymmscfin,highsymmscfout):
      dp=obtaindipolediffperiodtwo(highsymmscfin,highsymmscfout,highsymmscfout,self.scfin,self.zerofile,self.zerofile);
      return dp;
  def obtainenergy(self,scfoutlist):
      energylist=[];
      for i in scfoutlist:
        syscmd='grep "End final" '+i+' | wc -l';
        re=os.popen(syscmd).read();
        re=int(re);
        if re==1:
          syscmd='grep ! '+i+' | tail -1 | grep -Eo "[-]{1}[0-9]{0,9}[\.]{1}[0-9]{0,9}"'
          energy=os.popen(syscmd).read();
          energylist.append(float(energy));
        else:
          print(i,'CANNOT FIND CONVERGED RESULT');
          sys.exit();
      return energylist;
pwEscf=scanE('./PWOUT','ph.out0','dyn.out0','ite.out0','ite0')
path="/workspace/jiahaoz/BiFeO3/EfieldEnergy/05_phase/";
scfoutlist=[];
scfinlist=[];
pathlist=[];
for i in range(-7,8,1):
  ptemp=path+"E{0:02d}".format(i)+'/bto.out';
  scfoutlist.append(ptemp);
  ptemp=path+"E{0:02d}".format(i)+'/bto.in';
  scfinlist.append(ptemp);
  ptemp=path+"E{0:02d}".format(i)+'/';
  pathlist.append(ptemp);
elist=pwEscf.obtainenergy(scfoutlist)
pwEscf.updatescflist(pathlist,scfinlist,scfoutlist)
#pwEscf.scanP(scfinlist,scfoutlist,)
