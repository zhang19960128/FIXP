import numpy as np
import math
import re
class Penergy:
    def __init__(self):
        pass;
    def obtainefield(self,scfin):
        result=np.zeros(3);
        try:
            f=open(scfin,'r');
            lines=f.readlines();
            f.close();
        except IOError:
            return result;
        ebefore=np.zeros(3);
        for i in range(len(lines)):
            for j in range(3):
                if lines[i].find("efield_cart("+str(j+1)+")")!=-1:
                    ebefore[j]=float(lines[i].split("=")[-1]);
        return ebefore;
    def obtainenergy(self,scfout):
        result=-1;
        try:
            f=open(scfout,'r')
            lines=f.readlines();
            f.close();
        except IOError:
            return result;
        number=[];
        for i in range(len(lines)):
            if lines[i].find("!")!=-1:
                number=re.findall("[-]{0,1}[0-9]{0,10}[\.]{1}[0-9]{0,10}",lines[i]);
        if len(number)!=0:
            result=float(number[0]);
        return result
    def obtaindipolescf(self,file):
        berryout=open(file,'r');
        lines=berryout.readlines();
        berryout.close();
        epolar=np.zeros(3);
        apolar=np.zeros(3);
        total=np.zeros(3);
        for i in range(len(lines)):
            if(lines[i].find("Electronic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    epolar[j]=float(lines[i+1+j].split()[1]);
            elif(lines[i].find("Ionic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    apolar[j]=float(lines[i+1+j].split()[1]);
        for i in range(3):
            total[i]=epolar[i]+apolar[i];
        return total;
    def obtainvolume(self,scfin):
        try:
            f=open(scfin,'r');
            lines=f.readlines();
            f.close();
        except IOError:
            return -1;
        volume=np.zeros([3,3]);
        for i in range(len(lines)):
            if lines[i].find("CELL_PARAMETERS")!=-1:
                for j in range(3):
                    for k in range(3):
                        volume[j,k]=float(lines[i+j+1].split()[k]);
        return np.linalg.det(volume);
    def scan(self,path,ran):
        tick=[];
        j=len(tick);
        for i in range(1,ran+1):
            tick.append(self.obtainenergy(path+"ite.out"+str(i)));
        for i in range(len(tick)):
            if np.abs(tick[i]-(-1))<1e-8:
                j=i;
                break;
        energy=self.obtainenergy(path+"ite.out"+str(j));
        efield=self.obtainefield(path+"ite"+str(j));
        dp=self.obtaindipolescf(path+'ite.out'+str(j));
        return [efield,energy,dp,j];
    def converg(self,path,ran):
        tick=[];
        for i in range(1,ran+1):
            tick.append(self.obtainenergy(path+"ite.out"+str(i)));
        for i in range(len(tick)):
            if np.abs(tick[i]-(-1))<1e-8:
                j=i;
                break;
        energy=self.obtainenergy(path+"ite.out"+str(j));
        efieldunit=36.3609*10**10;
        volume=8.1435518265*8.2519273758*8.2056713104*10**-30;
        ctoe=6.24150975*10**18;
        etoc=1.602*10**-19;
        au=0.529*10**-10;
        rytoj=13.59*1.60218*10**-19;
        for i in range(1,j+1):
            energy=self.obtainenergy(path+"ite.out"+str(i));
            efield=self.obtainefield(path+'ite'+str(i));
            dp=self.obtaindipolescf(path+'ite.out'+str(j));
            plug=np.inner(dp,efield)*1.0/math.sqrt(2)*au*etoc*efieldunit/rytoj;
            print('raw print=',energy,'legend=',plug,'final=(mev)',(energy+plug)*1000*13.59);
    def scanP(self,path,Plist):
        length=len(Plist);
        HPlist=[];
        efieldunit=36.3609*10**10;
        volume=8.1435518265*8.2519273758*8.2056713104*10**-30;
        ctoe=6.24150975*10**18;
        etoc=1.602*10**-19;
        au=0.529*10**-10;
        rytoj=13.59*1.60218*10**-19;
        for i in range(0,length):
            result=self.scan(path+"P"+"%3.2f"%Plist[i]+"/",20);
            plug=np.inner(result[2],result[0])*1.0/math.sqrt(2)*au*etoc*efieldunit/rytoj;
            H0=(plug+result[1])*13.59*1000;
            HPlist.append(H0);
        return np.array(HPlist)-np.min(np.array(HPlist))
    def scanE(self,path,Plist):
        length=len(Plist);
        efieldunit=36.3609*10**10;
        elist=[];
        for i in range(0,length):
            result=self.scan(path+'P'+"%3.2f"%Plist[i]+'/',20);
            elist.append(result[0]*efieldunit/10**8);
        return elist;
if __name__!='__main__':
  exit();
else:
  from FIXP import pwout
Pe=Penergy();
Pstart=np.array([4.63807888e-05,-4.80890700e-07,3.8454e-01])
#pscan=np.arange(-0.30+0.000000001,0.15,0.05);
pscan=np.arange(-0.06,0.12,0.02);
pscan=pscan.tolist();
#pscan=[-0.20,-0.15]+pscan+[0.15,0.20]
result=Pe.scanP("/p/work1/jiahaoz/Polar_iteration_05/",pscan);
for i in range(len(result)):
    print("%3.2f"%pscan[i],result[i])
resultE=Pe.scanE("/p/work1/jiahaoz/Polar_iteration_05/",pscan);
for i in range(len(resultE)):
    print("%3.2f"%pscan[i],resultE[i])
print('===== converge test =======')
Pe.converg("/p/work1/jiahaoz/Polar_iteration_05/P-0.10/",20);
startingpoint=0;
pw=pwout("./P0.00/PWOUT","ph.out"+str(startingpoint),'dyn.out'+str(startingpoint),'ite.out'+str(startingpoint),'ite'+str(startingpoint));
pw.obtain(20);
dpe=[];
dpi=[];
for i in range(len(pscan)):
    tick=Pe.scan("/p/work1/jiahaoz/Polar_iteration_05/P"+"%3.2f"%pscan[i]+"/PWOUT/",20)[3];
    path="/p/work1/jiahaoz/Polar_iteration_05/P"+"%3.2f"%pscan[i]+"/PWOUT";
    scfbefore=path+'/'+'ite'+str(0);
    scfafter=path+'/'+'ite'+str(tick);
    scfbeforeout=path+'/'+'ite.out'+str(0);
    scfafterout=path+'/'+'ite.out'+str(tick);
    scfzerobeforeout=path+'/'+"itezerofield.out"+str(0);
    scfzeroafterout=path+'/'+"itezerofield.out"+str(tick);
    diffp=pw.obtaindipolediffperiodtwo(scfbefore,scfbeforeout,scfzerobeforeout,scfafter,scfafterout,scfzeroafterout);
    dpe.append(diffp[1]);
    dpi.append(diffp[2]);
for i in range(len(pscan)):
    print("%8.6f"%pscan[i],"%8.6f"%(dpe[i]*57.137)[0],"%8.6f"%(dpe[i]*57.137)[1],"%8.6f"%(dpe[i]*57.137)[2],"%8.6f"%(dpi[i]*57.137)[0],"%8.6f"%(dpi[i]*57.137)[1],"%8.6f"%(dpi[i]*57.137)[2])
