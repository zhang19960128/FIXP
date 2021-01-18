import numpy as np
import re
import math
import sys
class pwout:
    def __init__(self,filepath,phout,dynout,zero,scfin):
        self.path=filepath;
        self.phfile=self.path+"/"+phout;
        self.dynfile=self.path+'/'+dynout;
        self.zerofile=self.path+'/'+zero;
        self.scfin=self.path+'/'+scfin;
    def obtain(self,natom):
        au=0.52917721067121;
        self.natoms=natom;
        self.axis=np.zeros([3,3]);
        self.atomp=np.zeros([self.natoms,3]);
        self.atommass=np.zeros([self.natoms,1]);
        self.atomcharge=np.zeros([self.natoms,3,3]);
        self.atomdis=np.arange(self.natoms*3);
        self.dynmatrix=np.zeros([self.natoms,self.natoms,3,3]);
        self.edie=np.zeros([3,3]);
        self.force=np.zeros([self.natoms,3]);
        self.polar=np.zeros(3);
        self.readscf();
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        self.obtainforce(self.zerofile);
        self.vol=np.linalg.det(pw.axis)/au/au/au;
    def obtainph(self,file):
        phout=open(file,'r');
        #post-processing ph.out
        lines=phout.readlines();
        for i in range(len(lines)):
            if lines[i].find("site n.  atom      mass           positions (alat units)")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        self.atommass[j]=float(lines[i+j+1].split()[2]);
            if lines[i].find("Dielectric constant in cartesian axis")!=-1:
                for j in range(3):
                    for k in range(3):
                        self.edie[j][k]=float(lines[i+j+2].split()[k+1]);
            if lines[i].find("Effective charges (d P / du) in cartesian axis")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        for m in range(3):
                            self.atomcharge[j][k][m]=lines[i+2+4*j+k+1].split()[m+2];
        phout.close();
    def obtaindyn(self,file):
        dyn=open(file);
        # post-processing BFO.dyn
        lines=dyn.readlines();
        for i in range(len(lines)):
            if lines[i].find("Dynamical  Matrix in cartesian axes")!=-1:
                for j in range(self.natoms):
                    for k in range(self.natoms):
                        for m in range(3):
                            self.dynmatrix[j][k][0][m]=float(lines[i+4+4*(j*self.natoms+k)+1].split()[m*2]);
                            self.dynmatrix[j][k][1][m]=float(lines[i+4+4*(j*self.natoms+k)+2].split()[m*2]);                    
                            self.dynmatrix[j][k][2][m]=float(lines[i+4+4*(j*self.natoms+k)+3].split()[m*2]);
        for i in range(self.natoms):
            for j in range(self.natoms):
                #be careful about w and f, those are circle frequence and normal frequency.
                self.dynmatrix[i][j]=self.dynmatrix[i][j];
        dyn.close();
    def obtainforce(self,file): 
        # post-processing the scf.out
        scfout=open(file,'r');
        lines=scfout.readlines();
        for i in range(len(lines)):
            if lines[i].find("Forces acting on atoms (cartesian axes, Ry/au):")!=-1:
                for j in range(self.natoms):
                    self.force[j][0]=float(lines[i+j+2].split()[6]);
                    self.force[j][1]=float(lines[i+j+2].split()[7]);
                    self.force[j][2]=float(lines[i+j+2].split()[8]);
        scfout.close();
    def obtaindipolescf(self,file):
        berryout=open(file,'r');
        lines=berryout.readlines();
        epolar=np.zeros(3);
        apolar=np.zeros(3);
        total=np.zeros(3);
        Atoau=0.52917721067121;
        for i in range(len(lines)):
            if(lines[i].find("Electronic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    epolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
            elif(lines[i].find("Ionic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    apolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
        for i in range(3):
            total[i]=epolar[i]+apolar[i];
        # now the units of dipole is e * bohr   
        print('the polarization now is: ',total/np.linalg.det(self.axis/Atoau)*57.137)
        return total/np.linalg.det(self.axis/Atoau);
    def obtaindipolescftotal(self,file):
        berryout=open(file,'r');
        lines=berryout.readlines();
        epolar=np.zeros(3);
        apolar=np.zeros(3);
        total=np.zeros(3);
        Atoau=0.52917721067121;
        for i in range(len(lines)):
            if(lines[i].find("Electronic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    epolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
            elif(lines[i].find("Ionic Dipole on Cartesian axes")!=-1):
                for j in range(3):
                    apolar[j]=float(lines[i+1+j].split()[1])/math.sqrt(2);
        return[epolar,apolar];
    def setpspcharge(self):
        atomcharge=np.zeros([self.natoms,3,3]);
        for i in range(int(self.natoms*1/5)):
            for j in range(3):
                atomcharge[i,j,j]=16;# charge of Fe
        for i in range(int(self.natoms*1/5),int(self.natoms*2/5)):
            for j in range(3):
                atomcharge[i,j,j]=15;# charge of Bi
        for i in range(int(self.natoms*2/5),int(self.natoms)):
            for j in range(3):
                atomcharge[i,j,j]=6;
        return atomcharge;
    def obtainefield(self,scfin):
        f=open(scfin,'r');
        lines=f.readlines();
        f.close();
        ebefore=np.zeros(3);
        for i in range(len(lines)):
            for j in range(3):
                if lines[i].find("efield_cart("+str(j+1)+")")!=-1:
                    ebefore[j]=float(lines[i].split("=")[-1]);
        return ebefore;
    def obtaindiff(self,bscfout,ascfout):
        d1=self.obtaindipolescftotal(ascfout);
        d2=self.obtaindipolescftotal(bscfout);
        de=d2[0]-d1[0];
        di=d2[1]-d1[1];
        return [de,di];
    def obtaindipolediffperiodtwo(self,bscfin,bscfout,bscfzeroout,ascfin,ascfout,ascfzeroout):
#       print(bscfzeroout,'-----Obtain zero file Dipole change:');
        dp=self.obtaindiff(ascfzeroout,bscfzeroout);
#       print(dp);
#       print(ascfzeroout,'---------------------------------------------------------------');
        epsilzero=8.8541878*10**-12;
        angstrom=10**-10;
        efieldqe=36.3609*10**10;
        echarge=1.60217662*10**-19;
        au=0.529*10**-10;
        Atoau=0.529;
        dipoleqe=epsilzero*angstrom**3*efieldqe/echarge/au;
        atombefore=self.readposition(bscfin,self.natoms);
        atomafter=self.readposition(ascfin,self.natoms);
        period=np.zeros(3);
        for i in range(3):
            period[i]=self.axis[i][i]/0.529;
        diffp=atomafter-atombefore;
        dipoleion=np.zeros(3);
        for i in range(self.natoms):
            dipoleion=dipoleion+np.matmul(self.atomcharge[i,0:3,0:3],diffp[i,0:3])/Atoau;
        total=dp[0]+dp[1];
        print('-Correct Finite position difference change of Polariztion,Delta(#2,#1)---')
        print("Before Correction: ",total,'File:',ascfzeroout,bscfzeroout);
        for i in range(3):
            total[i]=total[i]-round(total[i]/period[i])*period[i];
            temp=period[i]*round(dipoleion[i]/period[i]);
            total[i]=total[i]+temp;
            total[i]=total[i]-dipoleion[i]-round((total[i]-dipoleion[i])/period[i])*period[i]+dipoleion[i]
        print("After Correction: ",total,'Reference: ',dipoleion);
        print('---------------------------------------------------------------------------');
        print('-Correct Finite EField Dipole change by Polarization Quanta, Delta(#2+E,#2+0)-',ascfout,ascfzeroout);     
        dpE2=self.obtaindiff(ascfout,ascfzeroout);
        efield=self.obtainefield(ascfin);
        polarestE2=np.matmul(self.edie,efield)*np.linalg.det(self.axis)*dipoleqe;# units e*au;
        totalE2=dpE2[0]+dpE2[1];
        print('Before Correction: ',totalE2);
        for i in range(3):
            totalE2[i]=totalE2[i]-round(totalE2[i]/period[i])*period[i];
            temp=period[i]*round(polarestE2[i]/period[i]);
            totalE2[i]=totalE2[i]+temp;
            totalE2[i]=totalE2[i]-polarestE2[i]-round((totalE2[i]-polarestE2[i])/period[i])*period[i]+polarestE2[i];
        print('After Correction: ',totalE2,"Reference:",polarestE2,'Efield',efield);
        print('---------------------------------------------------------------------------');
        print('-Correct Finite EField Dipole change by Polarization Quanta,Delta(#1+E,#1+0)-',bscfout,bscfzeroout);
        dpE1=self.obtaindiff(bscfout,bscfzeroout);
        efield=self.obtainefield(bscfin);
        # estimate the electrical polarization change;
        polarestE1=np.matmul(self.edie,efield)*np.linalg.det(self.axis)*dipoleqe;# units e*au;
        totalE1=dpE1[0]+dpE1[1];
        print('Before Correction: ',totalE1);
        for i in range(3):
            totalE1[i]=totalE1[i]-round(totalE1[i]/period[i])*period[i];
            temp=period[i]*round(polarestE1[i]/period[i]);
            totalE1[i]=totalE1[i]+temp;
            totalE1[i]=totalE1[i]-polarestE1[i]-round((totalE1[i]-polarestE1[i])/period[i])*period[i]+polarestE1[i];
        print('After Correction: ',totalE1,"Reference",polarestE1,'Efield',efield);
        print('----------------------------------------------------------------------------');
        return (totalE2+total-totalE1)/np.linalg.det(self.axis/Atoau);
    def readscf(self):
        scfinput=open(self.scfin,'r');
        lines=scfinput.readlines();
        length=len(lines);
        for i in range(length):
            if lines[i].find("CELL_PARAMETERS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(3):
                    for k in range(3):  
                        self.axis[j][k]=float(lines[i+j+1].split()[k]);
            if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(self.natoms):
                    for k in range(3):
                        self.atomp[j,k]=float(lines[i+j+1].split()[k+1]);
    def readposition(self,scfin,natoms):
        scfinput=open(scfin,'r');
        lines=scfinput.readlines();
        length=len(lines);
        atomp=np.zeros([natoms,3]);
        for i in range(length):
            if lines[i].find("ATOMIC_POSITIONS")!=-1 and lines[i].find("(angstrom)")!=-1:
                for j in range(natoms):
                    for k in range(3):
                        atomp[j,k]=float(lines[i+j+1].split()[k+1]);
        return atomp;
    def writenewscf(self,Efield,atomposition,filename):
        # the Efield units now is MV/cm
        # the atomposition units is the angstrom;
        changeunits=Efield*10**6/10**(-2)/(36.3509*10**10);
        scffiles=open(self.scfin,'r');
        newfilename=open(self.path+'/'+filename,'w');
        lines=scffiles.readlines();
        tick1=0;
        tick2=0;
        for i in range(len(lines)):
            if lines[i].find("efield_cart")!=-1:
                tick1=tick1+1;
                if tick1>1:
                    continue;
                else:
                    for j in range(3):
                        newfilename.write("efield_cart"+"("+str(j+1)+")="+':12.7f'.format(changeunits[j])+"\n");
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2=i;
                for j in range(self.natoms):
                    temp="";
                    for t in range(3):
                        temp=temp+" "+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+"\n");
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            elif lines[i].find("verbosity")!=-1:
                continue;
            elif lines[i].find("wf_collect")!=-1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();
    def writenewscfnoe(self,atomposition,filename):
        scffiles=open(self.scfin,'r');
        newfilename=open(self.path+'/'+filename,'w');
        lines=scffiles.readlines();
        tick1=0;
        tick2=0;
        skiplines=0;
        for i in range(len(lines)):
            if skiplines >0:
                skiplines=skiplines-1;
                continue;
            if lines[i].find("efield_cart")!=-1:
                continue;
            elif lines[i].find("lelfield")!=-1:
                continue;
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2=i;
                for j in range(self.natoms):
                    temp="";
                    for t in range(3):
                        temp=temp+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+'\n');
            elif lines[i].find("K_POINTS")!=-1:
                newfilename.write("K_POINTS automatic\n");
                newfilename.write("4 4 4 1 1 1\n")
                skiplines=1;
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();  
    def writenewscfnoedipole(self,atomposition,filename):
        scffiles=open(self.scfin,'r');
        newfilename=open(self.path+'/'+filename,'w');
        lines=scffiles.readlines();
        tick1=0;
        tick2=0;
        skiplines=0;
        for i in range(len(lines)):
            if skiplines >0:
                skiplines=skiplines-1;
                continue;
            if lines[i].find("efield_cart")!=-1:
                newfilename.write("efield_cart(1)=0.0\n");
                newfilename.write("efield_cart(2)=0.0\n");
                newfilename.write("efield_cart(3)=0.0\n");
                skiplines=skiplines+3;
                continue;
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2=i;
                for j in range(self.natoms):
                    temp="";
                    for t in range(3):
                        temp=temp+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+'\n');
            elif lines[i].find("K_POINTS")!=-1:
                newfilename.write("K_POINTS automatic\n");
                newfilename.write("6 6 6 1 1 1\n")
                skiplines=1;
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close(); 
    def solve(self,force,deltaP):
        #solve the linear equation Ax=y;
        self.A=np.zeros([3*self.natoms+3,3*self.natoms+3]);
        self.x=np.zeros(3*self.natoms+3);
        self.y=np.zeros(3*self.natoms+3);
        for i in range(self.natoms):
            for j in range(self.natoms):
                self.A[3*i:3*(i+1),3*j:3*(j+1)]=np.copy(self.dynmatrix[i,j,0:3,0:3]);
        for i in range(self.natoms):
            self.A[3*i:3*(i+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]*0.0003884098);
        for i in range(self.natoms):
            self.A[3*(self.natoms):3*(self.natoms+1),3*i:3*(i+1)]=np.copy(-1*self.atomcharge[i,0:3,0:3]);
        self.A[3*(self.natoms):3*(self.natoms+1),3*self.natoms:3*(self.natoms+1)]=np.copy(-1*self.edie[0:3,0:3]*0.0000154961*self.vol);
        for i in range(self.natoms):
            for j in range(3):
                self.y[i*3+j]=force[i][j];
        self.y[3*self.natoms+0]=self.vol*deltaP[0];
        self.y[3*self.natoms+1]=self.vol*deltaP[1];
        self.y[3*self.natoms+2]=self.vol*deltaP[2];
        self.x=np.matmul(np.linalg.inv(self.A),self.y);
        return np.matmul(np.linalg.inv(self.A),self.y);
    def iteratetwo(self,times,step):
        startP=self.obtaindipolescf(self.zerofile); #units e/bohr^2;
        endP=startP+np.array([0,0,1])*step/57.137;#units e/bohr^2;step units is C/m^2
        deltaP=startP-endP;
        print('the aim is: ',endP*57.137,'Please also be notifying that forces should also be zero');
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        # starting the first guess: th  ink the polarization increase partly come from the electrical contribution
        efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.0000154961;
        efieldzero=efield*1.0;#only count as 10% from the electrical contribution, the first step.
        au=0.52917721067121;
        self.writenewscf(efield,self.atomp,"ite1");
        self.writenewscfnoe(self.atomp,"itenoe1");
        self.writenewscfnoedipole(self.atomp,'itenoedipole'+str(1));
        efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.0000154961;
        efieldzero=efield*1.0;#only count as 10% from the electrical contribution, the first step.
        au=0.52917721067121;
        self.writenewscf(efield,self.atomp,"ite1");
        self.writenewscfnoe(self.atomp,"itenoe1");
        self.writenewscfnoedipole(self.atomp,'itenoedipole'+str(1));
        accup=np.zeros([self.natoms,3])
        accue=np.zeros(3);
        accupolar=np.zeros(3);
        for i in range(1,times):
            self.obtainforce(self.path+'/'+'ite.out'+str(i));
            scfbefore=self.path+'/'+'ite'+str(i-1);
            scfafter=self.path+'/'+'ite'+str(i);
            scfbeforeout=self.path+'/'+'ite.out'+str(i-1);
            scfafterout=self.path+'/'+'ite.out'+str(i);
            scfzerobeforeout=self.path+'/'+"itezerofield.out"+str(i-1);
            scfzeroafterout=self.path+'/'+"itezerofield.out"+str(i);
            diffp=self.obtaindipolediffperiodtwo(scfbefore,scfbeforeout,scfzerobeforeout,scfafter,scfafterout,scfzeroafterout);
            self.obtainph(self.path+'/'+'ph.out'+str(i));
            self.obtaindyn(self.path+'/'+'dyn.out'+str(i));
#            deltaP=diffp+startP-endP;
            accupolar=accupolar+diffp;
            print("The polarization difference is:(C/m^2) ",accupolar*57.137)
            deltax=self.solve(self.force,accupolar+startP-endP);
            accup=accup+au*np.reshape(deltax[0:self.natoms*3],[self.natoms,3]);
            accue=accue+deltax[self.natoms*3:(self.natoms+1)*3];
            posit=accup+self.atomp;
            efield=accue+efieldzero;
            self.writenewscf(efield,posit,'ite'+str(i+1));
            self.writenewscfnoe(posit,'itenoe'+str(i+1));
            self.writenewscfnoedipole(posit,'itenoedipole'+str(i+1));
        print('====================================================')
    def iterateintermediate(self,startingpoint,times,step):
        startP=self.obtaindipolescf(self.zerofile); #units e/bohr^2;
        endP=startP+np.array([0,0,1])*step/57.137;#units e/bohr^2;step units is C/m^2
        deltaP=startP-endP;
        print('the aim is: ',endP*57.137,'Please also be notifying that forces should also be zero');
        self.obtainph(self.phfile);
        self.obtaindyn(self.dynfile);
        # starting the first guess: th  ink the polarization increase partly come from the electrical contribution
        efield=np.matmul(-1*np.linalg.inv(self.edie*self.vol),deltaP*self.vol)/0.0000154961;
        startingefield=self.obtainefield(self.path+'/'+'ite'+str(startingpoint));
        efieldzero=efield+startingefield;#only count as 10% from the electrical contribution, the first step.
        au=0.52917721067121;
        self.writenewscf(efieldzero,self.atomp,"ite"+str(startingpoint+1));
        self.writenewscfnoe(self.atomp,"itenoe"+str(startingpoint+1));
        self.writenewscfnoedipole(self.atomp,'itenoedipole'+str(startingpoint+1));
        accup=np.zeros([self.natoms,3])
        accue=np.zeros(3);
        accupolar=np.zeros(3);
        for i in range(startingpoint+1,times,1):
            self.obtainforce(self.path+'/'+'ite.out'+str(i));
            scfbefore=self.path+'/'+'ite'+str(i-1);
            scfafter=self.path+'/'+'ite'+str(i);
            scfbeforeout=self.path+'/'+'ite.out'+str(i-1);
            scfafterout=self.path+'/'+'ite.out'+str(i);
            scfzerobeforeout=self.path+'/'+"itezerofield.out"+str(i-1);
            scfzeroafterout=self.path+'/'+"itezerofield.out"+str(i);
            diffp=self.obtaindipolediffperiodtwo(scfbefore,scfbeforeout,scfzerobeforeout,scfafter,scfafterout,scfzeroafterout);
            self.obtainph(self.path+'/'+'ph.out'+str(i));
            self.obtaindyn(self.path+'/'+'dyn.out'+str(i));
#            deltaP=diffp+startP-endP;
            accupolar=accupolar+diffp;
            print("The polarization difference is:(C/m^2) ",accupolar*57.137)
            deltax=self.solve(self.force,accupolar+startP-endP);
            accup=accup+au*np.reshape(deltax[0:self.natoms*3],[self.natoms,3]);
            accue=accue+deltax[self.natoms*3:(self.natoms+1)*3];
            posit=accup+self.atomp;
            efield=accue+efieldzero;
            self.writenewscf(efield,posit,'ite'+str(i+1));
            self.writenewscfnoe(posit,'itenoe'+str(i+1));
            self.writenewscfnoedipole(posit,'itenoedipole'+str(i+1));
startingpoint=0;
pw=pwout("./PWOUT","ph.out"+str(startingpoint),'bto.dyn'+str(startingpoint),'ite.out'+str(startingpoint),'ite'+str(startingpoint));
pw.obtain(5);
pw.iterateintermediate(int(sys.argv[1]),startingpoint,float(sys.argv[2]));
