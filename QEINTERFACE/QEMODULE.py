import numpy as np
import sys
class ABIIO:
    def __init__(self, natom, scfEtemplate, scfnoEtemplate):
        self.scfEtemplate = scfEtemplate;
        self.scfnoEtemplate = scfnoEtemplate;
        self.natom = natom;
        self.axis = np.zeros((3, 3));
        self.atomp = np.zeros((natom, 3, 3));
        self.force = np.zeros((natom, 3, 3));

    def readscf(self):
        scfinput = open(self.scfin,'r');
        lines = scfinput.readlines();
        length = len(lines);
        for i in range(length):
            if lines[i].find("CELL_PARAMETERS") != -1:
                if lines[i].find("(angstrom)") != -1:
                    for j in range(3):
                        for k in range(3):  
                            self.axis[j][k] = float(lines[i + j + 1].split()[k]);
                else:
                    sys.exit("!!ERROR!! I ONLY SUPPORT ANGSTROM");
            if lines[i].find("ATOMIC_POSITIONS") != -1:
                if lines[i].find("(angstrom)") != -1:
                    for j in range(self.natoms):
                        for k in range(3):
                            self.atomp[j,k] = float(lines[i + j + 1].split()[k + 1]);
                else:
                    sys.exit('!!ERROR!! I ONLY SUPPORT ANGSTROM')

    def obtainforce(self, filename): 
        scfout = open(filename, 'r'); # post-processing the scf.out
        lines = scfout.readlines();
        for i in range(len(lines)):
            if lines[i].find("Forces acting on atoms (cartesian axes, Ry/au):") != -1:
                for j in range(self.natom):
                    self.force[j][0] = float(lines[i + j + 2].split()[6]);
                    self.force[j][1] = float(lines[i + j + 2].split()[7]);
                    self.force[j][2] = float(lines[i + j + 2].split()[8]);
        scfout.close();

    def obtaindipolescf(self, filename):
        berryout = open(filename, 'r');
        lines = berryout.readlines();
        epolar = np.zeros(3);
        apolar = np.zeros(3);
        total = np.zeros(3);
        Atoau = 0.52917721067121;
        for i in range(len(lines)):
            if(lines[i].find("Electronic Dipole on Cartesian axes") != -1):
                for j in range(3):
                    epolar[j] = float(lines[i + 1 + j].split()[1]) / math.sqrt(2);
            elif(lines[i].find("Ionic Dipole on Cartesian axes") != -1):
                for j in range(3):
                    apolar[j] = float(lines[i + 1 + j].split()[1]) / math.sqrt(2);
        for i in range(3):
            total[i] = epolar[i] + apolar[i];
        # now the units of dipole is e * bohr   
        print('the polarization now is: ',total/np.linalg.det(self.axis/Atoau)*57.137)
        return total/np.linalg.det(self.axis/Atoau);

    def writescfE(self, Efield, atomposition, outputfilename):
        changeunits = Efield*10**6/10**(-2)/(36.3509*10**10); # the Efield units now is MV/cm the atomposition units is the angstrom;
        scffiles = open(self.scfEtemplate, 'r');
        newfilename = open(outputfilename, 'w');
        lines = scffiles.readlines();
        tick1 = 0;
        tick2 = 0;
        for i in range(len(lines)):
            if lines[i].find("efield_cart") != -1:
                tick1 = tick1 + 1;
                if tick1 > 1:
                    continue;
                else:
                    for j in range(3):
                        newfilename.write("efield_cart"+"("+str(j+1)+")="+'{:12.7f}'.format(changeunits[j])+"\n");
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)") != -1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2 = i;
                for j in range(self.natoms):
                    temp = "";
                    for t in range(3):
                        temp = temp+" " + '{:12.8f}'.format(atomposition[j, t]);
                    newfilename.write(lines[j+i+1].split()[0]+" "+temp+"\n");
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            elif lines[i].find("verbosity") != -1:
                continue;
            elif lines[i].find("wf_collect") != -1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();

    def writenewscfnoe(self, atomposition, outputfilename):
        scffiles = open(self.scfnoEtemplate,'r');
        newfilename = open(outputfilename,'w');
        lines = scffiles.readlines();
        tick1 = 0;
        tick2 = 0;
        skiplines = 0;
        for i in range(len(lines)):
            if skiplines > 0:
                skiplines = skiplines - 1;
                continue;
            if lines[i].find("efield_cart") != -1:
                continue;
            elif lines[i].find("lelfield") != -1:
                continue;
            elif lines[i].find("ATOMIC_POSITIONS (angstrom)") != -1:
                newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
                tick2 = i;
                for j in range(self.natoms):
                    temp = "";
                    for t in range(3):
                        temp = temp+'{:12.8f}'.format(atomposition[j,t]);
                    newfilename.write(lines[j + i + 1].split()[0] + " " + temp + '\n');
            elif lines[i].find("K_POINTS") != -1:
                newfilename.write("K_POINTS automatic\n");
                newfilename.write("4 4 4 1 1 1\n")
                skiplines = 1;
            elif tick2+self.natoms >= i and tick2 >1:
                continue;
            else:
                newfilename.write(lines[i]);
        newfilename.close();
        scffiles.close();

class PHIO:
    def __init__(self, natom):
        self.natom = natom;
        self.edie = np.zeros((3, 3));
        self.atomcharge = np.zeros((natom, 3, 3));

    def obtainph(self, filename):
        phout = open(filename, 'r');
        lines = phout.readlines();
        for i in range(len(lines)):
            if lines[i].find("site n.  atom      mass           positions (alat units)") != -1:
                for j in range(self.natoms):
                    for k in range(3):
                        self.atommass[j] = float(lines[i+j+1].split()[2]);
            if lines[i].find("Dielectric constant in cartesian axis") != -1:
                for j in range(3):
                    for k in range(3):
                        self.edie[j][k] = float(lines[i + j + 2].split()[k + 1]);
            if lines[i].find("Effective charges (d P / du) in cartesian axis") != -1:
                for j in range(self.natoms):
                    for k in range(3):
                        for m in range(3):
                            self.atomcharge[j][k][m] = lines[i + 2 + 4 * j + k + 1].split()[m + 2];
        self.edie = self.edie-np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
        phout.close();

    def obtaindyn(self, filename):
        dyn = open(filename, 'r');
        lines = dyn.readlines();
        for i in range(len(lines)):
            if lines[i].find("Dynamical  Matrix in cartesian axes") != -1:
                for j in range(self.natoms):
                    for k in range(self.natoms):
                        for m in range(3):
                            self.dynmatrix[j][k][0][m] = float(lines[i + 4 + 4 * (j * self.natoms + k) + 1].split()[m * 2]);
                            self.dynmatrix[j][k][1][m] = float(lines[i + 4 + 4 * (j * self.natoms + k) + 2].split()[m * 2]);                    
                            self.dynmatrix[j][k][2][m] = float(lines[i + 4 + 4 * (j * self.natoms + k) + 3].split()[m * 2]);
        for i in range(self.natoms):
            for j in range(self.natoms):
                self.dynmatrix[i][j] = self.dynmatrix[i][j]; #be careful about w and f, those are circle frequence and normal frequency.
        dyn.close();
