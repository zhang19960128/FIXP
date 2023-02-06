import numpy as np
import sys
class ABIIO:
    def __init__(self, natom, scfEtemplate, scfnoEtemplate):
        self.scfEtemplate = scfEtemplate;
        self.scfnoEtemplate = scfnoEtemplate;
        self.natom = natom;
        self.axis = np.zeros((3, 3));
        self.atomp = np.zeros((natom, 3));
        self.force = np.zeros((natom, 3));

    def readscf(self):
        scfinput = open(self.scfnoEtemplate,'r');
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
                    for j in range( self.natom ):
                        for k in range(3):
                            self.atomp[j, k] = float(lines[i + j + 1].split()[k + 1]);
                else:
                    sys.exit('!!ERROR!! I ONLY SUPPORT ANGSTROM')
        scfinput.close();

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
        Ry = 13.605;
        au = 0.5291;
        self.force = self.force * Ry / au;

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
                    epolar[j] = float(lines[i + 1 + j].split()[1]) / np.sqrt(2);
            elif(lines[i].find("Ionic Dipole on Cartesian axes") != -1):
                for j in range(3):
                    apolar[j] = float(lines[i + 1 + j].split()[1]) / np.sqrt(2);
        for i in range(3):
            total[i] = epolar[i] + apolar[i];
        # now the units of dipole is e * bohr
        echarge = 1.60217663 * 10**(-19);
        bohrR = 5.291772 * 10**(-11);
        polarization = total / np.linalg.det(self.axis/Atoau) * echarge * bohrR / (bohrR)**3;
        print('the polarization now is: ', polarization)
        return polarization;

    def writescfE(self, Efield, atomposition, outputfilename):
        changeunits = Efield * 10**6 / 10**(-2) / (36.3509*10**10); # the Efield units now is MV/cm the atomposition units is the angstrom;
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
        self.atommass = np.zeros(natom);
        self.dynmatrix = np.zeros((natom, natom, 3, 3), dtype = float);

    def obtainph(self, filename):
        phout = open(filename, 'r');
        lines = phout.readlines();
        for i in range(len(lines)):
            if lines[i].find("site n.  atom      mass           positions (alat units)") != -1:
                for j in range(self.natom):
                    for k in range(3):
                        self.atommass[j] = float(lines[i+j+1].split()[2]);
            if lines[i].find("Dielectric constant in cartesian axis") != -1:
                for j in range(3):
                    for k in range(3):
                        self.edie[j][k] = float(lines[i + j + 2].split()[k + 1]);
            if lines[i].find("Effective charges (d Force / dE) in cartesian axis with asr applied:") != -1:
                for j in range(self.natom):
                    for k in range(3):
                        for m in range(3):
                            self.atomcharge[j][k][m] = lines[i + 2 + 4 * j + k].split()[m + 2];
        self.edie = self.edie-np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]);
        phout.close();

    def obtaindyn(self, filename):
        dyn = open(filename, 'r');
        lines = dyn.readlines();
        for i in range(len(lines)):
            if lines[i].find("Dynamical  Matrix in cartesian axes") != -1:
                for j in range(self.natom):
                    for k in range(self.natom):
                        for m in range(3):
                            self.dynmatrix[j][k][0][m] = float(lines[i + 4 + 4 * (j * self.natom + k) + 1].split()[m * 2]);
                            self.dynmatrix[j][k][1][m] = float(lines[i + 4 + 4 * (j * self.natom + k) + 2].split()[m * 2]);
                            self.dynmatrix[j][k][2][m] = float(lines[i + 4 + 4 * (j * self.natom + k) + 3].split()[m * 2]);
        dyn.close(); # the raw output is not mass weighted the following lines will help determining the mass-weighted dynamical matrix and help.
        Ry = 13.605;
        au = 0.5291;
        self.dynmatrix = self.dynmatrix * Ry / au / au; # Ry = 13.605, au = 0.5291
        self.dynmatrix_massweighted = np.copy(self.dynmatrix);
        for i in range(self.natom):
            for j in range(self.natom):
                self.dynmatrix_massweighted[i][j] = self.dynmatrix[i][j] / np.sqrt(self.atommass[i] * self.atommass[j]);
        self.dynmatrix1D = np.zeros((3 * self.natom, 3 * self.natom));
        self.dynmatrix_massweighted1D = np.zeros((3 * self.natom, 3 * self.natom));
        for i in range(self.natom):
            for j in range(self.natom):
                for m in range(3):
                    for n in range(3):
                        self.dynmatrix1D[3 * i + m][3 * j + n] = self.dynmatrix[i][j][m][n];#be careful about w and f, those are circle frequence and normal frequency. (http://ilan.schnell-web.net/physics/rydberg.pdf) Units of mass weighted Dynmatrix should be  Ry/au/au.
                        self.dynmatrix_massweighted1D[3 * i + m, 3 * j + n] = self.dynmatrix_massweighted[i][j][m][n]; # units are Ry / au / au / proton mass

    def obtainborn(self, filename):
        f = open(filename, 'r');
        line = f.readlines();
        for i in range(len(line)):
            if lines[i].find("Effective charges (d Force / dE) in cartesian axis with asr applied") != -1:
                for j in range(self.natom):
                    for k in range(4):
                        if k != 0:
                            lin = line[i + 4 * j + k + 1].split();
                            for n in range(3):
                                self.atomcharge[j][k - 1][n] = lin[2 + n];
                        else:
                            continue;
        f.close();
