import numpy as np
class iteration:
    def __init__(self, natom, goal, Pth, Fth):
        self.goal = goal; # goal of total polarization,
        self.Pth = Pth; # Polarization difference should be converged;
        self.Fth = Fth; # force should be equilibriumed after iteration;
        self.natom = natom;
        self.Echain = []; # E field has the unit Mv/cm.
        self.Pchain = []; # Position has the unit angstrom.
        self.Polarization = []; # Polarization has the unit C/m^2;

    def loading(self, Efield, Position, axis, Polar):
        self.Echain.append(Efield); # E field has the unit Mv/cm
        self.Pchain.append(Position); # Position has the unit angstrom
        self.Polarization.append(Polar); # Polarization sequence. units C/m^2
        self.axis = axis; # Axis for volume has the units Angstrom.

    def nextsolution(self, force, dynmatrix, atomcharge, edie, Polarnow, tick):
        self.Polarization.append(Polarnow);
        columbtoelectron = 1 / (1.60217663 * (10**-19));
        epsilon0 = 8.85418 * ( 10**-12 );
        DP = Polarnow - self.goal;
        period = self.axis / np.linalg.det(self.axis) * ( 10**(-10) ) * 1.60217663 * ( 10**(-19) ) / ( 10**(-30) ) # Units C / m^2.
        period = np.array([period[0][0], period[1][1], period[2][2]]);
        DP = DP - np.round(DP / period) * period;
        DP = DP * columbtoelectron / ( 10**10 )**2; # Units of DP is now e / A^2, we use the Efield as V/A
        Amatrix = np.zeros((3 * self.natom + 3, 3 * self.natom + 3), dtype = float);
        x = np.zeros(3 * self.natom + 3);
        y = np.zeros(3 * self.natom + 3);
        for i in range(self.natom):
            for j in range(self.natom):
                Amatrix[(3 * i) : (3 * i + 3), (3 * j) : (3 * j + 3)] = np.copy(dynmatrix[i][j]);
        for i in range(self.natom):
            Amatrix[(3 * i) : (3 * i + 3), (3 * self.natom) : (3 * self.natom + 3)] = np.copy( -1 * atomcharge[i]);
        for i in range(self.natom):
            Amatrix[(3 * self.natom) : (3 *self.natom + 3), (3 * i) : (3 * i + 3)] = np.copy( -1 * atomcharge[i]);
        epsilon0 = epsilon0 * columbtoelectron / 10**( 10 ); # now the units are e V-1 / A
        unitconversion = epsilon0;
        Amatrix[(3 * self.natom) : (3 * self.natom + 3), (3 * self.natom) : (3 * self.natom + 3)] = np.copy(-1 * edie * unitconversion * np.linalg.det(self.axis));
        np.savetxt("Amatrix.dat{0:d}".format(tick), Amatrix, fmt="%10.5f", delimiter=' ')
        for i in range(self.natom):
            for j in range(3):
                y[i * 3 + j] = force[i][j];
        y[(3 * self.natom) : (3 * self.natom + 3)] = np.linalg.det(self.axis) * DP;
        y = y.reshape((self.natom * 3 + 3, 1));
        np.savetxt("y.dat{0:d}".format(tick), y, fmt="%10.5f", delimiter=' ');
        x = np.matmul(np.linalg.inv(Amatrix), y);
        x = x.reshape(3 * self.natom + 3)
        dPosition = x[0 : (3 * self.natom)];
        dEfield = x[(3 * self.natom) : (3 * self.natom + 3)];
        dEfield = dEfield * 100; # 100 is to convert V/ang -> Mv/cm
        self.Echain.append(self.Echain[-1] + dEfield);
        newP = self.Pchain[-1] + dPosition.reshape((self.natom, 3));
        for i in range(self.natom - 1, -1, -1):
            newP[i] = newP[i] - newP[0];
        self.Pchain.append(newP);
