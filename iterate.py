import numpy as np
class iteration:
    def __init__(self, natom, goal, Pth, Fth):
        self.goal = goal; # goal of difference polarization,
        self.Pth = Pth; # Polarization difference should be converged;
        self.Fth = Fth; # force should be equilibriumed after iteration;
        self.natom = natom;
        self.Efieldchain = []; # E field has the unit Mv/cm.
        self.Positchain = []; # Position has the unit angstrom.
        self.Polarization = []; # Polarization has the unit C/m^2;

    def loading(self, Efield, Position, axis, Polar):
        self.Echain.append(Efield); # E field has the unit Mv/cm
        self.Pchain.append(Position); # Position has the unit angstrom
        self.Polarization.append(Polar); # Polarization sequence. units C/m^2
        self.axis = axis; # Axis for volume has the units Angstrom.

    def nextsolution(self, force, dynmatrix, atomcharge, edie, Polarnow):
        self.Polarizaion.append(Polarnow);
        DP = Polarnow - self.goal;
        period = self.axis *( 10**(-10) )* 1.60217663 * ( 10**(-19) ) / ( 10**(-30) ) # Units C / m^2.
        DP = DP - np.round(DP / period) * period;
        Amatrix = np.zeros((3 * self.natom + 3, 3 * self.natom + 3), dtype = float);
        x = np.zeros(3 * self.natom + 3);
        y = np.zeros(3 * self.natom + 3);
        for i in range(self.natom):
            for j in range(self.natom):
                Amatrix[(3 * i) : (3 * i + 3), (3 * j) : (3 * j + 3)] = np.copy(dynmatrix[i][j]);
        for i in range(self.natom):
            Amatrix[(3 * i) : (3 * i + 3), (3 * self.natom) : (3 * self.natom + 3)] = np.copy( -1 * atomcharge[i] * 0.01); # Together makes the unit to eV / A when times the E field with Mv/cm.
        for i in range(self.natom):
            Amatrix[(3 * self.natom) : (3 *self.natom + 3), (3 * i) : (3 * i + 3)] = np.copy( -1 * atomcharge[i]);
        columbtoelectron = 1 / (1.60217663 * (10**-19));
        unitconversion = 10**6 / (10**-2) * (8.85418 * 10**-12) * columbtoelectron / ( 10**10 )**2;
        Amatrix[(3 * self.natom) : (3 * self.natom + 3), (3 * self.natom) : (3 * self.natom + 3)] = np.copy(-1 * edie * unitconversion * np.linalg.det(self.axis));
        for i in range(self.natom):
            for j in range(3):
                y[i * 3 + j] = force[i][j];
        y[(3 * self.natom) : (3 * self.natom + 3)] = np.linalg.det(self.axis) * DP * columbtoelectron / ( 10**10 )**2;
        x = np.matmul(np.linalg.inv(Amatrix), y);
        dPosition = x[0 : (3 * self.natom)];
        dEfield = x[(3 * self.natom) : (3 * self.natom + 3)];
        self.Echain.append(self.Echain + dEfield);
        self.Pchain.append(self.Pchain + dPosition);
