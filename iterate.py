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
        self.axis = axis; # Axis for volume has the units Angstrom

    def nextsolution(self, dynmatrix):
        Amatrix = np.zeros((3 * self.natom + 3, 3 * self.natom + 3), dtype = float);
        x = np.zeros(3 * self.natom + 3);
        y = np.zeros(3 * self.natom + 3);
        for i in range(self.natom):
            for j in range(self.natom):
                Amatrix[(3 * i) : (3 * i + 3), (3 * j) : (3 * j + 3)] = np.copy(dynmatrix[i][j]);
        for i in range(self.natom):
            