class G4Material():
    def __init__(self, label='', density=0.0, matid=0):
        self.label = label
        self.density = density
        self.matid = matid

    def __repr__(self):
        return '<Material: "{!s}", {:d}, {:f}>'.format(self.label, self.matid, self.density)

    def get_material_def(self):
        """returns entry suitable for voxel description in geometry.text file (format: [dens, matid] )"""
        return "{:0.3f} {:d}\n".format(self.density, self.matid)

mat_map = {
    "air":       G4Material('Air',               density=0.00129, matid=0),
    "lung_in":   G4Material('Lungs (inhale)',    density=0.217,   matid=1),
    "lung_ex":   G4Material('Lungs (exhale)',    density=0.508,   matid=2),
    "adipose":   G4Material('Adipose',           density=0.967,   matid=3),
    "breast":    G4Material('Breast',            density=0.99,    matid=4),
    "water":     G4Material('Water',             density=1.0,     matid=5),
    "muscle":    G4Material('Muscle',            density=1.061,   matid=6),
    "liver":     G4Material('Liver',             density=1.071,   matid=7),
    "bone_trab": G4Material('Bone (trabecular)', density=1.159,   matid=8),
    "bone_comp": G4Material('Bone (compact)',    density=1.575,   matid=9),
}
