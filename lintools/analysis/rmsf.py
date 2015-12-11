import MDAnalysis
from MDAnalysis.analysis import rms

class RMSF_measurements(object):
    """Measures RMSF of ligand atoms over a single trajectory."""
    def __init__(self,topol_object):
        self.ligand_rmsf = {}
        self.universe = topol_object
        self.measure_ligand_rmsf()
        self.universe.rmsf_analysis = True
        self.min_value = min(self.ligand_rmsf.values())
        self.max_value = max(self.ligand_rmsf.values())
    def measure_ligand_rmsf(self):
        R = MDAnalysis.analysis.rms.RMSF(self.universe.ligand_no_H)
        R.run()
        rmsf_list = R.rmsf.tolist()
        for i in range(self.universe.ligand_no_H.n_atoms):
            res_tuple = (i,rmsf_list[i])
            self.ligand_rmsf[i] = rmsf_list[i]