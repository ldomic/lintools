import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.align import *
from collections import defaultdict
from timeit import default_timer as timer


class RMSF_measurements(object):
    """Measures RMSF of ligand atoms over a single trajectory."""
    def __init__(self,topology_data_object, topology, trajectory,ligand_name):
        self.ligand_rmsf = defaultdict(int)
        self.topology_data = topology_data_object
        self.trajectory = trajectory
        self.topology = topology
        self.measure_ligand_rmsf(ligand_name)
        self.min_value = min(self.ligand_rmsf.values())
        self.max_value = max(self.ligand_rmsf.values())
    def measure_ligand_rmsf(self,ligand_name):
        i=0
        rmsf_list={}
        start = timer()
        for traj in self.trajectory:
            self.topology_data.universe.load_new(traj)
            reference = MDAnalysis.Universe(self.topology)
            MDAnalysis.analysis.align.rms_fit_trj(self.topology_data.universe, reference, filename='test.xtc',select='protein')
            aligned_universe = MDAnalysis.Universe(self.topology,"test.xtc")
            ligand_noH = aligned_universe.select_atoms('(segid '+str(ligand_name.segids[0])+' and resid '+str(ligand_name.resids[0])+')')
            R = MDAnalysis.analysis.rms.RMSF(ligand_noH)
            R.run()
            rmsf_list[i] = R.rmsf.tolist()

            i+=1

        for index,atom in enumerate(self.topology_data.universe.ligand_noH.atoms):
            for traj in rmsf_list:
                self.ligand_rmsf[atom.name] += rmsf_list[traj][index]

        self.ligand_rmsf  = {k:v/len(self.trajectory) for k,v in self.ligand_rmsf.items()}
        print "RMSF: "+str(timer()-start)