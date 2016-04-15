import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
from topol import Topol_Data
import numpy as np

class RMSF_measurements(object):
    """Measures RMSF of ligand atoms over a single trajectory."""
    def __init__(self,topol_object, topology, trajectory, ligand_name, offset, output_name):
        self.ligand_rmsf = {}
        self.universe = topol_object
        self.measure_ligand_rmsf(topology, trajectory, ligand_name,offset,output_name)
        self.universe.rmsf_analysis = True
        self.min_value = min(self.ligand_rmsf.values())
        self.max_value = max(self.ligand_rmsf.values())
    def measure_ligand_rmsf(self,topology, trajectory, ligand_name,offset,output_name):
        i=0
        rmsf_to_average = {}
        for traj in trajectory:
            rmsf_to_average[i]=[]
            mobile = MDAnalysis.Universe(topology,traj)
            reference = MDAnalysis.Universe(topology)
            MDAnalysis.analysis.align.rms_fit_trj(mobile, reference, filename='test.xtc',select='segid '+str(ligand_name.segids[0])+' and resid '+str(ligand_name.resids[0]))
            md_sim = Topol_Data(topology,"test.xtc",ligand_name,offset)
            ligand_no_H=md_sim.universe.select_atoms('segid '+str(ligand_name.segids[0])+' and resid '+str(ligand_name.resids[0])+" and not name H*")
            R = MDAnalysis.analysis.rms.RMSF(ligand_no_H)
            R.run()
            rmsf_list = R.rmsf.tolist()
            rmsf_to_average[i]=rmsf_list
            i+=1
        for k in range(ligand_no_H.n_atoms):
            self.ligand_rmsf[k] = np.array([atom[k] for atom in rmsf_to_average.values()]).mean()
            
        with open(str(output_name)+"_rmsf_results.dat","w+") as rmsf_file:
            rmsf_file.write("Atom Number  Atom Name  RMSF \n")
            for i in range(ligand_no_H.n_atoms):
                results = str(ligand_no_H.atoms[i].index)+"   "+str(ligand_no_H.atoms[i].name)+"   "+str(rmsf_list[i]) +"\n"
                rmsf_file.write(results)