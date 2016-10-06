import MDAnalysis
from collections import Counter
from timeit import default_timer as timer
from utils.progressbar import ProgressBar
import sys



class Residence_time(object):
    __version__ = "09.2016"
    __authors__ = "Laura Domicevica","Tom Newport"
    """
    This module analyses the residence time of each protein residue. If a residue in within
    the cutoff distance, it will be added to the list, but only those residues that spend a
    significant portion of simulation time in the vicinity of the ligand will be plotted at
    the end of the analysis.
    
    Takes:
        * topology data object * - information about system (lintools.Data object)
        * trajectory * - list of trajectories
        * start_frame_num * - list^^ of frame numbers for start of analysis (Opt)
        * end_frame_num * - list^^ of frame numbers for end of analysis (Opt) 
        * skip * - list^^ of numbers of how many frames should be skipped for this analysis (Opt)
    
    
    ^^ The reason these values are lists is because several trajectories can be submitted for 
    analysis and different values could be required for each simulation. Therefore, similarly 
    as with trajectories, start, end and skip variables are submitted as lists with values 
    corresponding for each trajectory.
    
    Example: trajectory = ["1.xtc","2.xtc"] #Both are 1000 frames, but the user only wants to 
             analyse second half the the second trajectory
             start = [0(for the first traj),500(for the second traj)]   
             Other values can be left as default.
    """
    def __init__(self, topology_data_object, trajectory,start_frame_num,end_frame_num,skip, topology,mol2,ligand,offset):
        self.residue_counts = {}
        self.frequency ={}
        self.frame_count = []
        self.topology_data = topology_data_object
        self.trajectory = trajectory
        self.start_frame_num = start_frame_num
        self.end_frame_num = end_frame_num
        self.skip = skip
        self.topology = topology
        self.mol2_file = mol2
        self.ligand = ligand
        self.offset = offset
    def measure_residence_time(self, cutoff = 3.5):
        """
        This function measures the residence time of all residues within the cutoff distance from
        the ligand.
        Takes:
            * cutoff * - cutoff distance in angstroms
        Output:
            * self.residue_counts * - dictionary of all residues that have made contact with the
            ligand (i.e. have been within the cutoff distance) and the number of frames of the
            contact
        """
        i=0
        for traj in self.trajectory:
            frame_dict = {}
            start = timer()
            prog = ProgressBar("Residence_time")
            self.topology_data.load_trajectory(traj)
            for frame in self.topology_data.universe.trajectory[self.start_frame_num[i]:self.end_frame_num[i]:self.skip[i]]:  
                selection = self.topology_data.universe.select_atoms("protein and around "+str(cutoff)+' (segid '+str(self.topology_data.universe.ligand.segids[0])+' and resid '+str(self.topology_data.universe.ligand.resids[0])+')')
                #selection by residues - they can be accessed later to make sure on which chain they are
                #also counter treats each residue as single object - i.e. no problems with multiple chains where
                #there are several residues with the same resid and resname - they will be treated seperately
                prog.update(frame.frame,self.topology_data.universe.trajectory.n_frames)
                residue_list = [residue for residue in selection.residues]
                frame_dict[frame.time]=set(residue_list)
                
            # How many frames each individual simulation has? The simulations can be of differing length    
            self.frame_count.append(len(frame_dict))
            prog.finish(self.topology_data.universe.trajectory.n_frames)
            self.residue_counts[i]=Counter([item for sublist in frame_dict.values() for item in sublist])
            end = timer()

            sys.stdout.write("Residence time: "+ str(end-start))
            #Assert whether there is anything to plot for this simulation
            assert len(self.residue_counts[i])!=0, "Nothing to draw for this ligand (residue number: "+ self.topology_data.universe.ligand.resids[0] +" on the chain "+ self.topology_data.universe.ligand.segids[0] +") - check the position of your ligand within the trajectory."
            i+=1
    def define_residues_for_plotting_traj(self, analysis_cutoff):
        """
        Since plotting all residues that have made contact with the ligand over a lenghty 
        simulation is not always feasible or desirable. Therefore, only the residues that 
        have been in contact with ligand for a long amount of time will be plotted in the
        final image. 
        
        The function first determines the fraction of time each residue spends in the 
        vicinity of the ligand for each trajectory. Once the data is processed, analysis
        cutoff decides whether or not these residues are plotted based on the total 
        frequency this residue has spent in the vicinity of the ligand. The analysis
        cutoff is supplied for a single trajectory and is therefore multiplied.
        
        Takes:
            * analysis_cutoff * - a fraction (of time) a residue has to spend in the 
            vicinity of the ligand for a single traj
        Output:
            * self.frequency * - frequency per residue per trajectory
            * topol_data.dict_of_plotted_res * - the residues that should be plotted in 
            the final image with the frequency for each trajectory (used for plotting)
        """
        self.residue_counts_fraction = {}
        #Calculate the fraction of time a residue spends in each simulation
        for traj in self.residue_counts:
            self.residue_counts_fraction[traj] = {residue:float(values)/self.frame_count[traj] for residue,values in self.residue_counts[traj].items()}
        
        for traj in self.residue_counts_fraction:
            for residue in self.residue_counts_fraction[traj]:
                try:
                    self.frequency[residue].append(self.residue_counts_fraction[traj][residue])
                except KeyError:
                    self.frequency[residue] = []
                    self.frequency[residue].append(self.residue_counts_fraction[traj][residue])

        self.topology_data.dict_of_plotted_res = {i:self.frequency[i] for i in self.frequency if sum(self.frequency[i])>(int(len(self.trajectory))*analysis_cutoff)}

        assert len(self.topology_data.dict_of_plotted_res)!=0,"Nothing to draw for this ligand:(residue number: "+ str(self.topology_data.universe.ligand.resids[0]) +" on the chain "+ str(self.topology_data.universe.ligand.segids[0]) +") - try reducing the analysis cutoff."

