
class Analysis_Data(object):
    def __init__(self, topology_data_object, lintools_object):
        self.analysis_data = {}
        self.topology_data = topology_data_object
        self.lintools = lintools_object

    def get_input_data(self):
        self.analysis_data["input data"]={
                        "topology": self.topology_data.topology,
                        "trajectory": self.topology_data.trajectory,
                        # TODO: add more things
        }
    def get_per_atom_data_raw(self):
        self.analysis_data["raw"]={atom:{
                    # Residence
                    "contacts":# TODO: function()

        } for atom in self.lintools.res_time.contacts_per_atom.keys() }

    def get_atom_residence_data(self):
        num_records = np.sum([1 for frame in self.timeseries])
        dtype = [("frame",float),("time",float),("ligand atom id",int),
                ("ligand atom name","|U4"),("distance",float),
                ("resid",int),("resname","|U4"),("segid","|U8") ]
        out = np.empty((num_records,),dtype=dtype)
        cursor=0
        for contact in self.timeseries:
            out[cursor] = (contact.frame, contact.time,contact.ligandatomid,contact.ligandatomname,contact.distance,
                           contact.resid,contact.resname,contact.segid)
            cursor+=1
        return out.view(np.recarray)
