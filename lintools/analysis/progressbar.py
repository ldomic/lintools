import sys

class ProgressBar(object):
    __author__ = "Laura Domicevica"
    __version__ = "09.2016"
    """
    Shows the progress of analysis in frame rate and percentage. 
    In August 2016, used for residence time analysis, since others are shown by MDAnalysis.
    """
    def __init__(self, title):
        self.title = title
    def update(self, frame_n, total_n_frames):
        percentage = frame_n/float(total_n_frames)*100
        percent_text = str(int(percentage))+"%"
        sys.stdout.write(self.title +"   frame    "+ str(frame_n)+"/"+str(total_n_frames)+"   ["+percent_text+"]"+ "\r")
        sys.stdout.flush()
    def finish(self, total_n_frames):
        sys.stdout.write(self.title +"   frame    "+ str(total_n_frames)+"/"+str(total_n_frames)+"   [100%]"+ "\r")