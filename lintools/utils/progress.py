import sys


class ProgressBar(object):
    def __init__(self, title):
        self.title = title
    def update(self, frame_n, total_n_frames):
        percentage = frame_n/float(total_n_frames)*100
        percent_text = str(int(percentage))+"%"
        if frame_n==1 or frame_n%11==0:
            sys.stdout.write("   "+ self.title +"   frame    "+ str(frame_n)+"/"+str(total_n_frames)+"   ["+percent_text+"]"+ "\r")
            sys.stdout.flush()
        elif frame_n==total_n_frames+1:
            sys.stdout.write("   "+ self.title +"   frame    "+ str(frame_n)+"/"+str(total_n_frames)+"   ["+percent_text+"]"+ "\r")
