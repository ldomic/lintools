import sys


class ProgressBar(object):
    def __init__(self, title):
        self.title = title
    def update(self, frame_n, total_n_frames):
        percentage = frame_n/float(total_n_frames)*100
        percent_text = str(int(percentage))+"%"
        if frame_n==1 or frame_n%10==0:
            sys.stdout.write("   "+ self.title +"   frame    "+ str(frame_n)+"/"+str(total_n_frames)+"   ["+percent_text+"]"+ "\r")
            sys.stdout.flush()
    def finish(self, total_n_frames):
        sys.stdout.write("   "+ self.title +"   frame    "+ str(total_n_frames)+"/"+str(total_n_frames)+"   [100%]"+ "\r")
