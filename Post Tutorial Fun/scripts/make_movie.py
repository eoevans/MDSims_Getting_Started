import cv2
import numpy as np
import os
import sys
from PIL import Image

frames_folder = 'MDSim/data/simulation_results/C10H22/CHARMM/run1/frames/'
im_pillow = np.array(Image.open(os.path.join(frames_folder, 'frame.0000.rgb')))
im_bgr = cv2.cvtColor(im_pillow, cv2.COLOR_RGB2BGR)
height, width, layers = im_bgr.shape
fps = 24 

video = cv2.VideoWriter('MDSim/data/simulation_results/C10H22/CHARMM/run1/cool_movie.avi',cv2.VideoWriter_fourcc(*'DIVX'), fps, (width,height))
 
for image in os.scandir(frames_folder): 
    im_pillow = np.array(PIL.Image.open(os.path.join(frames_folder, image)))
    im_bgr = cv2.cvtColor(im_pillow, cv2.COLOR_RGB2BGR)
    video.write(im_bgr)
      
# Deallocating memories taken for window creation
cv2.destroyAllWindows() 
video.release()  # releasing the video generated