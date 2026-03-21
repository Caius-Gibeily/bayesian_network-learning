# -*- coding: utf-8 -*-

from scenedetect import detect, ContentDetector
import numpy as np
import os
workingdir = "bn_analysis/"
os.chdir(workingdir)

video_dir = r"clips/"
clips = [335, 340, 348, 351, 352]

for c,clip in enumerate(clips):    
    scene_list = detect(f"{video_dir}0{clip}PEER_movie.mov", 
                        ContentDetector(threshold = 15,
                                        min_scene_len=45))
    frames = [cut[0].frame_num+1 for cut in scene_list] 
    np.savetxt(f"clip{clip}.csv", frames)
    