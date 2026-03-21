# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 19:27:24 2025

@author: CGIBEIL
"""

import torch
import clip
import os, glob
from torchvision.datasets import CIFAR100
import numpy as np

clip.available_models()

device = "cuda" if torch.cuda.is_available() else "cpu"
model, preprocess = clip.load('ViT-L/14', device)

os.chdir("C:\\Users\\cgibeil@emory.edu\\OneDrive - Emory\\Documents\\PhD\\Projects\\P1-CondProbs\\Core_scripts\\Cellular_targets_pipeline\\Marina\\CLIP")

cifar100 = CIFAR100(root=os.path.expanduser("~/.cache"), download=True, train=False)


image, class_id = cifar100[3637]
image_input = preprocess(image).unsqueeze(0).to(device)

from PIL import Image

## imdir
imdir = "C:\\Users\\cgibeil@emory.edu\\OneDrive - Emory\\Documents\\PhD\\Projects\\P1-CondProbs\\Core_scripts\\Cellular_targets_pipeline\\Marina\\CLIP\\patches\\"
os.chdir(imdir)

images = glob.glob("*.png")

embeddings = []
with torch.no_grad():
    for img_path in images:
        image = Image.open(img_path).convert("RGB")
        image_input = preprocess(image).unsqueeze(0).to(device)
        image_features = model.encode_image(image_input)
        embeddings.append(image_features.cpu())
        print(img_path)

face_probe = clip.tokenize(["image of a face"]).to(device)
body_probe = clip.tokenize(["image of a body"]).to(device)
background_probe = clip.tokenize(["image of background"]).to(device)
foreground_probe = clip.tokenize(["image of foreground object"]).to(device)


embeddings = torch.cat(embeddings, dim=0)  # shape: [num_images, embedding_dim]

embeddings_np = embeddings.numpy()
embeddings_norm = embeddings_np / np.linalg.norm(embeddings_np, axis=1, keepdims=True)

os.chdir("C:\\Users\\cgibeil@emory.edu\\OneDrive - Emory\\Documents\\PhD\\Projects\\P1-CondProbs\\Core_scripts\\Cellular_targets_pipeline\\Marina\\CLIP")
np.savetxt("embeddings_all_ViTL.csv",embeddings_norm)

     
images = [i.replace('-','_').replace('.png', '') for i in images]
clip_target = np.array([i.split("_") for i in images])
     
from umap import UMAP
reducer = UMAP(
       n_neighbors=10,  
       min_dist=0.1,        
       metric='correlation',
       random_state=402,
       n_components = 3
   )
embeddings_3d = reducer.fit_transform(embeddings_norm)

import pandas as pd
embeddings_labelled = pd.DataFrame(np.append(clip_target.astype(float)[:,0:2], embeddings_3d, axis=1))
embeddings_labelled.columns = ["clip", "target", "x", "y", "z"]
embeddings_labelled.to_csv("all_embeddings2.csv")

## Text probes
probes = [face_probe, body_probe, background_probe, foreground_probe]

txt_embeddings = []
for i,probe in enumerate(probes):
    text_features = model.encode_text(probe)
    txt_embeddings.append(text_features.cpu())
    
txt_embeddings = torch.cat(txt_embeddings, dim=0)  # shape: [num_images, embedding_dim]
txt_embeddings_np = txt_embeddings.detach().numpy()
txt_embeddings_norm = txt_embeddings_np / np.linalg.norm(txt_embeddings_np, axis=1, keepdims=True)





embeddings_1 = pd.DataFrame(np.append(embeddings_labelled[embeddings_labelled["clip"]==335][["target"]], embeddings_3d, axis=1))
embeddings_1.columns = ["target", "x", "y", "z"]
embeddings_1.to_csv("embeddings_335.csv")
