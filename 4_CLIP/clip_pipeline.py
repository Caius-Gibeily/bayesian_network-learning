# -*- coding: utf-8 -*-
"""
@author: CGIBEIL
"""
import os
import glob
import torch
import clip
import numpy as np
import pandas as pd
from PIL import Image
from umap import UMAP
from sklearn.metrics.pairwise import cosine_similarity

def compute_patch_embeddings(model, preprocess, image_paths, device, batch_size=64):
    """

    Parameters
    ----------
    model : Which Image model to use - ViNet, 
    preprocess : Preprocessed model 
    image_paths : path to patches
    device : hardware device: cpu/CUDA-supported gpu
    batch_size : number of patches to load and process together

    Returns
    -------
    Numpified embeddings (N patches x M features)

    """
    all_embs = []
    with torch.no_grad():
        for i in range(0, len(image_paths), batch_size):
            batch_files = image_paths[i:i+batch_size]
            imgs = [preprocess(Image.open(p).convert("RGB")) for p in batch_files]
            batch = torch.stack(imgs).to(device, non_blocking=True)
            emb = model.encode_image(batch).cpu()
            emb /= emb.norm(dim=1, keepdim=True)
            all_embs.append(emb)
            print(f"Encoded batch {i//batch_size+1}/{int(np.ceil(len(image_paths)/batch_size))}")
            del batch, emb  # free GPU memory
            torch.cuda.empty_cache()
    return torch.cat(all_embs, dim=0).numpy()

def compute_text_embeddings(model, device, probe_list, batch_size=64):
    """
    Similar to compute_patch_embeddings but for text probes

    Parameters
    ----------
    model : Which Image model to use - ViNet, 
    device : hardware device: cpu/CUDA-supported gpu
    probe_list : vector list of probes (strings) 
    batch_size : Number of probes to consider at a time

    Returns
    -------
    Numpified embeddings (N probes x M features)

    """
    embs = []
    for i in range(0, len(probe_list), batch_size):
        toks = clip.tokenize(probe_list[i:i+batch_size]).to(device)
        with torch.no_grad():
            t_emb = model.encode_text(toks).cpu()
            t_emb = t_emb / t_emb.norm(dim=1, keepdim=True)
            embs.append(t_emb)
    return torch.cat(embs, dim=0).numpy()


device = "cuda" if torch.cuda.is_available() else "cpu"
model, preprocess = clip.load("ViT-L/14", device=device)


base_dir = r"C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\Projects\P1-CondProbs\Core_scripts\Cellular_targets_pipeline\Marina\CLIP"

# Patch root directory
img_dir = os.path.join(base_dir, "patches")
image_paths = sorted(glob.glob(os.path.join(img_dir, "*.png")))

# Compute patch embeddings
image_embs = compute_patch_embeddings(model, preprocess, 
                                      image_paths, 
                                      device, batch_size=64)
#image_embs_np = image_embs.numpy()

# Text probes
probes = [
    # Social gaze and attention
    "toddler looking at face",
    "toddler looking at object",
    "toddler making eye contact",
    "toddler looking away",
    "toddler following gaze",
    
    # Interaction and play
    "toddler playing alone",
    "toddlers playing together",
    "toddler showing object",
    
    # Emotion and expression
    "toddler smiling",
    "toddler crying",
    "toddler frowning",
    "toddler surprised",
    "toddler neutral face",
    "toddler happy face",
    # Motion and salience
    "moving object",
    "moving person",
    "bright object",
    "dark object",
    "hand reaching",
    
    # Social vs nonsocial contrast
    "face close-up",
    "body close-up",
    "object close-up",
    "single toddler center",
    "voice off-screen"
]

text_embs = compute_text_embeddings(model, device, 
                                       probes, batch_size=64)

# Compute similarity betweeb text and image probes
similarity = cosine_similarity(image_embs, 
                               text_embs)

# UMAP dimensionality reduction
reducer = UMAP(n_neighbors=5, min_dist=0.1, 
               metric="cosine", random_state=42, n_components=3)

embeddings_probe = reducer.fit_transform(similarity)

# Convert to dataframe
similarity_df = pd.DataFrame(similarity, columns=[f"{t}" for t in probes])



file_ids = [os.path.basename(p).replace("-", "_").replace(".png", "") for p in image_paths]
clip_target = np.array([f.split("_") for f in file_ids], dtype=float)

# Attach data and x-y-z embedding positions
df = pd.DataFrame({
    "clip_id": clip_target[:, 0],
    "target_id": clip_target[:, 1],
    "frame": clip_target[:,2],
    "x": embeddings_probe[:, 0],
    "y": embeddings_probe[:, 1],
    "z": embeddings_probe[:, 2],
})
df = pd.concat([df, similarity_df], axis=1)

# Save the data
csv_path = os.path.join(base_dir, "probe_umap.csv")
df.to_csv(csv_path, index=False, float_format="%.6f")
