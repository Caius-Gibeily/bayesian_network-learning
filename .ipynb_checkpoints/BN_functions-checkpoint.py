# -*- coding: utf-8 -*-
"""
@author: CGIBEIL

Bayesian Network Functions - to complement 
"""


import numpy as np
import pandas as pd
import pickle as pc
import string
import matplotlib.pyplot as plt
from pgmpy.models import DiscreteBayesianNetwork
import itertools
# Epoch labels
epochs = list(string.ascii_uppercase)
#%% Pickle functions 
def save_pickle(savedir,var):
    """
    Parameters
    ----------
    savedir : base directory to save pickle (.pkl)
    var : object to save

    Returns
    -------
    None.
    """
    
    with open(f"{savedir}.pkl", "wb") as f:
        pc.dump(var, f, protocol=pc.HIGHEST_PROTOCOL)
        
        
def load_pickle(datadir):
    """
    Parameters
    ----------
    datadir : directory + filename to load from

    Returns
    -------
    var : variable to assign unpickled file to 

    """
    with open(f"{datadir}.pkl","rb") as file:
        var = pc.load(file)
    return var
        
#%% Translation functions
def adj2DAG(adj):
    """

    Parameters
    ----------
    adj : input adjacency table (1 = edge, 
                                 0 = no edge)

    Returns
    -------
    dag : converted directed acyclic graph (dag) object

    """
    from pgmpy.base import DAG
    import string
    epochs = list(string.ascii_uppercase)[0:adj.shape[0]]

    edges = [(parent,child) for p,parent in enumerate(epochs) 
             for c,child in enumerate(epochs) if adj[p, c] == 1]
    dag = DAG()
    dag.add_nodes_from(epochs)  
    dag.add_edges_from(edges)
    return dag
    
def DAG2adj(model):    
    """
    Function: 
        Performs the inverse of the adj2DAG, taking a DAG object and converting it 
        into a DAG
    Parameters
    ----------
    model : network model (pgmpy.base.DAG or 
                           pgmpy.models.DiscreteBayesianNetwork)
    Returns
    -------
    adj : 

    """
    adj = np.zeros([len(model.nodes),len(model.nodes)])
    epochs = pd.Index(model.nodes)
    
    # Where 1 indicates an edge between two variables 
    for i,epoch in enumerate(model.edges()):
      
        adj[(epochs.get_loc(epoch[0]),
             epochs.get_loc(epoch[1]))] = 1                            
    return adj


#%% Bootstrap functions
def get_bootstrapped_samples(data,nboot=1000,thresh=0.54):
    """
    Parameters
    ----------
    data : input state-matrix data (Pandas object)
    nboot : optional, number of bootstraps to run
    thresh : edge frequency threshold
    Returns
    -------
    adjAll : thresholded matrices
    adjs : bootstrapped adjacency matrices

    """
    nparts = data.shape[0]
    adjs = {}
    adjAll = np.zeros([data.shape[1],data.shape[1]])
    for i in range(nboot):
        bootID = np.random.choice(nparts,nparts,replace=True)
        sampleboot = data.iloc[bootID,:]
        
        _,adjs[i],_ = learn_structure(sampleboot,"AICScore")
    
    adjAll = threshold_bootstrapAdj(adjs, thresh)
    return adjAll,adjs

def bootstrap_clips(rootdir, clips,nboot,thresh,score="AICScore"):
   """
    Parameters
    ----------
    datdir : directory and root file name

    without the  
    
    clips : list of movie clips
    
    nboot : number of bootstraps to run
    
    thresh : the threshold to apply on frequencies of edges across bootstrapped 
    networks
    
    score : Scoring method to use for structure learning. The default is "AICScore".

    Returns
    -------
    adjBootTD : All n bootstrapped adjacency matrices 
    for TD group (returned as dict)
    
    adjBootASD : All n bootstrapped adjacency matrices 
    for ASD group (returned as dict)

    """
   adjBootTD = {}
   adjBootASD = {}
   for i,m in enumerate(clips):
        adjBootTD[m] = {}
        adjBootASD[m] = {}
        
        arrayTD =  pd.read_csv(f"{rootdir}TD_clip{str(m)}.csv")
        
        _, adjBootTD[m] = get_bootstrapped_samples(arrayTD,nboot,thresh)
        
        arrayASD =  pd.read_csv(f"{rootdir}ASD_clip{str(m)}.csv")
        _, adjBootASD[m] = get_bootstrapped_samples(arrayASD,nboot,thresh)

   return adjBootTD, adjBootASD 



def threshold_bootstrapAdj(adjBoots,thresh):
    """
    Parameters
    ----------
    adjBoots : input dictionary of bootstrapped adjacency matrices
    thresh : edge frequency cutoff to use

    Returns
    -------
    adjAll : thresholded bootstrapped adjacency matrices
    """
    adjAll = np.zeros([adjBoots[0].shape[0],adjBoots[0].shape[1]])
    for i, adj in enumerate(adjBoots.values()):
        adjAll = adjAll + adj
    adjAll = adjAll/len(adjBoots)
    adjAll = adjAll>=thresh
    return adjAll

#%% Supplementary figure 7: edge cutoff validation
def collate_data(datadir,dir_adjBoot,clips,thresh):
    """
    Parameters
    ----------
    datadir : data directory of 
    state-epoch matrices
    clips : 
    cutoff : TYPE
        DESCRIPTION.

    Returns
    -------
    dat : TYPE
        DESCRIPTION.

    """
    dat = []
    
    
    adjBootTD = load_pickle(f"{dir_adjBoot}TD")
    adjBootASD = load_pickle(f"{dir_adjBoot}ASD")

    adjAllTD  =  threshold_bootstrapAdj(adjBootTD,thresh)
    adjAllASD =  threshold_bootstrapAdj(adjBootASD,thresh)
    
    adjAll = {"TD":adjAllTD,
              "ASD":adjAllASD}
    
    for i, clip in enumerate(clips):
        adjBootTD = threshold_bootstrapAdj(adjBootTD[clip],thresh)
        
        for dx in ["TD","ASD"]:
            data =  pd.read_csv(f"{datadir}{dx}_clip{str(clip)}",header=False)
            data.columns = epochs[0:data.shape[1]]
            
            DAG_thresh = adj2DAG(adjAll[dx][clip])
            model = learn_parameters(data, DAG_thresh)
            
            dat.append({"clip": clip, "dx": dx, "model": model, 
                        "adj": adjAll[dx][clip], "data": data })
            
         
    dat = pd.DataFrame(dat)
    return dat



#%% NMI calculation and plotting
def calc_nmi(parent,child):
    """
    Parameters
    ----------
    parent : column of states in the upstream epoch
    child :  column of states in the downstream epoch

    Returns
    -------
    The pairwise normalised mutual information (NMI) between upstream and 
    downstream epochs
    """
    
    states_p, counts_p = np.unique(parent, return_counts=True)
    p_p = counts_p / counts_p.sum()
    
    states_c, counts_c = np.unique(child, return_counts=True)
    p_c = counts_c / counts_c.sum()
    
    # Mutual Information
    mi = 0.0
    for i, p_pi in enumerate(p_p):
        for j, p_cj in enumerate(p_c):
            joint = np.mean(
                (parent == states_p[i]) &
                (child == states_c[j])
            )
            if joint > 0:
                mi += joint * np.log2(joint / (p_pi * p_cj))
    
    # Entropies
    H_p = -np.sum(p_p * np.log2(np.clip(p_p, 1e-12, 1)))
    H_c = -np.sum(p_c * np.log2(np.clip(p_c, 1e-12, 1)))

    # Normalised MI
    denom = (H_p + H_c)

    if denom == 0:
        return 0.0
    return 2 * mi / denom


def calc_mi_model(data,model):
    """
    Parameters
    ----------
    data : input state-matrix dataframe (Pandas object)
    model : model graph structure (either a pgmpy.base.DAG or 
                                   pgmpy.models.DiscreteBayesianNetwork)

    Returns
    -------
    mis : all NMIs for each pair of epochs
    """
    # iterate through each node pair (parent_node,child_node)
    mis = {}
    for p, pnode in enumerate(model.nodes):
        for c,cnode in enumerate(model.nodes):
            if p == c:
                continue
            mis[(pnode,cnode)] = calc_nmi(data[pnode],data[cnode])
    return mis
    
def find_ci_edges(model,data):
    
    """
    Function: 
        The recovered BN structure makes a statement regarding which variables 
        are conditionally separated (d-separation criteria). If a conditional dependency 
        between two unconnected variables is explained by the network structure, the
        two variables are said to be conditionally independent. Otherwise, the two variables are said to 
        have a conditional dependency not explained by the model. 
        
    Parameters
    ----------
    model : BN model of type pgmpy.models.DiscreteBayesianNetwork
    data : state-matrix 

    Returns
    -------
    ci_array : array of implied pairwise conditional independencies between epochs.
    Where an unexplained conditional dependency arises, this is marked as True
    """
    
    from pgmpy.metrics import implied_cis
    from pgmpy.estimators.CITests import chi_square
    
    cis = implied_cis(model,data,ci_test=chi_square)
    cis = cis[cis["p-value"]<0.05] 
    ci_array = np.zeros([len(model.nodes),len(model.nodes)])
    ci_array[cis["u"].to_numpy().astype(int),
             cis["v"].to_numpy().astype(int)] = 1
    return ci_array
    

def plot_mi_adjacency(model, adj,ci_array, mi_dict, 
                      bootstrapped,cmap = "Blues", fig=None,ax=None):
    """
    Function: 
        We have derived (1) bootstrapped graphs (bootstrapped), (2) an edge frequency 
        threshold to incorporate edges that recur across more than the proportion 
        of bootstrapped samples set by the threshold, (3) the normalised mutual information 
        (NMI) and implied conditional indepedendencies between each pair of variables
        
        This function then plots these pieces of information together
    Parameters
    ----------
    model : Bayesian Network model (pgmpy.models.DiscreteBayesianNetwork)
    adj : Adjacency matrix
    ci_array : Implied conditional independencies (1 = implied independence given structure violated)
    mi_dict : Collection of normalised mutual information (NMI) values calculated per pair of epochs
    bootstrapped : Bootstrapped networks
    cmap : TYPE, optional
        "Blues" for TD group, "Greens" for ASD group
    fig : fig on which to plot
    ax : figure axes

    Returns
    -------
    None.

    """
    from matplotlib import patches
    import string
    
    nodes = list(model.nodes())
    n = len(nodes)
    adj_matrix = np.zeros((n, n))

    # Fill adjacency matrix with MI values
    for (parent, child), mi in mi_dict.items():
        if parent in nodes and child in nodes:
            i, j = nodes.index(parent), nodes.index(child)
            adj_matrix[i, j] = mi
    
    bootAdj = sum(bootstrapped.values())/len(bootstrapped)
    print(bootAdj)
    
    if fig == None:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.set_aspect("equal")

    # Coordinates
    x, y = np.meshgrid(range(n), range(n))
    x, y, mi_values = x.flatten(), y.flatten(), adj_matrix.flatten()

    # Keep only upper triangle (lower triangle contains forbidden edges)
    mask = np.triu(np.ones((n, n)), k=0).flatten().astype(bool)
    x, y, mi_values = x[mask], y[mask], mi_values[mask]
    
    # bootstrapped edge frequencies
    bs_values = bootAdj.flatten()
    bs_values = bs_values[mask] 

    # Normalise circle sizes
    if np.nanmax(mi_values) > 0:
        sizes = 800 * (mi_values / np.nanmax(mi_values))
    else:
        sizes = np.zeros_like(mi_values)

    # Plot circles
    sc = ax.scatter(
        x,
        y,
        s=sizes,
        c=bs_values,
        cmap=cmap,
        edgecolors="none",
        alpha=0.9
    )

    # Plot shaded areas corresponding to incorporated edges at the selected cutoff
    parents, children = adj.nonzero()
 
    # Draw light gridlines for upper triangle
    for i in range(n):
        for j in range(i, n):
            ax.plot([j - 0.5, j + 0.5], [i - 0.5, i - 0.5], color="lightgray", lw=0.8)
            ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color="lightgray", lw=0.8)
            ax.plot([j - 0.5, j - 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.8)
            ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.8)

    # Hide lower triangle with white polygon
    triangle = plt.Polygon(
        [[-0.5, n - 0.5], [-0.5, -0.5], [n - 0.5, n - 0.5]],
        closed=True,
        facecolor="white",
        edgecolor="none",
        zorder=3,
    )
    ax.add_patch(triangle)

    # Add uniform bold border enclosing upper triangle
    border = patches.Polygon(
        [[-0.5, -0.5], [n - 0.5, -0.5], [n - 0.5, n - 0.5]],
        closed=False,
        fill=False,
        edgecolor="black",
        linewidth=2.5,
        capstyle='projecting',
        joinstyle='miter',
        zorder=4,
    )
    ax.add_patch(border)

    # Add bold diagonal
    ax.plot(
        [-0.5, n - 0.5],
        [-0.5, n - 0.5],
        color="black",
        lw=2.5,
        solid_capstyle="projecting",
        zorder=4,
    )

    
    # Ticks only on top
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    
    
    ax.set_xticklabels(epochs[:n], fontsize=35)
    ax.set_yticklabels(epochs[:n],fontsize=35)
    ax.xaxis.tick_top()

    # Clean up
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.invert_yaxis()

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
    sc.set_clim(vmin=0, vmax=1)
    cbar.ax.tick_params(labelsize=20)
    
    rgba = np.zeros((*adj.shape, 4))
    rgba[adj == 1] = [1.0, 0.75, 0.8, 1.0]
    rgba[ci_array==1] = [0.42,0.65,0.73,1]
    ax.imshow(
        rgba,
        alpha=0.5
        )
    
    
#%% Network learning
def generate_forbidden_edges(columns):
    """
    Parameters
    ----------
    columns : column names of state-matrix. Epochs are labelled alphabetically 

    Returns
    -------
    forbidden : TYPE
        DESCRIPTION.

    """
    forbidden = []
    for i in range(len(columns)):
        for j in range(i+1, len(columns)):
            forbidden.append((columns[j], columns[i]))  # No future -> past edges
    return forbidden

def learn_structure(data,score="AICScore"): 
    """
    Parameters
    ----------
    data : input data (as str if directory or as numpy array/Pandas dataframe)
    score : TYPE, The default structure learning score to use. The default is 
    Akaike Information Criterion (AIC)

    Returns
    -------
    best_model : the learned structure
    adj : recovered adjacency matrix
    data : mirror output of data

    """
    
    from pgmpy.estimators import BIC,HillClimbSearch,AIC,BDeu
    from pgmpy.estimators import ExpertKnowledge

    if isinstance(data,str): # if data input is a directory
        data = pd.read_csv(data+".csv",header=None )
    else:
        data = pd.DataFrame(data).astype(str)

    epochs = list(string.ascii_uppercase)
    data.columns  = epochs[0:data.shape[1]]

    # Restrict structures to only 
    black_list= generate_forbidden_edges(data.columns)
    forbidden = ExpertKnowledge(forbidden_edges=black_list)
    if score == "BICScore":
        scoring_method = BIC(data)
    elif score == "AICScore":
        scoring_method = AIC(data)
    elif score == "BDeu":
        scoring_method = BDeu(data)

    state_names={data.columns[i]:np.unique(data.iloc[:,i]) for i in range(data.shape[1])}
    hc = HillClimbSearch(data,state_names=state_names)
    best_model = hc.estimate(scoring_method=scoring_method,
                             expert_knowledge=forbidden,show_progress=True) # ,tabu_length=10
    
    adj = DAG2adj(best_model)

    return best_model, adj,data


def learn_parameters(data,dag): 
    """
    Parameters
    ----------
    data : state-epoch matrix
    dag : directed acyclic graph object

    Returns
    -------
    model : model graph and parameters (pgmpy.models.DiscreteBayesianNetwork)
    """
    
    from pgmpy.estimators import MaximumLikelihoodEstimator, ExpectationMaximization
    
    model = DiscreteBayesianNetwork(dag)
    model.fit(data, estimator=MaximumLikelihoodEstimator)

    return model


#%% Model evaluation
def js_distance(p,q):
    from scipy.spatial import distance
    if np.isnan(p).any() or np.isnan(q).any():
        return np.nan
    else: 
        return distance.jensenshannon(p, q)

def compare_recovered_ppds(model1,model2, fun=np.mean):
    from pgmpy.inference import CausalInference
    variables = model1.nodes()
    variables2 = model2.nodes()
    if len(variables) != len(variables2):
        raise ValueError("Unequal number of variables")
    #model1Inf = VariableElimination(model1)
    #model2Inf = VariableElimination(model2)
    model1Inf = CausalInference(model1)
    model2Inf = CausalInference(model2)
    jsds = []
    for i, node in enumerate(variables):
        parents = model2.get_parents(node)
        par_states = [model2.states[k] for k in variables if k in parents]
        combs = itertools.product(*par_states)
        
        for c, comb in enumerate(combs):

            evidence = {parents[p]: state_val for p,state_val in enumerate(comb)}
            print(evidence)
            query_fit = model2Inf.query(variables=[node],do=evidence)
            query_ground = model1Inf.query(variables=[node],do=evidence)
            diff = len(query_fit.values) - len(query_ground.values)
            if diff == 0:   
                jsds.append(js_distance(query_fit.values, query_ground.values))
            elif diff < 0:
                query_fit.values = np.append(query_fit.values,np.repeat(0,abs(diff)))
                jsds.append(js_distance(query_fit.values, query_ground.values))
            elif diff > 0:
                query_ground.values = np.append(query_ground.values,np.repeat(0,abs(diff)))
                jsds.append(js_distance(query_fit.values, query_ground.values))                
    jsd_fn = fun(jsds)
    return jsd_fn


def hamming_distance(adj1,adj2,directed=True):

    if not directed:
        T = np.triu(adj1 + adj1.T, 1)
        P = np.triu(adj2 + adj2.T, 1)
        diff = (T != P).astype(int)
        possible = T.size // 2
        raw = diff.sum()
    else:
        diff = (adj1 != adj2).astype(int)
        raw = diff.sum()
        possible = adj1.size
    return raw / possible

#%% Model prediction

def get_all_predictions(clips,data,nreps):
    
    """
    Function: 
        Collate predictive accuracies for each query participant across all epochs
        under three different models:
            M0 - null model, where predictions are drawn from a posterior categorical 
            distribution informed by state frequencies in each epoch
            M1 - first-order only model, predictions are drawn from a reduced model 
            containing only the first-order edges of each BN
            M2 - full model, containing all edges learned via the structure learning
            and bootstrapping steps
        The posterior probabilities of entering each state are then inferred 
        for each participant.
        The function using parallel programming to accelerate the runtime of the 
        function
    Downstream: 
        Bayes factors are computed between M1 and M0 and M2 and M0, providing the
        relative evidence favouring 
        the full (M2) or partial (M1) models in predicting subsequent states over
        the null model.
        
    Parameters
    ----------
    clips : List of movie clips to perform the prediction over
    data : dataframe object collating state-matrix, adjacency matrix, 
    model parameters, dx and clip id (Pandas object)
    
    nreps : Number of likelihood weighted samples to take

    Returns
    -------
    preds: A pandas dataframe of epoch-wise prediction probabilities for each clip
    and participant who viewed that clip 

    """
    from tqdm import tqdm
    preds = []
    for i,clip in enumerate(clips):
        
        it = tqdm(range(nreps), desc="Likelihood-weighted sampling")
        for rep in it:
            seed = np.random.randint(0,10000)
            
            # TD #################################################################
            dat_loop = data[(data["dx"]=="TD") & (data["clip"]==clip)]
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"evidence", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "evidence",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"rand", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "rand",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"forder", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "forder",
                                                   seed = seed),
                          "seed": seed})
          
            
            # AD ############################################################
            dat_loop = data[(data["dx"]=="ASD") & (data["clip"]==clip)]
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"evidence", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "evidence",
                                                   seed = seed),
                          "seed": seed})
            
            
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"rand", 
                          "pred_array":  predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "rand",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"forder", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "forder",
                                                   seed = seed),
                          "seed": seed})
            
        
            

            #################################################################
    preds = pd.DataFrame(preds)
    return preds
        

def predictBayes(data, adj, model, size=10000, n_jobs=-1,
            cond = "evidence",alpha=1,seed = None):
    """
    Parameters
    ----------
    data : state-matrix dataframe object 
    adj : adjacency matrix 
    model : the BN model for which prediction is performed
    size : size of likelihood-weighted sample to pull for each prediction
    n_jobs : optional, specifies the number of computing cores to use
    cond : optional, which model setup to use, M0, M1 or M2
    alpha : the prior set on each state. By default = 1
    seed : optional, seed is set `get_all_predictions` for reproducibility

    Returns
    -------
    pred : set of all predictions for the query movie clip and group
    """
    from pgmpy.factors.discrete import State
    from joblib import Parallel, delayed
    from pgmpy.sampling import BayesianModelSampling
    

    pred = np.zeros([data.shape[0], data.shape[1]])
    model_preds= []
    if cond == "forder": # M1 condition with only first-order edges
        # Convert model adjacency matrix into a numpy array
        adj = DAG2adj(model)
        # Extract only the first-order edges from the model
        adj[~np.eye(adj.shape[0], dtype=bool,k=1)] = 0
        # Convert adjacency matrix back into DAG
        dag = adj2DAG(adj)

        # LOO parameter learning. Using the same bootstrapped structure, 
        # refit the model parameters leaving each participant out a time 
        for i in range(data.shape[0]): # for each participant
            
            data_train = data.drop(i) # drop that participant from the training 
            model_preds.append(learn_parameters(data_train, dag))

        
        
    elif cond == "evidence" or cond == "rand": # M2 or M0 conditions

        for i in range(data.shape[0]):
            data_train = data.drop(i)
            model_preds.append(learn_parameters(data_train, adj2DAG(adj)))

    def predict_row(i, samp, model_rec, adj_rec,data_train): # predict each 
        # participant's state sequence
        

        rng = np.random.default_rng(seed)
        
        inference = BayesianModelSampling(model_rec)
        
        row_pred = np.zeros(data.shape[1])
        
        epochs = data_train.columns

        for j,node in  enumerate(model_rec.nodes): # For each epoch, 

            if cond == "rand": # For M0, find the frequencies of all the states in the epoch
                states,prior_dist = np.unique(data_train.iloc[:,j],
                                          return_counts= True)
                #Specify alpha as the prior set on the null model for each state 
                #Dirichlet prior is the conjugate prior of the categorical distribution. The 
                #Updated probabilities are given for each state i as 
                #n_{state,i} + alpha / sum(n_{states}) + (alpha)
                
                null_posterior = (prior_dist+alpha)/(np.sum(prior_dist) + states.size * alpha)
                
                # which state did the participant enter
                
                # log-likelihood
                samp_state = str(samp.iloc[j])
                row_pred[j] = np.log10(null_posterior[states==samp_state]+1e-16) # 1e-16 added to avoid log0
 
            
            else: # if not M0 (M1/M2)
                # for the first epoch with no history
                if j == 0: #len(parents) == 0:
                    
                    # generate a simulated dataframe given the model parameters.
                    df_sim = inference.likelihood_weighted_sample(
                         size=size, show_progress=False
                    )
                
                else:
                    # for subsequent epochs
                    history = epochs[0:j]
                    parent_states = samp.iloc[:]
                    evidence = [
                        State(var=str(pid), state=parent_states.iloc[n])
                        for n, pid in enumerate(history)
                    ]
                     
                    # generate likelihood-weighted samples after setting the participant's 
                    # state history 
    
                    df_sim = inference.likelihood_weighted_sample(
                         size=size, show_progress=False, evidence=evidence
                    )
                    
               
                samples = df_sim.iloc[:, j].astype(str)
                weights = df_sim["_weight"].values
                
                post_dist = {}
                for s, w in zip(samples, weights):
                    post_dist[s] = post_dist.get(s, 0.0) + w

                
                Z = sum(post_dist.values())
                if Z <= 0 or not np.isfinite(Z):
                    # Default to null if likelihoods low
                    states,prior_dist = np.unique(data_train.iloc[:,j],
                                              return_counts= True)
                    # as before for M0
                    null_posterior = (prior_dist+alpha)/(np.sum(prior_dist) + 
                                                         states.size * alpha)
                    
                    row_pred[j] = np.log10(null_posterior[states==samp_state]+1e-16)
                    
                    continue                
                
                else:
                    for s in post_dist:
                        post_dist[s] /= Z
                
                samp_state = str(samp.iloc[j])
                prob = post_dist.get(samp_state, 0.0)
  
                row_pred[j] = np.log10(prob+1e-16)
                
        return i, row_pred

    # Parallel processing 
    results = Parallel(n_jobs=n_jobs, verbose=False)(
        delayed(predict_row)(
            i, samp, model_preds[i], adj, data_train)
        
    for i, samp in data.iterrows()
    )

    for i, row_pred in results:

        pred[i, :] = row_pred

    return pred



