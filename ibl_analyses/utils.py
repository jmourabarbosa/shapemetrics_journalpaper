# %%
# -*- coding: utf-8 -*-
"""
@author: Amin
"""

import numpy as np

import ray
from netrep.metrics import LinearMetric
from netrep.metrics import GaussianStochasticMetric

from scipy.stats import rankdata
from sklearn.metrics import pairwise_distances

import jax

# %%
@ray.remote
def _stochastic_metrics_pair(pair,alpha=2.,niter=1000):
    Xi,Xj = pair

    metric = GaussianStochasticMetric(alpha,niter=niter)
    metric.fit(Xi,Xj)
    dist_neural = metric.score(Xi,Xj)
    return dist_neural

# %%
def ssd(pairs,alpha=2.,niter=1000):
    refs = [
        _stochastic_metrics_pair.remote(pair,alpha,niter)
            for pair in pairs
    ]
    D = ray.get(refs)
    D = np.array(D)

    return D

# %%
@ray.remote
def _deterministic_metrics_pair(pair,alpha=0.):
    Xi,Xj = pair

    metric = LinearMetric(alpha=alpha,center_columns=True,score_method='euclidean')
    metric.fit(Xi,Xj)
    dist_neural = metric.score(Xi,Xj)
    return dist_neural

# %%
def dsd(pairs,alpha=0.):
    refs = [
        _deterministic_metrics_pair.remote(pair,alpha)
            for pair in pairs
    ]
    D = ray.get(refs)
    D = np.array(D)

    return D

# %%
def delay_embedding(data,tau,delta=1):
    K,C,T,N = data.shape
    embedding = np.zeros((K,C,T-(tau-1)*delta,N*tau))

    for d in range(tau):
        embedding[:,:,:,d*N:d*N+N] = data[:,:,(tau-1-d)*delta:T-d*delta]

    return embedding

# %%
def create_adjacency(x):
    idx = rankdata(x, method='dense',axis=0)-1
    dist = pairwise_distances(idx,metric='l1')
    dist[dist != 1] = 0
    return dist

# %%
def split_data_cv(data,props,seeds):
    # props: train, validation, test
    # seeds: test, validation
    # data: y (possibly mu, sigma, F, mu_g, sigma_g)

    assert 'train' in props.keys() and 'test' in props.keys() and 'validation' in props.keys()
    assert props['train'] + props['test'] + props['validation'] == 1
    assert 'test' in seeds.keys() and 'validation' in seeds.keys()
    assert 'y' in data.keys()
     
    N,M,D = data['y'].shape
    
    trial_indices = jax.random.permutation(
        jax.random.PRNGKey(seeds['test']),
        np.arange(N)
    )

    test_trials = trial_indices[-int(props['test']*N):]

    train_trials = jax.random.choice(
        jax.random.PRNGKey(seeds['validation']),
        shape=(int(N*props['train']),),
        a=trial_indices[:-int(props['test']*N)],
        replace=False
    ).sort()

    validation_trials = np.setdiff1d(trial_indices[:-int(props['test']*N)],train_trials).tolist()

    out = {}
    for k in data.keys():
        out[k+'_train'] = data[k][train_trials,...]
        out[k+'_test'] = data[k][test_trials,...]
        out[k+'_validation'] = data[k][validation_trials,...]
    
    return out