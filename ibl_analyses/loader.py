# %%
# -*- coding: utf-8 -*-
"""
@author: Amin
"""

from ibllib.atlas.regions import BrainRegions
from ibllib.atlas import AllenAtlas

from brainbox.io.one import SpikeSortingLoader

from one.api import ONE

import brainbox

from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict
import warnings
import ray

import utils

warnings.filterwarnings(
    "ignore", 
    category=DeprecationWarning
) 

import yaml
import sys

sys.path.insert(1,'../')

# %%
class IBLSession:
    def __init__(self,params):
        # params: file,tag,probe,areas,eid,
        # align_to='response',pre_time=0,post_time=.2,n_bins
        # train_trial_prop, train_condition_prop, seed
        
        self.eid = params['eid']
        self.params = params

        cache_dir = Path(params['file'])
        cache_dir = cache_dir / params['tag']
        self.one = ONE(
            base_url="https://openalyx.internationalbrainlab.org", 
            username='intbrainlab',
            password="international", 
            silent=True, 
            cache_dir=cache_dir
        )
        self.brain_atlas = AllenAtlas()

        # self.load_session()
        # self.load_session_data()

    def check_rep(self,spikes,clusters):
         # spikes.clusters provides ids which match clusters['cluster_id']
        assert np.all(np.sort(np.unique(spikes.clusters)) == np.sort(np.unique(clusters['cluster_id']))) 
        assert len(clusters['cluster_id']) == len(clusters['acronym'])
        assert len(clusters['cluster_id']) == len(clusters['z'])
        assert len(clusters['cluster_id']) == len(clusters['firing_rate'])
        assert len(clusters['cluster_id']) == len(np.unique(clusters['cluster_id']))
        return


    def load_session(self):
        trials = self.one.load_object(self.eid,'trials',collection='alf')
        try:
            sl = SpikeSortingLoader(
                eid=self.eid, 
                pname=self.params['probe'], 
                one=self.one, 
                atlas=self.brain_atlas
            )
            spikes, clusters, channels = sl.load_spike_sorting()
            clusters = sl.merge_clusters(spikes, clusters, channels)
            self.check_rep(spikes, clusters)
            probe_data = dict(spikes=spikes, clusters=clusters, channels=channels)
        except:
            probe_data = []
        
        self.data = {'trials':trials,self.params['probe']:probe_data}
        return self.data
    
    
    def load_session_data(self):
        
        spikes, clusters, channels = self.data[self.params['probe']].values()

        trials = self.data['trials']

        bin_size = (self.params['post_time']-self.params['pre_time'])/self.params['n_bins']

        acronym_allen = self.data[self.params['probe']]['clusters']['acronym']
        acronym_beryl = BrainRegions().acronym2acronym(acronym_allen, 'Beryl')

        acronym_bool = np.isin(acronym_beryl, self.params['areas'])

        trial_indices = trials.probabilityLeft != .5 # exclude 50-50 trials
        

        y, t = brainbox.singlecell.bin_spikes2D(
            spike_times=spikes.times, 
            spike_clusters=spikes.clusters, 
            cluster_ids=clusters.cluster_id,
            align_times=trials.stimOn_times if self.params['align_to'] == 'stim' else trials.response_times, 
            pre_time=self.params['pre_time'],
            post_time=self.params['post_time'],
            bin_size=bin_size
        )
        
        y = y[:,acronym_bool,:]

        reaction_times = trials['response_times'] - trials['stimOn_times']
        correct = trials['feedbackType']
        
        contrast = np.diff(
            np.nan_to_num(np.c_[trials['contrastLeft'], trials['contrastRight']])
        )*100

        reaction_times = reaction_times[trial_indices]
        correct = correct[trial_indices]
        y = y[trial_indices]
        contrast = contrast[trial_indices]

        x, indices, counts = np.unique(
            contrast,axis=0,
            return_inverse=True,
            return_counts=True
        )
        n_trials = min(counts)
        y = np.array([
            [y[j][:,t_]
            for j in np.where(indices==i)[0][:n_trials].tolist()] 
            for i in range(x.shape[0]) for t_ in range(len(t))]
        ).transpose(1,0,2)
        

        reaction_times = np.array([
            [reaction_times[j]
            for j in np.where(indices==i)[0][:n_trials].tolist()] 
            for i in range(x.shape[0])]
        ).T
        correct = np.array([
            [correct[j]
            for j in np.where(indices==i)[0][:n_trials].tolist()] 
            for i in range(x.shape[0])]
        ).T
        
        
        x = np.array([[x_,t_] for x_ in x.squeeze() for t_ in t])

        y = np.sqrt(
            y[:,:,np.argsort(y.mean(0).var(0))[::-1]]
        )

        if self.params['n_trials'] is not None:
            y = y[:self.params['n_trials']]
        
        if self.params['n_neurons'] is not None:
            y = y[:,:,:self.params['n_neurons']]

        self.data = utils.split_data_cv(
            {
                'y':y,
                'x':x,
                'reaction_times':reaction_times[:self.params['n_trials'],:],
                'correct':correct[:self.params['n_trials'],:],
             },
            self.params['props'],
            self.params['seeds']
        ) 
        
    
    def load_train_data(self):
        return self.data['x_train'], self.data['y_train'], self.data['reaction_times_train'], self.data['correct_train']
    
    def load_test_data(self):
        return self.data['x_test'], self.data['y_test'], self.data['reaction_times_test'], self.data['correct_test']
    
    def load_validation_data(self):
        return self.data['x_validation'], self.data['y_validation'], self.data['reaction_times_validation'], self.data['correct_validation']

    
# %%
@ray.remote
class IBLSessionRemote(IBLSession):
    pass

# %%
class IBLDataLoader:
    def load_sessions(self):
        refs = [sess.load_session.remote() for sess in self.sessions]
        self.data = ray.get(refs)

        refs = [sess.load_session_data.remote() for sess in self.sessions]
        ray.get(refs)

    def __init__(self,params,eids=None):
        # params: file,tag,probe,sessions,areas

        cache_dir = Path(params['file'])
        cache_dir = cache_dir / params['tag']
        one = ONE(
            base_url="https://openalyx.internationalbrainlab.org", 
            username='intbrainlab',
            password="international", 
            silent=True, 
            cache_dir=cache_dir
        )

        bwm_sessions = one.alyx.rest(
            'sessions', 'list', dataset_types='spikes.times', tag=params['tag']
        )

        df = pd.DataFrame(bwm_sessions)
        if eids is None:
            eids = list(df['id']) if params['sessions'] is None else [list(df['id'])[s] for s in params['sessions']]

        self.eids = eids
        self.probe = params['probe']
        self.areas = params['areas']

        
        self.sessions = [
            IBLSessionRemote.remote(
                {**params,**{'eid':eid}}
            ) for eid in eids
        ]
        self.load_sessions()

        valid = [i for i in range(len(self.data)) if bool(self.data[i][params['probe']])]
        self.eids = [eids[i] for i in valid]
        self.sessions = [self.sessions[i] for i in valid]
        self.data = [self.data[i] for i in valid]



    def load_train_data(self):
        return zip([sess.load_train_data() for sess in self.sessions])
        
    
    def load_test_data(self):
        return zip([sess.load_test_data() for sess in self.sessions])
        
    
    def load_validation_data(self):
        return zip([sess.load_validation_data() for sess in self.sessions])
