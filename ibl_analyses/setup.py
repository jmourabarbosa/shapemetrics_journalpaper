
good_eids = [
    'e2b845a1-e313-4a08-bc61-a5f662ed295e',
    #'111c1762-7908-47e0-9f40-2f2ee55b6505',
    '2bdf206a-820f-402f-920a-9e86cd5388a4',
    '5dcee0eb-b34d-4652-acc3-d10afc6eae68',
    'c7bf2d49-4937-4597-b307-9f39cb1c7b16',
    '824cf03d-4012-4ab1-b499-c83a92c5589e',
    '51e53aff-1d5d-4182-a684-aba783d50ae5',
    'c51f34d8-42f6-4c9c-bb5b-669fd9c42cd9',
    '0802ced5-33a3-405e-8336-b65ebc5cb07c',
    'd2832a38-27f6-452d-91d6-af72d794136c',
    '88224abb-5746-431f-9c17-17d7ef806e6a',
    '72cb5550-43b4-4ef0-add5-e4adfdfb5e02',
    '0a018f12-ee06-4b11-97aa-bbbff5448e9f',
    '7af49c00-63dd-4fed-b2e0-1b3bd945b20b',
    '73918ae1-e4fd-4c18-b132-00cb555b1ad2',
    'd0ea3148-948d-4817-94f8-dcaf2342bbbe',
    'c4432264-e1ae-446f-8a07-6280abade813',
    '746d1902-fa59-4cab-b0aa-013be36060d5',
    '54238fd6-d2d0-4408-b1a9-d19d24fd29ce',
    '4a45c8ba-db6f-4f11-9403-56e06a33dfa4',
    '754b74d5-7a06-4004-ae0c-72a10b6ed2e6',
    'd23a44ef-1402-4ed7-97f5-47e9a7a504d9',
    '15763234-d21e-491f-a01b-1238eb96d389',
    #'71e55bfe-5a3a-4cba-bdc7-f085140d798e', # load_session_data() breaks
    '7f6b86f9-879a-4ea2-8531-294a221af5d0',
    '61e11a11-ab65-48fb-ae08-3cb80662e5d6',
    'c7248e09-8c0d-40f2-9eb4-700a8973d8c8',
    'f312aaec-3b6f-44b3-86b4-3a0c119c0438',
    'dda5fc59-f09a-4256-9fb5-66c67667a466',
    'db4df448-e449-4a6f-a0e7-288711e7a75a',
    'ecb5520d-1358-434c-95ec-93687ecd1396',
    'ee40aece-cffd-4edb-a4b6-155f158c666a',
    'b03fbc44-3d8e-4a6c-8a50-5ea3498568e0',
    'a8a8af78-16de-4841-ab07-fde4b5281a03',
    # 'e9b57a5a-b06d-476d-ad20-7ec42a16f5f5',
    '4b7fbad4-f6de-43b4-9b15-c7c7ef44db4b'
]

params = {
    'file': '../data/',
    'tag': '2022_Q2_IBL_et_al_RepeatedSite',
    'probe': 'probe00',
    'sessions': [0,5,6],
    'areas': ['CA1','DG','LP','PO','VISa'],
    'props':{'train':.5,'test':.5,'validation':0},
    'seeds':{'train':0,'test':1,'validation':2},
    'n_neurons': None, # all neurons
    'n_trials': None, # all trials
    'pre_time':0,
    'post_time':.4,
    'n_bins': 10,
    'align_to': 'response',
    'train_trial_prop':.9, 
    'train_condition_prop':1, 
    'seed':0,
    'verbose': True
}