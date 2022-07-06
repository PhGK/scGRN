import pandas as pd
import torch as tc
import numpy as np
import numpy.ma as ma

def load_data(name):
    normalize_data = False



    if name == 'epi2000_bih':
        data = pd.read_csv('./data/epi_top2000.csv')
        data = data[(data['cell_type_epi'] != 'Neuroendocrine') & (data['patient_id'] != 'p028') & (data['patient_id'] != 'p029')] #filter patients without tumor tissue


        cell_type = data['cell_type_epi']
        patient_id = data['patient_id']



        all_patient_ids = np.array(patient_id.unique())
        all_patient_ids.sort()
        all_cell_types = np.array(cell_type.unique())
        all_cell_types.sort()

        print(all_cell_types, all_patient_ids)

        cell_onehot = tc.tensor(np.array(cell_type)[:,None] == all_cell_types[None,:])  
        patient_onehot = tc.tensor(np.array(patient_id)[:,None] == all_patient_ids[None,:])  

        patient_cell_onehot = tc.cat((patient_onehot, cell_onehot), axis = 1)

        print(patient_cell_onehot.shape)
        nsamples = data.shape[0]
        
        sel_genes = np.array(pd.read_csv('./data/selected_genes.csv', header=None)) # select only genes that are protein coding

        sel_columns = [gene for gene in data.columns if gene in sel_genes]
                

        sc_data = data[sel_columns].iloc[:,:800]    #select genes that have the highest expression (genes are already sorted from left to right)
        
        feature_names = list(sc_data.columns) + list(all_patient_ids) + list(all_cell_types)
        nfeatures = len(feature_names)
        sample_names = np.array(list(data.iloc[:, 1]))
        np.random.seed(0)
        permutation = np.random.permutation(sc_data.shape[0]).astype(int)
        randomized_data = tc.tensor(np.array(sc_data))[permutation, :]

        patient_cell_shuffled = patient_cell_onehot[permutation,:]         
        sample_names = sample_names[permutation]
        normalize_data = True




    if name == 'random':
        randomized_data = tc.randn(1000,10)
        sample_names = ['sample' + str(i) for i in range(1000)]
        feature_names = ['feature' + str(i) for i in range(10)]

    if normalize_data:
        meanv, sdv, minv, maxv = randomized_data.mean(axis=0, keepdim=True), randomized_data.std(axis=0, keepdim=True), \
                                 randomized_data.min(axis=0)[0], tc.abs(randomized_data).max(axis=0)[0]
        


        randomized_data = (randomized_data) / sdv

    train_set, test_set = randomized_data[:randomized_data.shape[0]//10*9,:], randomized_data[randomized_data.shape[0]//10*9:,:]
    train_names, test_names = sample_names[:randomized_data.shape[0]//10*9], sample_names[randomized_data.shape[0]//10*9]

    train_onehot, test_onehot = patient_cell_shuffled[:randomized_data.shape[0]//10*9,:], patient_cell_shuffled[randomized_data.shape[0]//10*9:,:]

    print(train_set.shape, test_set.shape)
    return train_set.float(), test_set.float(), feature_names, train_names, test_names, train_onehot, test_onehot




load_data('epi2000_bih')
