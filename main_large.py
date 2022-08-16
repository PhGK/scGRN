from LRPmodel_cuda import train, calc_all_paths, Model
#from dataloading_simple import Dataset_train, Dataset_LRP
from dataloading_simple import Dataset_train, Dataset_LRP
from data import load_data
import torch as tc
import sys
from joblib import Parallel, delayed
import os


data_type =  sys.argv[1]
train_new_network = sys.argv[2] == 'new'
print(sys.argv[3])
sample_id = int(sys.argv[3])
PATH = sys.argv[4]
cuda = sys.argv[5] == 'cuda'
nepochs = int(sys.argv[6])
njobs = int(sys.argv[7])

model_path = PATH + '/results/models/'
result_path= PATH + '/results/LRP_values_raw/'

print('data_type:', data_type)

train_data, test_data, featurenames, train_names, test_names, train_onehot, test_onehot = load_data(data_type)
print(train_onehot.shape, test_onehot.shape)
input_size = train_data.shape[1] + train_onehot.shape[1]
if train_new_network:
    print(train_data.shape)

    model = Model(input_size*2, input_size, hidden = (train_data.shape[1])*10, hidden_depth=1) # best: hidden_depth=3, train_data.shape[1]*3
    testlosses, epoch_list, network_list= train(model, train_data, test_data, train_onehot, test_onehot, epochs=nepochs, lr = 0.01, noise_level=0.0, weight_decay = 0.0, batch_size = 5, device = tc.device("cuda:0" if cuda else "cpu"))

    if not os.path.exists(model_path):
        os.makedirs(model_path)
    #model.cpu()

    mindex = tc.argmin(testlosses)
    min_network = network_list[mindex]
   
    #tc.save(model.cpu(), model_path + data_type + '.pt')
    tc.save(min_network, model_path+data_type+ '.pt')
    print('saved model from epoch ' + str(epoch_list[mindex]) + ' with testloss ' + str(testlosses[mindex]))

else:
    print('loading old model')
    #model = tc.load(model_path + data_type + '.pt', map_location=tc.device('cpu'))
    model = Model(input_size*2, input_size, hidden = (train_data.shape[1])*10, hidden_depth=1)
    model.load_state_dict(tc.load(model_path + data_type + '.pt', map_location = tc.device('cpu')))
    model.eval()


    print('training finished, starting LRP')

    calc_all_paths(model, train_data, train_onehot, sample_id, train_names[sample_id], featurenames=featurenames, data_type=data_type,
                   result_path = result_path, device=tc.device('cpu'))


