import torch as tc
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import copy

from torch.utils.data import Dataset, DataLoader
from torch.optim.lr_scheduler import ExponentialLR
from itertools import permutations
import pandas as pd

from dataloading_simple import Dataset_train, Dataset_LRP
import os


class LogCoshLoss(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, y_t, y_prime_t):
        ey_t = y_t - y_prime_t
        return tc.mean(tc.log(tc.cosh(ey_t + 1e-12)))


class LRP_Linear(nn.Module):
    def __init__(self, inp, outp, gamma=0.01, eps=1e-5):
        super(LRP_Linear, self).__init__()
        self.A_dict = {}
        self.linear = nn.Linear(inp, outp)
        nn.init.xavier_uniform_(self.linear.weight, gain=nn.init.calculate_gain('relu'))
        self.gamma = tc.tensor(gamma)
        self.eps = tc.tensor(eps)
        self.rho = None
        self.iteration = None

    def forward(self, x):

        if not self.training:
            self.A_dict[self.iteration] = x.clone()
        return self.linear(x)

    def relprop(self, R):
        device = next(self.parameters()).device

        A = self.A_dict[self.iteration].clone()
        A, self.eps = A.to(device), self.eps.to(device)

        Ap = A.clamp(min=0).detach().data.requires_grad_(True)
        Am = A.clamp(max=0).detach().data.requires_grad_(True)


        zpp = self.newlayer(1).forward(Ap)  
        zmm = self.newlayer(-1, no_bias=True).forward(Am) 

        zmp = self.newlayer(1, no_bias=True).forward(Am) 
        zpm = self.newlayer(-1).forward(Ap) 

        with tc.no_grad():
            Y = self.forward(A).data

        sp = ((Y > 0).float() * R / (zpp + zmm + self.eps * ((zpp + zmm == 0).float() + tc.sign(zpp + zmm)))).data # new version
        sm = ((Y < 0).float() * R / (zmp + zpm + self.eps * ((zmp + zpm == 0).float() + tc.sign(zmp + zpm)))).data

        (zpp * sp).sum().backward()
        cpp = Ap.grad
        Ap.grad = None
        Ap.requires_grad_(True)

        (zpm * sm).sum().backward()
        cpm = Ap.grad
        Ap.grad = None
        Ap.requires_grad_(True)

        (zmp * sp).sum().backward()
        cmp = Am.grad

        Am.grad = None
        Am.requires_grad_(True)

        (zmm * sm).sum().backward()
        cmm = Am.grad
        Am.grad = None
        Am.requires_grad_(True)


        R_1 = (Ap * cpp).data
        R_2 = (Ap * cpm).data
        R_3 = (Am * cmp).data
        R_4 = (Am * cmm).data


        return R_1 + R_2 + R_3 + R_4

    def newlayer(self, sign, no_bias=False):

        if sign == 1:
            rho = lambda p: p + self.gamma * p.clamp(min=0) # Replace 1e-9 by zero
        else:
            rho = lambda p: p + self.gamma * p.clamp(max=0) # same here

        layer_new = copy.deepcopy(self.linear)

        try:
            layer_new.weight = nn.Parameter(rho(self.linear.weight))
        except AttributeError:
            pass

        try:
            layer_new.bias = nn.Parameter(self.linear.bias * 0 if no_bias else rho(self.linear.bias))
        except AttributeError:
            pass

        return layer_new


class LRP_ReLU(nn.Module):
    def __init__(self):
        super(LRP_ReLU, self).__init__()
        self.relu = nn.ReLU()

    def forward(self, x):
        return self.relu(x)

    def relprop(self, R):
        return R



class Model(nn.Module):
    def __init__(self, inp, outp, hidden, hidden_depth):
        super(Model, self).__init__()
        self.layers = nn.Sequential(LRP_Linear(inp, hidden), LRP_ReLU())
        for i in range(hidden_depth):
            self.layers.add_module('LRP_Linear' + str(i + 1), LRP_Linear(hidden, hidden))
            self.layers.add_module('LRP_ReLU' + str(i + 1), LRP_ReLU())

        self.layers.add_module('LRP_Linear_last', LRP_Linear(hidden, outp))

    def forward(self, x):
        return self.layers.forward(x)

    def relprop(self, R):
        assert not self.training, 'relprop does not work during training time'
        for module in self.layers[::-1]:
            R = module.relprop(R)
        return R


def train(neuralnet, train_data, test_data, train_onehot,test_onehot, epochs, lr, noise_level, weight_decay = 0.0, batch_size=50, device=tc.device('cpu')):
    nsamples, nfeatures = train_data.shape
    optimizer = tc.optim.SGD(neuralnet.parameters(), lr=lr, momentum=0.9, weight_decay=weight_decay) 
    scheduler = ExponentialLR(optimizer, gamma = 0.99)

    criterion = LogCoshLoss() #nn.MSELoss()
    testlosses, epoch_list, network_list = [], [], []

    neuralnet.train().to(device)

    for epoch in range(epochs):
        if epoch<20:
            optimizer.param_groups[0]['lr']=0.001*(epoch+1)

        trainset = Dataset_train(train_data, train_onehot, noise_level=0.0)
        trainloader = DataLoader(trainset, batch_size=batch_size, shuffle=True)

        
        with tc.no_grad():
            pass#for param in neuralnet.parameters():
            #    param.add_(tc.randn(param.size(), device=device)*noise_level)
      
        for masked_data, mask, full_data in trainloader:
            masked_data = masked_data.to(device)
            mask = mask.to(device)
            full_data = full_data.to(device)
            #print(masked_data, full_data)

            optimizer.zero_grad()
            pred = neuralnet(masked_data)
            loss = criterion(pred[mask==0], full_data[mask==0]) 
            loss.backward()
            optimizer.step()
        scheduler.step()
            

        if epoch%10==0:
            print(optimizer.param_groups[0]['lr'])
            neuralnet.eval()
            testset = Dataset_train(test_data,test_onehot, noise_level=0)
            traintestset = Dataset_train(train_data,train_onehot, noise_level=0)
            testloader = DataLoader(testset, batch_size=test_data.shape[0], shuffle=False)
            traintestloader = DataLoader(traintestset, batch_size=test_data.shape[0], shuffle=False)

            for masked_data, mask, full_data in testloader:
                masked_data = masked_data.to(device)
                mask = mask.to(device)
                full_data = full_data.to(device)
                with tc.no_grad():
                    pred = neuralnet(masked_data)
                testloss = criterion(pred[mask==0], full_data[mask==0])
                testlosses.append(testloss)
                epoch_list.append(epoch)
                network_list.append(neuralnet.state_dict())
                break

            for masked_data, mask, full_data in traintestloader:
                masked_data = masked_data.to(device)
                mask = mask.to(device)
                full_data = full_data.to(device)
                with tc.no_grad():
                    pred = neuralnet(masked_data)
                traintestloss = criterion(pred[mask==0], full_data[mask==0])
                print(epoch, 'trainloss:', traintestloss, 'testloss:', testloss)
                break
    return tc.tensor(testlosses), epoch_list, network_list




def compute_LRP(neuralnet, test_set, test_onehot, target_id, sample_id, batch_size, device):
    criterion = nn.MSELoss()
    testset = Dataset_LRP(test_set,test_onehot, target_id, sample_id)
    testloader = DataLoader(testset, batch_size=batch_size, shuffle=True)

    neuralnet.to(device).eval()

    masked_data, mask, full_data = next(iter(testloader))
    masked_data, mask, full_data = masked_data.to(device), mask.to(device), full_data.to(device)
    pred = neuralnet(masked_data)

    error = criterion(pred.detach()[:,target_id], full_data.detach()[:,target_id]).cpu().numpy()
    y = full_data.detach()[:,target_id].cpu().mean().numpy()
    y_pred = pred.detach()[:,target_id].cpu().mean().numpy()

    R = tc.zeros_like(pred)
    R[:,target_id] = pred[:,target_id].clone()
    #R = R.to(device)
    a = neuralnet.relprop(R)
    LRP_sum = (a.sum(dim=0))

    LRP_unexpanded = 0.5 * (LRP_sum[:LRP_sum.shape[0] // 2] + LRP_sum[LRP_sum.shape[0] // 2:])


    mask_sum = mask.sum(dim=0).float()

    LRP_scaled = LRP_unexpanded/mask_sum
    LRP_scaled = tc.where(tc.isnan(LRP_scaled),tc.tensor(0.0).to(device), LRP_scaled)
     
    full_data_sample = full_data[0,:].cpu().detach().numpy().squeeze()
    return LRP_scaled.cpu().numpy(), error, y , y_pred, full_data_sample


def calc_all_paths(neuralnet, test_data,test_onehot, sample_id, sample_name, featurenames, data_type, result_path, device = tc.device('cuda:0')):
    #if not os.path.exists(result_path + data_type + '/'):
    #    os.makedirs(result_path + data_type + '/')
    end_frame = []
    print(sample_id)
    for target in range(test_data.shape[1]):
        LRP_value, error, y, y_pred, full_data_sample = compute_LRP(neuralnet, test_data, test_onehot, target, sample_id, batch_size = 100, device = device)

        frame = pd.DataFrame({'LRP': LRP_value, 'source_gene': featurenames, 'target_gene': featurenames[target] ,'sample_name': sample_name, 'error':error, 'y':y, 'y_pred':y_pred, \
            'inpv': full_data_sample})
       

        end_frame.append(frame)
        end_result_path = result_path + data_type + '/'  + 'LRP_' + str(sample_id) + '_'+ sample_name + '.csv'
        if not os.path.exists(result_path + data_type):
            os.makedirs(result_path + data_type)

    end_frame = pd.concat(end_frame, axis=0)
    end_frame.to_csv(end_result_path)

