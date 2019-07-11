import numpy as np
import math
import torch 
import torchvision.datasets as dsets
import torchvision.transforms as transforms
import torchvision
import torch.nn as nn
from torch.autograd import Variable

import sys
from time import time

class Autoencoder(nn.Module):
    def __init__(self, in_dim=15, h_dim=4, c_dim=1):
        super(Autoencoder, self).__init__()

        self.encoder = nn.Sequential(
            nn.Linear(in_dim, h_dim),
            nn.ReLU()
            )

        self.encoder2 = nn.Sequential(
            nn.Linear(h_dim, 1),
            nn.ReLU()
            )

#        self.decoder = nn.Sequential(
#            nn.Linear(h_dim+c_dim, in_dim),
#            nn.Sigmoid()
#            )
        self.decoder = nn.Linear(1+c_dim, in_dim)

    def forward(self, input, context):
        latent = self.encoder(input)
        combined = torch.cat((self.encoder2(latent), context), dim=1)
        out = self.decoder(combined)
        return out


def to_var(x):
    if torch.cuda.is_available():
        x = x.cuda()
    return Variable(x)

def main(data_loader, args):
    ae = Autoencoder(in_dim=args.in_dim, h_dim=args.hidden_size)

    if torch.cuda.is_available():
        print("Use CUDA")
        ae.cuda()

    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(ae.parameters(), lr=0.001)

    for epoch in range(args.num_epochs):
        t0 = time()
        for i, (X, context) in enumerate(data_loader):
            out = ae(X, context)
            loss = criterion(out, X)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if (i+1) % 100 == 0:
                print ('Epoch [%d/%d], Iter [%d/%d] Loss: %.4f Time: %.2fs' 
                    %(epoch+1, args.num_epochs, i+1, len(data_loader.dataset)//args.batch_size, loss.data.item(), time()-t0))

    return ae
    
        # save the reconstructed images
        #reconst_images = ae.encoder(fixed_x).data.numpy()
        #print(reconst_images)
        #reconst_images = reconst_images.view(reconst_images.size(0), 1, 28, 28)
        #torchvision.utils.save_image(reconst_images.data.cpu(), './data/reconst_images_%d.png' % (epoch+1))

class Options():
    def __init__(self):
        self.batch_size = 10
        self.num_epochs = 30
        self.in_dim = 14
        self.hidden_size = 4


def readContext(fl):
    with open(fl, 'r') as f:
        return [[math.log(int(l.strip().split('\t')[1]))] for l in f]

args = Options()

data = torch.tensor(np.loadtxt(sys.argv[1]), dtype=torch.float32)
context = torch.tensor(readContext(sys.argv[2]), dtype=torch.float32)
dataset = torch.utils.data.TensorDataset(data, context)
data_loader = torch.utils.data.DataLoader(dataset=dataset,
    batch_size=args.batch_size, shuffle=True)

ae = main(data_loader, args)
np.savetxt(sys.argv[3], ae.encoder(data).data.numpy())