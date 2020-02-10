import numpy as np
import matplotlib.pyplot as plt
import os


def plot_hist(ax, bins, edges, type, col):
    if type == 'bars':
        for i in np.arange(0,len(bins)):
            ax.bar(edges[i],bins[i],edges[i+1]-edges[i],align='edge', color = col)
    #if type == 'plot':
 
def create_freq_sum(bins, edges, is_normolized):
    res = np.zeros((len(bins),2))
    for i in np.arange(0,len(bins)):
        if i == 0:
            res[i,1] = bins[i]
        else:
            res[i,1] = (res[i-1,1] + bins[i])
        res[i,0] = (edges[i+1]+edges[i])/2;
    if is_normolized==False:
        s = np.sum(bins)
        res[:,1] = res[:,1]/s
    return res

def create_mass_distrib(n_bins, n_edges):
    res = n_edges[0:-1]**3*n_bins
    s = np.sum(res)
    return res/s
    #for i in np.arange(0,len(n_bins)):

def load_ipi_probe(path):
    files = (f for f in os.listdir(path) if f.endswith('.txt'))
    data = []
    for f in files:
        temp = np.loadtxt(path+f)
        for z in temp:
            if z>0.0:
                data.append(z*2.0)
    [n_bins, n_edges] = np.histogram(data, bins='doane')
    n_cum_data = create_freq_sum(n_bins, n_edges, False)
    return n_bins, n_edges,n_cum_data



    
def load_LMZ_probe(path):
    data = np.loadtxt(path)
    j = 0
    while data[j,0]<2.5:
        j = j+1
    n_bins = np.zeros((len(data)-j,1))
    n_edges = np.zeros((len(data)+1-j,1))
    k = 0
    for i in np.arange(j,len(data)):
        n_bins[k] = data[i,1]/100
        if k == 0:
            n_edges[k] = data[i,0]-(data[i+1,0]-data[i,0])/2
        n_edges[k+1] = data[i,0]+(data[i,0]-n_edges[k])
        k = k+1
    s = np.sum(n_bins)
    n_bins = n_bins/s
    n_cum_data = create_freq_sum(n_bins, n_edges, True)
    return n_bins, n_edges,n_cum_data



def load_footprint_probe_data(path, scale_coeff):
    data = np.loadtxt(path)
    [n_bins, n_edges] = np.histogram(data, bins='doane')
    n_edges = n_edges*scale_coeff
    n_cum_data = create_freq_sum(n_bins, n_edges, False)
    return n_bins, n_edges,n_cum_data
   


f, axs = plt.subplots(1,2)
n_bins, n_edges, n_cum_data = load_footprint_probe_data(r'F:\LMZ\test_foot.txt', 1.0)
m_bins = create_mass_distrib(n_bins, n_edges)
lmz_bins, lmz_edges, lmz_cum_data = load_LMZ_probe(r'F:\LMZ\probe_1.txt')
plot_hist(axs[0], m_bins,n_edges,'bars','r')
load_ipi_probe(path)
#plot_hist(axs[0], lmz_bins,lmz_edges,'bars','r')

foot_cum_x = (n_cum_data[:,0])/(n_cum_data[-1,0])
lmz_cum_x = (lmz_cum_data[:,0])/(lmz_cum_data[-1,0])
axs[1].plot(foot_cum_x,n_cum_data[:,1],'-ro')
axs[1].plot(lmz_cum_x,lmz_cum_data[:,1],'-go')
plt.show()