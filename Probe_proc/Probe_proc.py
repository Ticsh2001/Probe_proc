import numpy as np
import matplotlib.pyplot as plt
import os


def plot_hist(ax, bins, edges, type, col, name):
    if type == 'bars':
        ax.bar(edges[0],bins[0],edges[0+1]-edges[0],align='edge', edgecolor  = col, label = name, alpha = 1.0, fill = False,linewidth  = 2)
        for i in np.arange(1,len(bins)):
            ax.bar(edges[i],bins[i],edges[i+1]-edges[i],align='edge', edgecolor  = col, alpha = 1.0, fill = False,linewidth  = 2)
    if type == 'plot':
        ax.plot(edges[:-1],bins,label = name, color = col)
 
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
    res = np.zeros((len(n_bins),1))
    for i in np.arange(0,len(n_bins)):
        res[i] = n_edges[i]**3*n_bins[i]
    s = np.sum(res)
    return res/s
    #for i in np.arange(0,len(n_bins)):

def load_ipi_probe(path):
    files = (f for f in os.listdir(path) if f.endswith('.txt'))
    data = []
    for f in files:
        temp = np.loadtxt(path+'\\'+f)
        for z in temp:
            if z>0.0e-6:
                data.append(z*2.0*1000000)
                
    [n_bins, n_edges] = np.histogram(data, bins='doane')
    n_cum_data = create_freq_sum(n_bins, n_edges, False)
    return n_bins, n_edges,n_cum_data

def load_lmz_probe(path):
    data = []
    with open(path,'r') as f:
        for line in f:
            str = line.replace(',','.').replace(' ','').replace('\n','')
            if str.count('.') > 0:
                data.append(float(str))
        f.close()
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
   
def load_ipi_conditions(conditions, names, axs):
    j = 0
    colors = ['r','g','b','k']
    for cond in conditions:
        ipi_n_bins, ipi_n_edges, ipi_cum = load_ipi_probe(cond)
        ipi_m_bins = create_mass_distrib(ipi_n_bins, ipi_n_edges)
        print(np.sum(ipi_m_bins))
        plot_hist(axs[0], ipi_n_bins/np.sum(ipi_n_bins),ipi_n_edges/ipi_n_edges[-1],'bars',colors[j],names[j])
        axs[1].plot(ipi_cum[:,0],ipi_cum[:,1],'-'+colors[j]+'o', label = names[j])
        j = j + 1



f, axs = plt.subplots(1,2)

n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\08.02.2020\f3f4\amplitude.txt')
m_bins = create_mass_distrib(n_bins, n_edges)
plot_hist(axs[0], n_bins/np.sum(n_bins),n_edges[:]/n_edges[-1],'plot','r','lmz f3f4')
n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\08.02.2020\f5f6\amplitude.txt')
m_bins = create_mass_distrib(n_bins, n_edges)
plot_hist(axs[0], n_bins/np.sum(n_bins),n_edges[:]/n_edges[-1],'plot','g','lmz f5f6')
n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\08.02.2020\f9f12\amplitude.txt')
m_bins = create_mass_distrib(n_bins, n_edges)
plot_hist(axs[0], n_bins/np.sum(n_bins),n_edges[:]/n_edges[-1],'plot','b','lmz f9f12')

ipi_conditions = []
#ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f3f4_2_0')
#ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f3f4_new_filter')
#ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f9f12_2_0')

ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f3f4')
ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f5f6')
ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f9f12')

ipi_names = []
#ipi_names.append(r'f3f4 ipi_2_0')
ipi_names.append(r'f3f4 ipi')
ipi_names.append(r'f3f4 ipi_new')
ipi_names.append(r'f9f12 ipi')
load_ipi_conditions(ipi_conditions, ipi_names, axs)



#n_bins, n_edges, n_cum_data = load_footprint_probe_data(r'F:\LMZ\test_foot.txt', 1.0)
#m_bins = create_mass_distrib(n_bins, n_edges)
#lmz_bins, lmz_edges, lmz_cum_data = load_LMZ_probe(r'F:\LMZ\probe_1.txt')
#plot_hist(axs[0], m_bins,n_edges,'bars','r')
#plot_hist(axs[0], lmz_bins,lmz_edges,'bars','r')
#ipi_n_bins, ipi_n_edges, ipi_cum = load_ipi_probe(r'T:\POLIS\Зонд влажности\2020.02.08\Export\f3f4')
#ipi_m_bins = create_mass_distrib(ipi_n_bins, ipi_n_edges)
#plot_hist(axs[0], ipi_m_bins,ipi_n_edges,'bars','r')
#foot_cum_x = (n_cum_data[:,0])/(n_cum_data[-1,0])
#lmz_cum_x = (lmz_cum_data[:,0])/(lmz_cum_data[-1,0])
#axs[1].plot(foot_cum_x,n_cum_data[:,1],'-ro')
#axs[1].plot(lmz_cum_x,lmz_cum_data[:,1],'-go')
axs[0].legend()
axs[1].legend()
plt.show()