import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.neighbors import KernelDensity

def prepare_lmz_file(path):
    data = []
    files = (f for f in os.listdir(path) if f.endswith('.txt'))
    for file_p in files:
        fl = open(path+'\\'+file_p,'r')
        str = fl.read()
        str = str.replace('\n','\t')
        str = str.replace(' ','')
        str = str.replace('\t','\n')
        str = str.replace(',','.')
        fl.close()
        fl = open(path+'\\'+file_p+'1','w')
        fl.write(str)
        fl.close()
    files_new = (f for f in os.listdir(path) if f.endswith('.txt1'))
    for file_p in files_new:
        fl = open(path+'\\'+file_p,'r')
        for line in fl:
            str = line.replace('\n','')
            if str.count('.')>0:
                data.append(float(str))
        fl.close()
    return data



def ipi_adapt_data(bins, edges, data):
    vals = np.zeros((len(bins),1))
    for i in np.arange(0,len(bins)):
        vals[i] = (edges[i+1]+edges[i])/2.0
    kde = KernelDensity(bandwidth = 0.6,kernel = 'exponential')
    kde.fit(np.asarray(data)[:,None])
    new_bins = np.exp(kde.score_samples(np.asarray(vals[:,0])[:,None]))
    new_bins = new_bins*np.sum(bins)
    return new_bins



def analyze_extremum(bins, edges):
    i = np.argmax(bins)
    val = (edges[i]+edges[i+1])/2
    return i , val

def plot_hist(ax, bins, edges, type, col, name):
    if type == 'bars':
        ax.bar(edges[0],bins[0],edges[0+1]-edges[0],align='edge', edgecolor  = col, label = name, alpha = 1.0, fill = False,linewidth  = 2)
        for i in np.arange(1,len(bins)):
            ax.bar(edges[i],bins[i],edges[i+1]-edges[i],align='edge', edgecolor  = col, alpha = 1.0, fill = False,linewidth  = 2)
    if type == 'plot':
        ax.plot(edges[:-1],bins,label = name, color = 'k', linestyle='--')
 
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
    n = 0
    res = np.zeros((len(n_bins),1))
    for i in np.arange(0,len(n_bins)):
        res[i] = (n_edges[i]**3)*n_bins[i]
        n = n + n_bins[i]
    s = np.sum(res)
    d_m_av = (s/n)**(1/3)
    return res/s, d_m_av
    #for i in np.arange(0,len(n_bins)):

def load_ipi_probe(path,xmin, ymin,xmax,ymax):
    files = (f for f in os.listdir(path) if f.endswith('.txt'))
    data = []
    for f in files:
        temp = np.loadtxt(path+'\\'+f)
        for i in np.arange(0,len(temp)):
            if temp[i,0]>xmin and temp[i,0]<xmax and temp[i,1]>ymin and temp[i,1]<ymax and temp[i,2]>0.0:
                data.append(temp[i,2]*2.0*1000000.0)                
    [n_bins, n_edges] = np.histogram(data, bins=28)
    n_bins = ipi_adapt_data(n_bins,n_edges, data)
    #kde = KernelDensity(bandwidth = 1,kernel = 'exponential')
    #kde.fit(np.asarray(data)[:,None])
    #n_bins1 = np.exp(kde.score_samples(np.asarray(n_edges)[:,None]))
    #plt.plot(n_edges[0:-1],n_bins/np.sum(n_bins),'r')
    #plt.plot(n_edges,n_bins1/np.sum(n_bins1),'g')
    #plt.show()
    #print(np.sum(n_bins1))

    n_cum_data = create_freq_sum(n_bins, n_edges, False)
    return n_bins, n_edges,n_cum_data

def load_lmz_probe(path, min, max):
    data = []
    data = prepare_lmz_file(path)

   

    #with open(path,'r') as f:
    #    for line in f:
    #        str = line.replace(',','.').replace(' ','').replace('\n','')
    #        if str.count('.') > 0:
    #            data.append(float(str))
    #    f.close()

    [n_bins, n_edges] = np.histogram(data, bins='doane')
    s = np.sum(n_bins)
    i = 0
    j = len(n_bins)-1
    if min!=-1:
        
        while n_bins[i]/s<min and n_bins[i]<n_bins[i+1]:
            #n_bins.pop(i)
            #n_edges.pop(i)
            i = i + 1
    if max!=-1:
        while n_bins[j]/s < max:
            #n_bins.pop(j)
            #n_edges.pop(j+1)
            j = j - 1
    #n_bins = n_bins[i:j+1]
    #n_edges = n_edges[i:j+2]
    [n_bins, n_edges] = np.histogram(data, bins='doane', range = (n_edges[i],n_edges[j+1]))
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
   
def load_ipi_conditions(conditions, names, axs,xmin,ymin,xmax,ymax):
    j = 0
    colors = ['r','g','b','k']
    mins = []
    maxs = []
    extremum_i = []
    extremum_val = []
    for cond in conditions:
        ipi_n_bins, ipi_n_edges, ipi_cum = load_ipi_probe(cond,xmin,ymin,xmax,ymax)
        ipi_m_bins,ipi_d_m_av = create_mass_distrib(ipi_n_bins, ipi_n_edges)
        i, val = analyze_extremum(ipi_m_bins, ipi_n_edges)
        extremum_i.append(i)
        extremum_val.append(val)
        
        #plot_hist(axs[0], ipi_m_bins,ipi_n_edges/ipi_n_edges[-1],'bars',colors[j],names[j])
        plot_hist(axs[j], ipi_m_bins,ipi_n_edges,'bars',colors[j],names[j])
        #axs[0].plot((ipi_d_m_av/ipi_n_edges[-1],ipi_d_m_av/ipi_n_edges[-1]),(0,np.max(ipi_m_bins)),colors[j])
       # axs[1].plot(ipi_cum[:,0],ipi_cum[:,1],'-'+colors[j]+'o', label = names[j])
        j = j + 1
        mins.append(ipi_n_bins[0]/np.sum(ipi_n_bins))
        maxs.append(ipi_n_bins[-1]/np.sum(ipi_n_bins))
    return mins, maxs,extremum_i,extremum_val



f, axs = plt.subplots(1,3)

ipi_conditions = []

ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.11\2020.02.11_f5f6\Export')
ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.11\2020.02.11_f9f12\Export')
ipi_conditions.append(r'T:\POLIS\Зонд влажности\2020.02.11\2020.02.11_f3f4\Export')

ipi_names = []
ipi_names.append(r'f5f6 ipi')
ipi_names.append(r'f9f12 ipi')
ipi_names.append(r'f3f4 ipi')
#ipi_names.append(r'f9f12 ipi')
xmin = 0.12
xmax = 39.83
ymin = 6.76
ymax = 11.2
mins, maxs,index,extr_v = load_ipi_conditions(ipi_conditions, ipi_names, axs,xmin,ymin,xmax,ymax)


coeff = 13.75
n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\11.02.2020\f5f6',-1,-1)
m_bins,d_m_av = create_mass_distrib(n_bins, n_edges)
#plot_hist(axs[0],m_bins,n_edges[:]/n_edges[-1],'plot','r','lmz f3f4')
i, val = analyze_extremum(m_bins, n_edges)
c = extr_v[0]/val
print(c)
plot_hist(axs[0],m_bins,n_edges*coeff,'plot','r','lmz f5f6')
print(np.sum(m_bins))

n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\11.02.2020\f9f12',-1,-1)
m_bins,d_m_av = create_mass_distrib(n_bins, n_edges)
#plot_hist(axs[1], m_bins,n_edges[:]/n_edges[-1],'plot','g','lmz f5f6')
i, val = analyze_extremum(m_bins, n_edges)
c = extr_v[1]/val
print(c)
plot_hist(axs[1],m_bins,n_edges*coeff,'plot','g','lmz f9f12')

n_bins, n_edges,n_cum_data = load_lmz_probe(r'V:\Стенды КВП\Эксперименты\ЛМЗ тарировка зодна\11.02.2020\f3f4',-1,-1)
m_bins,d_m_av = create_mass_distrib(n_bins, n_edges)
#plot_hist(axs[2], m_bins,n_edges[:]/n_edges[-1],'plot','b','lmz f9f12')
i, val = analyze_extremum(m_bins, n_edges)
c = extr_v[2]/val
print(c)
plot_hist(axs[2],m_bins,n_edges*coeff,'plot','b','lmz f3f4')



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
axs[2].legend()
plt.show()