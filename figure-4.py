import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import warnings
warnings.filterwarnings("ignore")

def histogram_3d_PQD(data,colour,filename):
    
    raw_data = []
    with open(data) as f:
                raw_data.extend([tuple(map(float, i.split(' '))) for i in f])

    ps = [i[0] for i in raw_data]
    densities = [np.array(i[1:-1]) for i in raw_data]

    ps = ps[0:-1:10]
    densities = densities[0:-1:10]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect((12,12,4))
    ax.view_init(azim=-30)

    counter = 0
    binns = np.linspace(0,0.5,51)

    for c, z in zip([colour]*len(ps), ps):

        hist, bins = np.histogram(densities[counter], bins=binns, density = True)
        xs = (bins[:-1] + bins[1:])/2

        ax.bar(xs, hist, zs=z, width=0.02, zdir='y', color=c, ec=c, alpha=0.8,linewidth=0.0)
        counter += 1

    ax.set_xlabel('density ($\\rho$)', fontsize=20,labelpad=25)
    ax.set_ylabel('baseline birth probability ($p$)', fontsize=20,labelpad=25)
    ax.set_zlabel('PDF ($P(\\rho))$', fontsize=20,labelpad=15)
    ax.set_xticks([0.0,0.25,0.5])
    ax.set_yticks([0.4,0.45,0.5,0.55])
    ax.set_zticks([0,50,100])
    ax.xaxis.set_tick_params()
    ax.tick_params(axis='both', which='major', labelsize=15)

    #plt.tight_layout()
    plt.savefig(filename,dpi=300,bbox_inches='tight')

def find_bistable_region(data):
    
    raw_data = []
    with open(data) as f:
                raw_data.extend([tuple(map(float, i.split(' '))) for i in f])
    
    ps = [i[0] for i in raw_data]
    densities = [np.array(i[1:-1]) for i in raw_data]
    
    upper_arm = np.array([np.mean([j for j in i if j!=0]) for i in densities])
    upper_ps = np.array([ps[i] for i in range(len(ps)) if not np.isnan(upper_arm[i])])
    upper_arm = np.array([i for i in upper_arm if not np.isnan(i)])

    lower_arm = np.array([np.mean([j for j in i if j==0]) for i in densities])
    lower_ps = np.array([ps[i] for i in range(len(ps)) if not np.isnan(lower_arm[i])])
    lower_arm = np.array([i for i in lower_arm if not np.isnan(i)])
    
    return upper_ps[0], lower_ps[-1]

def draw_bifurcation_diagram(data,grid_size,paint):
    
    raw_data = []
    with open(data) as f:
                raw_data.extend([tuple(map(float, i.split(' '))) for i in f])
    
    ps = [i[0] for i in raw_data]
    densities = [np.array(i[1:-1]) for i in raw_data]
    
    upper_arm = np.array([np.mean([j for j in i if j!=0]) for i in densities])
    upper_ps = np.array([ps[i] for i in range(len(ps)) if not np.isnan(upper_arm[i])])
    upper_arm = np.array([i for i in upper_arm if not np.isnan(i)])
    
    upper_arm_errors = np.array([np.std([j for j in i if j!=0]) for i in densities])
    upper_arm_errors = np.array([i for i in upper_arm_errors if not np.isnan(i)])

    lower_arm = np.array([np.mean([j for j in i if j==0]) for i in densities])
    lower_ps = np.array([ps[i] for i in range(len(ps)) if not np.isnan(lower_arm[i])])
    lower_arm = np.array([i for i in lower_arm if not np.isnan(i)])
    
    plt.fill_between(upper_ps, upper_arm-upper_arm_errors, upper_arm+upper_arm_errors,alpha=0.2,facecolor=paint,linewidth=1, linestyle='dashdot')

    plt.plot(upper_ps,upper_arm,linewidth=4,color=paint,label="$L = $"+grid_size)
    plt.plot(lower_ps,lower_arm,linewidth=4,color=paint)

    plt.axvline(x=upper_ps[0],color='black',linestyle='--',linewidth=3,ymin=0.0, ymax=(upper_arm[0]/0.5))
    plt.axvline(x=lower_ps[-1],color='black',linestyle='--',linewidth=3,ymin=0.0, ymax=(upper_arm[np.where(upper_ps==lower_ps[-1])]/0.5))
    plt.xlim(0.4,0.6)
    plt.ylim(0.0,0.5)
    #plt.xlabel('baseline birth probability $(p)$',fontsize=20)
    #plt.ylabel('steady-state density $(\\rho^{*})$',fontsize=20)
    plt.xticks(fontsize=12)
    #plt.yticks(fontsize=12)
    plt.locator_params(axis='x', nbins=5)
    plt.yticks([0.25,0.5], fontsize=12)
    plt.legend(loc='upper left',fontsize=15)
    #plt.savefig('bif.png',dpi=300)

    histogram_3d_PQD('16.txt','blue','16_3d.png')
    histogram_3d_PQD('128.txt','red','128_3d.png')

    data= []
with open("128.txt") as f:
    data.extend([tuple(map(float, i.split(' '))) for i in f])

ps_128 = [i[0] for i in data]
densities_128 = [np.array(i[1:-1]) for i in data]

entire_128 = np.array([np.mean(i) for i in densities_128])
entire_128_errors = np.array([np.std(i) for i in densities_128])

fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(15,4))
plt.subplot(1,3,1)
draw_bifurcation_diagram('16.txt','16','blue')
plt.subplot(1,3,2)
draw_bifurcation_diagram('32.txt','32','orange')
plt.subplot(1,3,3)
plt.fill_between(ps_128, entire_128-entire_128_errors, entire_128+entire_128_errors,alpha=0.2,facecolor='red',linewidth=1, linestyle='dashdot')
plt.plot(ps_128,entire_128,color='red',linewidth=4,label='$L = 128$')
plt.xlim(0.4,0.6)
plt.ylim(-0.01,0.5)
plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
plt.locator_params(axis='x', nbins=5)
plt.yticks([0.25,0.5], fontsize=12)
plt.legend(loc='upper left',fontsize=15)
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("baseline birth probability $(p)$",fontsize=18,labelpad = 10)
plt.ylabel("steady-state density $(\\rho^{*})$",fontsize=18,labelpad=15)
plt.tight_layout()
plt.savefig('bif.png',dpi=300)
