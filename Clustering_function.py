import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
from astropy.table import Table


def find_center(Data, Flux, atom_max, min_samples, eps):
    """
    find the centers for each data set using DBSCAN.
    Data should be a n*2 array.
    """
    db = DBSCAN(eps, min_samples).fit(Data)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    atom = len(set(labels)) - (1 if -1 in labels else 0) # number of atoms
    print ('We find', atom, 'atoms')
    if atom > atom_max:
        print ('You need to increase the variable: min_samples')
        print (atom, 'is even bigger than the possible biggest atom number', atom_max)
    unique_labels = set(labels)
    noise_label = (labels == -1)
    noise = Data[noise_label & ~core_samples_mask]
    if len(noise) == 0:
        print ('There is no outlier.')
    else:
        print ('There are', len(noise), 'outliers. Please have a look at them.')
    xy = []
    flux = []
    centers = np.zeros((atom,2))
    for k in range(atom):
        class_member_mask = (labels == k)
        xy += [Data[class_member_mask & core_samples_mask]]
        flux += [Flux[class_member_mask & core_samples_mask]]
        centers[k] = np.mean(Data[class_member_mask & core_samples_mask], axis = 0)
    print ('The', atom, 'cluster centers are \n', centers)
    return atom, xy, flux, noise, labels, centers

def distance(point1, point2):
    """
    Euclidean distance from point1 and point2
    >>>distance([0.1, 0.1], [0.2, 0.1])
    >>>0.1
    >>>distance([0.2, 0.1], [0.2, 0.1])
    >>>0.0
    """
    point1[0] = np.floor(point1[0]) + 0.5
    point1[1] = np.floor(point1[1]) + 0.5
    point2[0] = np.floor(point2[0]) + 0.5
    point2[1] = np.floor(point2[1]) + 0.5
    return np.abs((point1[0] - point2[0]) + complex(0,1)*(point1[1]- point2[1]))

def angle_between(p1, p2):
    p1x = p1[0]
    p1y = p1[1]
    p2x = p2[0]
    p2y = p2[1]
    ang = np.arctan((p1y - p2y)/(p1x - p2x))
    return ang/np.pi*180

def angle_among(p1, p2, p3):
    p1x = p1[0]
    p1y = p1[1]
    p2x = p2[0]
    p2y = p2[1]
    p3x = p3[0]
    p3y = p3[1]
    ang1 = np.arctan((p1y - p2y)/(p1x - p2x))
    ang2 = np.arctan((p1y - p3y)/(p1x - p3x))
    ang3 = np.arctan((p2y - p3y)/(p2x - p3x))
    return (ang1 + ang2 + ang3)/3./np.pi*180

def brightest(flux,atom):
    avg_flux = []
    index = []
    for i in range(0,atom):
        avg_flux += [np.mean(flux[i])]
    reorder_flux = sorted(avg_flux, reverse=True)
    for i in range(0,atom):
        index += [avg_flux.index(reorder_flux[i])]
    return avg_flux, index

cluster_num = []
d = []
ang = []
filename = str(input('Enter the result file name:\n'))
min_samples = int(input('Enter the minimum group number (you can try to start with 10):\n'))
eps = float(input('Enter the searching distance for each group (you can try to start with 2):\n'))

## 1. Start the initial preprocess the input table###
print ('Processing ' + filename + '...')
infile = open(filename, 'r')
infile.readline()  # skip the first line
data = [line.split() for line in infile]
infile.close()
Data_original = [float(j) for j in data[0]]
for i in range(1,len(data)):
    Data_original = np.vstack([Data_original,[float(j) for j in data[i]]])
Data_c = Data_original
Data_vxy = Data_c[:,1:3]
Flux = Data_c[:,3]
n = 0
col_num = []
for i in range(1,len(Data_c)):
    if Data_c[i,0] == Data_c[i-1,0]:
        n += 1
    else:
        col_num += [n+1]
        n = 0
atom_max = int(np.max(col_num))

## 2. First Clustering process ##
atom, xy, flux, noise, labels, centers = find_center(Data_vxy, Flux, atom_max, min_samples, eps)
#col = ['b', 'g', 'r', 'c', 'm'] * atom
if atom > 0:
    fig1 = plt.figure(1)
    for i in range(atom):
        l1 = plt.scatter(xy[i][:,0], xy[i][:,1], c = 'k', s = 10)
        l2 = plt.scatter(centers[i][0], centers[i][1], c = 'r', s = 100, marker = '*')
    plt.grid(True)
    plt.legend((l1, l2), ('Clusters', 'Centers'), ncol=2, fontsize=8)
    plt.title('Clustering of the MCMC results for file ' + filename)
    plt.xlabel('Image plane x axis')
    plt.ylabel('Image plane y axis')
    plt.gca().set_aspect('equal', adjustable='box')
    figname = 'cluster_results_' + filename + '.png'
    fig1.savefig(figname)
    if atom == 2:
        d += [distance(centers[0], centers[1])]
        ang += [angle_between(centers[0], centers[1])]
    elif atom == 1:
        d += [0]
        ang += [0]
    else:
        d += [-100]
        ang += [-100]
    cluster_num += [len(centers)]

## 3. Second Clustering process ##
xy = np.asarray(xy)
Data_vxy = xy[0]
flux = np.asarray(flux)
Flux = flux[0]
for i in range(1,atom):
    Data_vxy = np.append(Data_vxy, xy[i], axis=0)
    Flux = np.append(Flux, flux[i], axis=0)

min_samples1 = int(input('Enter the minimum group number:\n'))
eps1 = float(input('Enter the searching distance for each group:\n'))
atom1, xy1, flux1, noise1, labels1, centers1 = find_center(Data_vxy, Flux, atom_max, min_samples1, eps1)

if atom1 > 0:
    if (atom1 < atom):
        file_num += [file_i]
    fig2 = plt.figure(2)
    for i in range(atom1):
        l3 = plt.scatter(noise1[:,0], noise1[:,1], c = 'r', s = 30)
        l4 = plt.scatter(xy1[i][:,0], xy1[i][:,1], s = 50)
        l5 = plt.scatter(centers1[:,0], centers1[:,1], c = 'k', s = 50, marker = '*')
    plt.grid(True)
    plt.legend((l3, l4, l5), ('Clusters', 'Centers', 'Noise'), ncol=3, fontsize=8)
    #plt.title('Clustering of the MCMC result ' + str(file_i) + ' for Basc')
    plt.xlabel('Image plane x axis')
    plt.ylabel('Image plane y axis')
    plt.gca().set_aspect('equal', adjustable='box')
    figname = 'Demo_results_' + filename + '.png'
    fig2.savefig(figname)
    fig3 = plt.figure(3)
    for i in range(atom1):
        l6 = plt.scatter(xy1[i][:,0], xy1[i][:,1], c = 'r', s = 10)
        l7 = plt.scatter(centers1[i][0], centers1[i][1], c = 'k', s = 100, marker = '*')
    plt.grid(True)
    plt.legend((l4, l5), ('Clusters', 'Centers'), ncol=2, fontsize=8)
    #plt.title('Clustering of the MCMC results ' + str(file_i) + ' for Basc')
    plt.xlabel('Image plane x axis')
    plt.ylabel('Image plane y axis')
    plt.gca().set_aspect('equal', adjustable='box')
    figname = 'cluster_results_' + filename + '_1.png'
    fig3.savefig(figname)

if atom1 == 2:
    d += [distance(centers1[0], centers1[1])]
    ang += [angle_between(centers1[0], centers1[1])]
    cluster_num += [len(centers1)]
elif atom == 2 and atom1 < 2:
    d += [distance(centers[0], centers[1])]
    ang += [angle_between(centers[0], centers[1])]
    cluster_num += [len(centers)]
elif atom == 1 and atom1 <= 1:
    d += [0]
    ang += [0]
    cluster_num += [len(centers1)]
else:
    d += [-100]
    ang += [-100]
    cluster_num += [len(centers1)]

plt.close('all')
#d = np.asarray(d).T
#ang = np.asarray(ang).T
#cluster_num = np.asarray(cluster_num).T
#print (no.shape, d.shape, ang.shape, cluster_num.shape)
#result = np.vstack((no, cluster_num, d, ang)).T
#np.savetxt('pipeline_result_3_2_retry.csv', result, delimiter=',')
#print (file_num)
