import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

def find_center(Data, atom_max, min_samples = 10, eps = 3):
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
    centers = np.zeros((atom,2))
    for k in range(atom):
        class_member_mask = (labels == k)
        xy += [Data[class_member_mask & core_samples_mask]]
        centers[k] = np.mean(Data[class_member_mask & core_samples_mask], axis = 0)
    print ('The', atom, 'cluster centers are \n', centers)
    return atom, xy, noise, labels, centers

def distance(point1, point2):
    """
    Euclidean distance from point1 and point2
    >>>distance([0.1, 0.1], [0.2, 0.1])
    >>>0.1
    >>>distance([0.2, 0.1], [0.2, 0.1])
    >>>0.0
    """
    return np.abs((point1[0] - point2[0]) + complex(0,1)*(point1[1]- point2[1]))

def angle_between(p1, p2):
    p1x = p1[0]
    p1y = p1[1]
    p2x = p2[0]
    p2y = p2[1]
    ang = np.arctan((p1y - p2y)/(p1x - p2x))
    return ang/np.pi*180

cluster_num = []
d = []
ang = []

for file_i in range(100):
    filename = 'results_' + str(file_i) + '.txt'
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
    n = 0
    col_num = []
    for i in range(1,len(Data_c)):
        if Data_c[i,0] == Data_c[i-1,0]:
            n += 1
        else:
            col_num += [n+1]
            n = 0
    atom_max = int(np.max(col_num))
    atom, xy, noise, labels, centers = find_center(Data_vxy, atom_max)
    plt.figure(file_i)
    for i in range(atom):
        plt.scatter(xy[i][:,0], xy[i][:,1], c = 'r', s = 10)
        plt.scatter(centers[i][0], centers[i][1], c = 'k', s = 100, marker = '*')
    figname = 'cluster_results_' + str(file_i) + '.eps'
    plt.savefig(figname)
    cluster_num += [len(centers)]
    if len(centers) == 2:
        d += [distance(centers[0], centers[1])]
        ang += [angle_between(centers[0], centers[1])]
    elif len(centers) == 1:
        d += [0]
        ang += [0]
    else:
        d += [-1]
        ang += [-1]
no = np.asarray([i for i in range(100)])
d = np.asarray(d).T
ang = np.asarray(ang).T
cluster_num = np.asarray(cluster_num).T
result = np.vstack((no, cluster_num, d, ang)).T
np.savetxt('pipeline_result.csv', result, delimiter=',')
