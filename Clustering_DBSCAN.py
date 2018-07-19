'''

   Author: Haoyang Ye (Zoe)

'''

import numpy as np
from sklearn.cluster import DBSCAN
from astropy.table import Table
import matplotlib.pyplot as plt


def find_center(Data, Flux, atom_max, min_samples, eps):
    """
    Find the cluster centers for given data using DBSCAN package.

    Data: input data, it should be a n*2 array
    atom_max: the MCMC results determine the maximum number of the atoms that can be
        found.
    min_samples: sets the value 'min_samples' in function dbscan, which determines
        the minimum grouping sample number
    eps: the maximum distance between two samples for them to be considered as in
        the same neighborhood
    atom: the number of clusters picked by dbscan
    xy: data arranged according to clustering results
    noise: data that are not selected into any clusters
    labels: flags/labels that determine which data belongs to which cluster
    centers: cluster centre location (x_i,y_i) of each cluster i
    widths: deviation of each clusters
    """
    db = DBSCAN(eps, min_samples).fit(Data)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    atom = len(set(labels)) - (1 if -1 in labels else 0) # number of atoms
    unique_labels = set(labels)
    noise_label = (labels == -1)
    noise = Data[noise_label & ~core_samples_mask]
    xy = []
    flux = []
    centers = np.zeros((atom,2))
    widths = np.zeros((atom,2))
    for k in range(atom):
        class_member_mask = (labels == k)
        xy += [Data[class_member_mask & core_samples_mask]]
        flux += [Flux[class_member_mask & core_samples_mask]]
        centers[k] = np.mean(Data[class_member_mask & core_samples_mask], axis = 0)
        widths[k] = np.std(Data[class_member_mask & core_samples_mask], axis = 0)
    print ('We find ', atom, ' atoms')
    if len(noise) == 0:
        print ('There is no outlier.')
    else:
        print ('There are ', len(noise), 'outliers.')
    print ('The', atom, 'cluster centers are \n', centers)
    print ('The deviation of the ', atom, 'clusters are \n', widths)
    plt.figure()
    plt.scatter(noise.transpose()[0], noise.transpose()[1], s = 20, c = 'k', label = "Outlier")
    for i in range(atom):
        plt.scatter(xy[i][:,0], xy[i][:,1], s = 100, label = "Cluster %d" %i)
    for i in range(atom):
        plt.scatter(centers[i][0], centers[i][1], s = 80, marker = '*', label = "Centre %d" %i)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Image plane x axis')
    plt.ylabel('Image plane y axis')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    return atom, xy, flux, noise, labels, centers, widths

def distance(point1, point2):
    """
    Euclidean distance from cluster centre 1 and cluster centre 2,
    the coordinates of the cluster centres should be located at the
    centre of the cell.
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
    """
    Angle in degrees between cluster centre 1 and cluster centre 2.
    """
    p1x = p1[0]
    p1y = p1[1]
    p2x = p2[0]
    p2y = p2[1]
    ang = np.arctan((p1y - p2y)/(p1x - p2x))
    return ang/np.pi*180

if __name__ == "__main__":
    # 1. All the input configuration
    cluster_num = []
    d = []
    ang = []
    filename = str(input('Enter the result file name:\n'))
    min_samples = int(input('Enter min_samples (you can try to start with 10):\n'))
    eps = float(input('Enter eps (you can try to start with 1):\n'))

    # 2. Clustering process
    # 2.1 Read in data
    print ('Processing ' + filename + '...')
    Data_c = Table.read(filename, format="ascii")
    Data_vxy = np.zeros(shape=(len(Data_c),2))
    Data_vxy[:,0]=Data_c['x'].data
    Data_vxy[:,1]=Data_c['y'].data
    Flux = Data_c['F'].data
    # 2.2 Determine atom_max
    n = 0
    col_num = []
    n_old = Data_c['k'][0]
    for n_new in Data_c['k'][1:]:
        if n_new==n_old:
            n += 1
        else:
            col_num += [n+1]
            n = 0
        n_old = n_new
    atom_max = int(np.max(col_num))
    print ('atom_max', atom_max)

    ## 2.2 Clustering process ##
    atom, xy, flux, noise, labels, centers, widths = find_center(Data_vxy, Flux, atom_max, min_samples, eps)
    if atom > atom_max:
        print ('You need to increase the variable min_samples or/and increase the variable eps')
        print (atom, 'is even bigger than the possible biggest atom number', atom_max)
    again = str(input('Do you want to change variables and try again? (y/n):'))
    print (again)
    while (again == 'Y' or again == 'y'):
        min_samples = int(input('Enter min_samples:\n'))
        eps = float(input('Enter eps:\n'))
        atom, xy, flux, noise, labels, centers, widths = find_center(Data_vxy, Flux, atom_max, min_samples, eps)
        again = str(input('Do you want to change variables and try again? (y/n):'))
        if (again == 'N' or again == 'n'):
            break
    print ('You chose eps = ', eps, ' and min_samples = ', min_samples, '\n')
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
