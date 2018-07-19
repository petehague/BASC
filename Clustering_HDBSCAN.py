'''
   Author: Haoyang Ye (Zoe)
'''

import numpy as np
import hdbscan
from astropy.table import Table

def atom_finder(Data, cluster_size, samples):
    """
    Find the cluster centers for given data using HDBSCAN package.

    Data: input data, it should be a n*2 array
    cluster_size: sets the value 'min_cluster_size' in function hdbscan,
        which is a relatively intuitive parameter to select, the
        bigger the value is, the less clusters hdbscan will pick
    samples: sets the value 'min_samples' in function hdbscan, which determines
        the minimum grouping sample number
    atom: the number of clusters picked by hdbscan
    xy: data arranged according to clustering results
    noise: data that are not selected into any clusters
    labels: flags/labels that determine which data belongs to which cluster
    centers: cluster centre location (x_i,y_i) of each cluster i
    """
    hdb = hdbscan.HDBSCAN(min_cluster_size = cluster_size, min_samples = samples,allow_single_cluster = True).fit(Data_vxy)
    labels = hdb.labels_
    atom = hdb.labels_.max() + 1
    unique_labels = set(labels)
    noise_label = (labels == -1)
    noise = Data[noise_label]
    xy = []
    centers = np.zeros((atom,2))
    widths = np.zeros((atom,2))
    for k in range(atom):
        class_member_mask = (labels == k)
        xy += [Data[class_member_mask]]
        centers[k] = np.mean(Data[class_member_mask], axis = 0)
        widths[k] = np.std(Data[class_member_mask], axis = 0)
    return atom, xy, noise, labels, centers, widths

def combine_atoms(centers, xy, atom):
    """
    If two cluster centres are in the same pixel, combine the two cluters as one cluster.
    """
    centers_int = np.floor(centers)
    xy_updated_list = []
    not_list = []
    for i in range(atom):
        xy_updated = xy[i]
        if i not in not_list:
            for j in range(i+1,atom):
                if (j not in not_list):
                    if np.array_equal(centers_int[j],centers_int[i]):
                        not_list += [j]
                        xy_updated = np.append(xy_updated,xy[j],axis = 0)
            xy_updated_list += [xy_updated]
    atom_updated = len(xy_updated_list)
    centers_updated = np.zeros((atom_updated,2))
    widths_updated = np.zeros((atom_updated,2))
    for i in range(atom_updated):
        centers_updated[i] = np.mean(xy_updated_list[i], axis = 0)
        widths_updated[i] = np.std(xy_updated_list[i], axis = 0)
    return atom_updated, xy_updated_list, centers_updated, widths_updated

def find_center_atom(Data, atom_max, cluster_size, samples):
    """
    Adjust the parameter 'cluster_size' to find cluster centers using
    function 'atom_finder'. If the atom number is found to be zero, 'cluster_size'
    needs to be decreased. If the atom number is more than 'atom_max',
    'cluster_size' needs to be increased.

    Data: input data, it should be a n*2 array
    atom_max: the MCMC results determine the maximum number of the atoms that can be
        found. atom_max helps adjust 'cluster_size', if the atom number is more
        than 'atom_max','cluster_size' needs to be increased.
    cluster_size: sets the value 'min_cluster_size' in function hdbscan,
        which is a relatively intuitive parameter to select, the
        bigger the value is, the less clusters hdbscan will pick
    samples: sets the value 'min_samples' in function hdbscan, which determines
        the minimum grouping sample number
    """
    atom, xy, noise, labels, centers, widths = atom_finder(Data, cluster_size, samples)
    atom, xy, centers, widths = combine_atoms(centers, xy, atom)
    for i in range(len(Data)):
        if atom > atom_max:
            cluster_size = cluster_size + 1
            print (cluster_size)
            atom, xy, noise, labels, centers, widths = atom_finder(Data, cluster_size, samples)
            atom, xy, centers = combine_atoms(centers, xy, atom)
        elif atom == 0:
            cluster_size = cluster_size - 1
            atom, xy, noise, labels, centers, widths = atom_finder(Data, cluster_size, samples)
            atom, xy, centers = combine_atoms(centers, xy, atom)
        else:
            break
    if len(noise) == 0:
        print ('There is no outlier.')
    else:
        print ('There are', len(noise), 'outliers. Please have a look at them.')
    print ('We find', atom, 'atoms, with min_cluster_size = ', cluster_size, '\n')
    print ('The', atom, 'cluster centers are \n', centers)
    print ('with a standard deviation of ', widths, '\n')
    return atom, xy, noise, labels, centers, widths


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
    import matplotlib.pyplot as plt
    # 1. All the input configuration
    cluster_num = []
    d = []
    ang = []
    filename = str(input('Enter the MCMC result file name:\n'))
    samples = int(input('Enter the minimum group number (you can try to start with 10):\n'))
    cluster_size = float(input('Enter the min_cluster_size (you can try to start with 30):\n'))

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
    # 2.3 Start HDBSCAN algorithm
    atom, xy, noise, labels, centers, widths = find_center_atom(Data_vxy, atom_max, int(cluster_size), int(samples))
    print ('We find', atom, 'atoms\n')
    # 2.4 Plot the results and save ############
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
    figname = filename + '_HDBSCAN_clustering_results.png'
    plt.savefig(figname)

    if atom == 2:
        d += [distance(centers[0], centers[1])]
        ang += [angle_between(centers[0], centers[1])]
        cluster_num += [len(centers)]
    elif atom >= 3:
        d += [100]
        ang += [100]
        cluster_num += [len(centers)]
    elif atom == 1:
        d += [0]
        ang += [0]
        cluster_num += [len(centers)]
    else:
        d += [-100]
        ang += [-100]
        cluster_num += [len(centers)]

    plt.close('all')
