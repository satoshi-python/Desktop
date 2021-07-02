
# load data
from msmbuilder.dataset import dataset
ds = dataset("tica_trajs.h5")

# plot histogram
import msmexplorer as msme
import numpy as np
txx = np.concatanate(ds)
f = msme.plot_histogram(txx)
f.savefig("histogram.pdf")

# plot free energy landscape
from msmbuilder.utils import load
msm = load('msm.pkl')
clusterer = load('clusterer.pkl')
assignments = clusterer.partial_transform(txx)
assignments = msm.partial_transform(assignments)

from matplotlib import pyplot as plt
msme.plot_free_energy(txx, obs=(0, 1), n_samples=10000,
                      pi=msm.populations_[assignments],
                      xlabel='tIC 1', ylabel='tIC 2')

plt.scatter(clusterer.cluster_centers_[msm.state_labels_, 0],
            clusterer.cluster_centers_[msm.state_labels_, 1],
            s=1e4 * msm.populations_,       # size by population
            c=msm.left_eigenvectors_[:, 1], # color by eigenvector
            cmap="coolwarm",
            zorder=3
           )
plt.colorbar(label='First dynamical eigenvector')
f = plt.tight_layout()
f.savefig("energy_landscape.pdf")
