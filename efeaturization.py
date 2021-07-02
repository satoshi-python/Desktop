
import MDAnalysis
import mdtraj as md
import numpy as np
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from multiprocessing import Pool
from msmbuilder.preprocessing import RobustScaler
from msmbuilder.decomposition import tICA

## load
def load():
    ## load
    frm = md.load("/home/xieqilin/test/test_all_0.1ns.xtc",top="/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb")
    #u = MDAnalysis.Universe("/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb","/lustre7/home/lustre3/satoshi/MED/aff4/test_all.trr")
    #frm = u.trajectory
    #xyz = [t[::10] for t in u]
    print(frm)
    return frm

#Featurization
def feat(xyz):   ##Featurize
    featurizer = DihedralFeaturizer(types=['phi', 'psi'])
    diheds = featurizer.fit_transform(xyz)
    print(xyz[0].xyz.shape)
    print(diheds[0].shape)
    return xyz, diheds

# Preprocessing
def preprocessing(diheds):
    scaler = RobustScaler()
    scaled_diheds = scaler.fit_transform(diheds)
    print(scaled_diheds[0].shape)
    return scaled_diheds

# Intermediate kinetic model: tICA
def tica(diheds):
    # fit and transform can be done in seperate steps:
    tica_model = tICA(n_components=2, lag_time=0)
    tica_model.fit(diheds)
    tica_trajs = tica_model.transform(diheds)
    print(tica_trajs[0].shape)
    return tica_trajs

def main():
    frm = load()
    xyz, diheds = feat(frm)
    scaled_diheds = preprocessing(diheds)
    tica_trajs = tica(diheds)

if __name__ == "__main__":
    main()
