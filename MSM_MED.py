
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.dataset import dataset
from msmbuilder.decomposition import tICA
from msmbuilder.preprocessing import RobustScaler
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.lumping import PCCAPlus
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump
import numpy as np
import mdtraj
import glob


def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr", "--trajectory", type=str, help="path of trr")
    argparser.add_argument("-pdb", "--protein", type=str, help="path of pdb")
    argparser.add_argument("-prob", "--probtxt", type=str, help="path of pdb")
    argparser.add_argument("-name", "--name_kei", type=str, help="name of file")
    return argparser.parse_args()

def mak_dihe():
    trr = ["/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/AFF/re_name/",
          "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/EAF/re_name/",
          "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/TAF/re_name/"]
    pdb = ["/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb",
           "/lustre7/home/lustre3/satoshi/MED/eaf1/HEN.pdb",
           "/lustre7/home/lustre3/satoshi/MED/taf7/HEN.pdb"]
    dir2 = ["aff2_diheds/","eaf2_diheds/","taf2_diheds/"]
    sca_dir = ["aff2_scaled_diheds/","eaf2_scaled_diheds/","taf2_scaled_diheds/"]
    clu_dir = ["aff2_kmeans/","eaf2_kmeans/","taf2_kmeans/"]
    tic = ["aff2_tica/","eaf2_tica/","taf2_tica/"]
    name = ["aff2", "eaf2", "taf2"]
    num = [260, 5, 94]
    for i in range(1,3):
        TRR = sorted(glob.glob(trr[i] + "*.trr"))
        TTTT = []
        for i1 in TRR:
            xyz = dataset(i1,
                        topology=pdb[i],
                        stride=2)
            TTTT.append(xyz)
        TTTT = mdtraj.join(TTTT)
        print(xyz)
        xyz = TTTT
        print(xyz)
        featurizer = DihedralFeaturizer(types=['phi', 'psi'])
        diheds = xyz.fit_transform_with(featurizer, dir2[i], fmt='dir-npy')
        scaler = RobustScaler()
        scaled_diheds = diheds.fit_transform_with(scaler,
                                                  sca_dir[i],
                                                  fmt='dir-npy')
        print(num[i], " :lagtime")
        tica_model = tICA(lag_time=num[i], n_components=4)
        tica_model = scaled_diheds.fit_with(tica_model)
        print(tica_model.score(diheds))
        tica_trajs = scaled_diheds.transform_with(tica_model, tic[i], fmt='dir-npy')
        clusterer = MiniBatchKMeans(n_clusters=20, random_state=42)
        clustered_trajs = tica_trajs.fit_transform_with(
            clusterer, clu_dir[i], fmt='dir-npy'
        )
        msm = MarkovStateModel(lag_time=1, n_timescales=20)
        msm.fit(clustered_trajs)
        print(msm.score(clustered_trajs))
        txx = np.concatenate(ds)
        assignments = clusterer.partial_transform(txx)
        assignments = msm.partial_transform(assignments)


def main():
    mak_dihe()


if __name__ == '__main__':
    main()
