
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.dataset import dataset
from msmbuilder.decomposition import tICA
from msmbuilder.preprocessing import RobustScaler
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.lumping import PCCAPlus
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump
import shutil


def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr", "--trajectory", type=str, help="path of trr")
    argparser.add_argument("-pdb", "--protein", type=str, help="path of pdb")
    argparser.add_argument("-prob", "--probtxt", type=str, help="path of pdb")
    argparser.add_argument("-name", "--name_kei", type=str, help="name of file")
    return argparser.parse_args()

def mak_dihe():
    trr = ["/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/AFF/TRR/",
          "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/EAF/TRR/",
          "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/TAF/TRR/"]
    pdb = ["/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb",
           "/lustre7/home/lustre3/satoshi/MED/eaf1/HEN.pdb",
           "/lustre7/home/lustre3/satoshi/MED/taf7/HEN.pdb"]
    dir2 = ["aff2_diheds/","eaf2_diheds/","taf2_diheds/"]
    sca_dir = ["aff2_scaled_diheds/","eaf2_scaled_diheds/","taf2_scaled_diheds/"]
    clu_dir = ["aff2_kmeans/","eaf2_kmeans/","taf2_kmeans/"]
    tic = ["aff2_tica/","eaf2_tica/","taf2_tica/"]
    name = ["aff2", "eaf2", "taf2"]
    for i in range(3):
        xyz = dataset(trr[i] + "*.trr",
                      topology=pdb[i],
                      stride=2)
        featurizer = DihedralFeaturizer(types=['phi', 'psi'])
        diheds = xyz.fit_transform_with(featurizer, dir2[i], fmt='dir-npy')
        scaler = RobustScaler()
        scaled_diheds = diheds.fit_transform_with(scaler,
                                                  sca_dir[i],
                                                  fmt='dir-npy')
        # f = open("{0}_tica.txt".format(name[i]))
        for j in range(1,1001):
            print(name[i], "lagtime:","j")
            tica_model = tICA(lag_time=j, n_components=4)
            tica_model = scaled_diheds.fit_with(tica_model)
            print(tica_model.score(diheds))
            # a = "lag_TIME:" + str(j) + str(tica_model.score(diheds)) + "\n"
            # f.write(str(tica_model.score(diheds)))
            # f.write("\n")
            # f.close()
            tica_trajs = scaled_diheds.transform_with(tica_model, tic[i], fmt='dir-npy')
            clusterer = MiniBatchKMeans(n_clusters=20, random_state=42)
            clustered_trajs = tica_trajs.fit_transform_with(
                clusterer, clu_dir[i], fmt='dir-npy'
            )
            # f = open("{0}_msm.txt".format(name[i]))
            for j in range(1,2):
                msm = MarkovStateModel(lag_time=j, n_timescales=20)
                msm.fit(clustered_trajs)
                print(msm.score(clustered_trajs))
                # a = "lag_time" + str(j) + " " + str(msm.score(clustered_trajs))
                # f.write(a)
                # f.close()
            shutil.rmtree(clu_dir[i])
            shutil.rmtree(tic[i])

def main():
    mak_dihe()


if __name__ == '__main__':
    main()
