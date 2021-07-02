import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from argparse import ArgumentParser
import glob
pyemma.config.show_progress_bars = 'True'


def project_and_cluster(trajfiles, featurizer, sparsify=False, tica=True, lag=100, scale=True, var_cutoff=0.95, ncluster=100):
    """
    Returns
    -------
    trans_obj, Y, clustering

    """
    X = coor.load(trajfiles, featurizer) #str, list of str or nested list (one level) of str
    if sparsify:
        X = remove_constant(X)
    if tica:
        trans_obj = coor.tica(X, lag=lag, var_cutoff=var_cutoff)
    else:
        trans_obj = coor.pca(X, dim=-1, var_cutoff=var_cutoff)
    Y = trans_obj.get_output()
    if scale:
        for y in Y:
            y *= trans_obj.eigenvalues[:trans_obj.dimension()]
    cl_obj = coor.cluster_kmeans(Y, k=ncluster, max_iter=3, fixed_seed=True)
    return trans_obj, Y, cl_obj

def remove_constant(X, threshold=0.001):
    if isinstance(X, np.ndarray):
        X = [X]
    Ds = [np.max(x, axis=0) - np.min(x, axis=0) for x in X]
    D = np.min(np.array(Ds), axis=0)
    Ivar = np.where(D > 0.001)[0]
    Y = [x[:, Ivar] for x in X]
    if len(Y) == 1:
        Y = Y[0]
    return Y

def eval_transformer(trans_obj):
    # Effective dimension (Really? If we just underestimate the Eigenvalues this value also shrinks...))
    print('Evaluating transformer: ', str(trans_obj.__class__))
    print('effective dimension', np.sum(1.0 - trans_obj.cumvar))
    print('eigenvalues', trans_obj.eigenvalues[:5])
    print('partial eigensum', np.sum(trans_obj.eigenvalues[:10]))
    print('total variance', np.sum(trans_obj.eigenvalues ** 2))
    print()

def plot_map(Y, sx=None, sy=None, tickspacing1=1.0, tickspacing2=1.0, timestep=1.0, timeunit='ns'):
    if not isinstance(Y, np.ndarray):
        Y = Y[0]
    if sx is None:
        sx = -np.sign(Y[0,0])
    if sy is None:
        sy = -np.sign(Y[0,1])
    Y1 = sx*Y[:, 0]
    min1 = np.min(Y1)
    max1 = np.max(Y1)
    Y2 = sy*Y[:, 1]
    min2 = np.min(Y2)
    max2 = np.max(Y2)
    # figure
    figure(figsize=(16,4))
    # trajectories
    subplot2grid((2,2), (0,0))
    plot(timestep*np.arange(len(Y1)), Y1)
    xlim(0, timestep*len(Y1))
    yticks(np.arange(int(min1), int(max1)+1, tickspacing1))
    ylabel('component 1')
    subplot2grid((2,2), (1,0))
    plot(timestep*np.arange(len(Y2)), Y2)
    xlim(0, timestep*len(Y2))
    ylabel('component 2')
    yticks(np.arange(int(min2), int(max2)+1, tickspacing2))
    xlabel('time / ' + timeunit)
    # histogram data
    subplot2grid((2,2), (0,1), rowspan=2)
    z,x,y = np.histogram2d(Y1, Y2, bins=50)
    z += 0.1
    # compute free energies
    F = -np.log(z)
    # contour plot
    extent = [x[0], x[-1], y[0], y[-1]]
    xticks(np.arange(int(min1), int(max1)+1, tickspacing1))
    yticks(np.arange(int(min2), int(max2)+1, tickspacing2))
    contourf(F.T, 50, cmap=plt.cm.nipy_spectral, extent=extent)
    xlabel('component 1')
    ylabel('component 2')


def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-name", "--name_kei", type=str, help="name of file")
    return argparser.parse_args()


def PATH_PDB(name):
    if name == "aff":
        path = "/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb"
    elif name == "eaf":
        path = "/lustre7/home/lustre3/satoshi/MED/eaf1/HEN.pdb"
    elif name == "taf":
        path = "/lustre7/home/lustre3/satoshi/MED/taf7/HEN.pdb"
    return path


def PATH_TRR(name):
    if name == "aff":
        TRR = "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/AFF/re_name/"
    elif name = "eaf":
        TRR = "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/EAF/re_name/"
    elif name = "taf":
        TRR = "/lustre7/home/lustre3/satoshi/cluster/prepare_msm/tsubame/TAF/re_name/"
    path = sorted(glob.glob(TRR + "*.trr"))
    return path


def main():
    args = get_option()
    top = PATH_PDB(str(args.name_kei))
    trajs = PATH_TRR(str(args.name_kei))
    feat_Cadist = coor.featurizer(top) #top:pdb-file
    feat_Cadist.add_distances_ca()
    tica_Cadist, tica_Y_Cadist, tica_cl_Cadist = project_and_cluster(trajs, feat_Cadist) #trajs:trr-file
    eval_transformer(tica_Cadist)
    pca_Cadist, pca_Y_Cadist, pca_cl_Cadist = project_and_cluster(trajs, feat_Cadist, tica=False)
    eval_transformer(pca_Cadist)
    feat_Cadist = coor.featurizer(top) #top:pdb-file
    feat_Cadist.add_distances_ca()
    tica_Cadist, tica_Y_Cadist, tica_cl_Cadist = project_and_cluster(trajs, feat_Cadist) #trajs:trr-file
    eval_transformer(tica_Cadist)
    pca_Cadist, pca_Y_Cadist, pca_cl_Cadist = project_and_cluster(trajs, feat_Cadist, tica=False)
    eval_transformer(pca_Cadist)
    feat_Cainvdist = coor.featurizer(top)
    pairs = feat_Cainvdist.pairs(feat_Cainvdist.select_Ca())
    feat_Cainvdist.add_inverse_distances(pairs)
    feat_vc0_75 = coor.featurizer(top)
    feat_vc0_75.add_residue_mindist(threshold=0.75)
    tica_vc0_75, tica_Y_vc0_75, tica_cl_vc0_75 = project_and_cluster(trajs, feat_vc0_75, sparsify=True)
    eval_transformer(tica_vc0_75)
    pca_vc0_75, pca_Y_vc0_75, pca_cl_vc0_75 = project_and_cluster(trajs, feat_vc0_75, sparsify=True, tica=False)
    eval_transformer(pca_vc0_75)
    labels = ['Ca coords', 'Ca dists', 'contacts 0.75']
    ticas = [tica_Ca, tica_Cadist, tica_vc0_75]
    tica_Ys = [tica_Y_Ca, tica_Y_Cadist, tica_Y_vc0_75]
    tica_cls = [tica_cl_Ca, tica_cl_Cadist, tica_cl_vc0_75]
    pcas = [pca_Ca, pca_Cadist, pca_vc0_75]
    pca_Ys = [pca_Y_Ca, pca_Y_Cadist, pca_Y_vc0_75]
    pca_cls = [pca_cl_Ca, pca_cl_Cadist, pca_cl_vc0_75]
    for i,Y in enumerate(pca_Ys):
        plot_map(Y, tickspacing1=2.0, tickspacing2=2.0, timestep=0.01, timeunit='$\mu$s')
        gcf().suptitle('PCA ' + labels[i], fontsize=20)
        savefig('figs/pca_' + labels[i].replace(' ','_') + '.png', bbox_inches='tight')
    for i,Y in enumerate(tica_Ys):
        plot_map(Y, tickspacing1=2.0, tickspacing2=1.0, timestep=0.01, timeunit='$\mu$s')
        gcf().suptitle('TICA ' + labels[i], fontsize=20)
        savefig('figs/tica_' + labels[i].replace(' ','_') + '.png', bbox_inches='tight')
    for i, tica_obj in enumerate(ticas):
        plot(np.cumsum(tica_obj.eigenvalues**2), linewidth=2, label=labels[i])
    text(1.2, 13, 'TICA ', fontweight='bold')
    xlabel('eigenvalue index')
    ylabel('total kinetic variance')
    legend(fontsize=14, loc=4)
    semilogx()
    savefig('figs/tica_totvar.png', bbox_inches='tight')
    for i, pca_obj in enumerate(ticas):
        plot(np.cumsum(pca_obj.eigenvalues), linewidth=2, label=labels[i], linestyle='dashed')
    text(1.2, 13, 'PCA ', fontweight='bold')
    xlabel('eigenvalue index')
    ylabel('total variance')
    legend(fontsize=14, loc=4)
    semilogx()
    savefig('figs/pca_totvar.png', bbox_inches='tight')
    maxlag = 500
    tica_itss = [msm.timescales_msm(cl_obj.dtrajs, lags=maxlag, nits=5, errors='bayes') for cl_obj in tica_cls]
    pca_itss = [msm.timescales_msm(cl_obj.dtrajs, lags=maxlag, nits=5, errors='bayes') for cl_obj in pca_cls]
    fig, axes = subplots(3, figsize=(5,12))
    #axes = axes.flatten()
    for i, its_obj in enumerate(tica_itss):
        mplt.plot_implied_timescales(its_obj, ax=axes[i], ylog=True, dt=0.01, units='$\mu$s')
        axes[i].text(0.2, 500, 'TICA ' + labels[i], fontweight='bold')
        axes[i].set_ylim(1, 1000)
    savefig('figs/tica_itss.png', bbox_inches='tight')
    fig, axes = subplots(3, figsize=(5,12))
    #axes = axes.flatten()
    for i, its_obj in enumerate(pca_itss):
        mplt.plot_implied_timescales(its_obj, ax=axes[i], ylog=True, dt=0.01, units='$\mu$s')
        axes[i].text(0.2, 500, 'PCA ' + labels[i], fontweight='bold')
        axes[i].set_ylim(1, 1000)
    savefig('figs/pca_itss.png', bbox_inches='tight')
    lags = tica_itss[0].lags
    colors = ['blue', 'green', 'red']
    kinvars = [[(M.eigenvalues()**2).sum() for M in its_obj.models] for its_obj in tica_itss]
    for i, kinvar in enumerate(kinvars):
        plot(0.01*lags, kinvar, linewidth=2, label='TICA '+labels[i], color=colors[i])
    kinvars = [[(M.eigenvalues()**2).sum() for M in its_obj.models] for its_obj in pca_itss]
    for i, kinvar in enumerate(kinvars):
        plot(0.01*lags, kinvar, linewidth=2, label='PCA '+labels[i], color=colors[i], linestyle='dashed')
    text(0.2, 40, 'MSM ', fontweight='bold')
    xlabel('lag time / $\mu$s')
    ylim(0, 45); ylabel('total kinetic variance')
    legend(fontsize=14)
    savefig('figs/total_kinetic_variance.png', bbox_inches='tight')


if __name__ == '__main__':
    main()
