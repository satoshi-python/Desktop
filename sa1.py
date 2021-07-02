
import numpy as np
import matplotlib
import seaborn as sns
import pandas as pd
import sys
import itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt

args = sys.argv


def NAME():
    filename = args[1]
    name = filename[:-4]
    x_a = []
    y_a = []
    z_a = []
    x_b = []
    y_b = []
    z_b = []
    for i in open(filename):
        f = i.split()
        if f[2] == "CA":
            f = i.split()
            if f[4] == "A":
                x_a.append(float(f[6]))
                y_a.append(float(f[7]))
                z_a.append(float(f[8]))
            if f[4] == "B":
                x_b.append(float(f[6]))
                y_b.append(float(f[7]))
                z_b.append(float(f[8]))
    CC = [[0 for i in range(len(x_b))]for j in range(len(x_a))]
    for i in range(len(x_a)):
        for j in range(len(x_b)):
            A = np.array((x_a[i], y_a[i], z_a[i]), dtype=float)
            B = np.array((x_b[j], y_b[j], z_b[j]), dtype=float)
            CC[i][j] += np.linalg.norm(A - B)
    return name, CC


def main():
    # CM = [[0 for i in range(92)]for j in range(22)]
    # CM1 = []
    OM, CM1 = NAME()
    CM = list(itertools.chain.from_iterable(CM1))
    MAX = np.amax(np.array(CM))
    MIN = np.amin(np.array(CM))
    print(MIN)
    # print(MAX,MIN)
    df = pd.DataFrame(CM1,
                      index=range(1, len(CM1)+1),
                      columns=range(1, len(CM1[0]) + 1))
    # df_mask = (df <= statistics.median(CM))
    ax = sns.heatmap(df, cmap="jet", vmin=MIN, vmax=MAX,
                     cbar_kws={'label': 'Angstrom'})
    # jet or OrRd
    ax.grid()
    plt.xlabel("residue of B")
    plt.ylabel("residue of A")
    # plt.figure(figsize=(20,10))
    plt.show()
    plt.savefig("photo/{0}.png".format(OM))
    plt.close()


if __name__ == "__main__":
    main()
