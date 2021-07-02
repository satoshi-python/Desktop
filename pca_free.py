
import pandas as pd
import seaborn as sns
import statistics
import math
import numpy as np
import matplotlib
import itertools
import os
import tool
import shutil
matplotlib.use('Agg')
import matplotlib.pyplot as plt

path = "/home/satoshi/ALL_PCA/txt/"
file_txt = ["aff.txt", "eaf.txt", "taf.txt",
            "affkai.txt", "eafkai.txt", "tafkai.txt"]


class FREE_ENERGY_LANDSCAPE:
    def ____init__(self):
        self.MED = []
        self.num_CAL_ENERGY = 0
        self.prob = [0]
        self.num_prob = []
        self.num_try1 = []
        self.num_try2 = []
        self.kai = []
        self.x = []
        self.y = []

    def REAS_result_txt(self, txt):
        x = []
        y = []
        kai = []
        for i in open(txt):
            f = i.split()
            x.append(float(f[0]))
            y.append(float(f[1]))
        self.x += x
        self.y += y
        self.prob_num_kai = len(x)
        kai.append(x)
        kai.append(y)
        self.MED.append(kai)
    # calculation prob

    def CAL1(self, prob_txt="a"):
        if prob_txt == "a":
            self.num_prob += [float(1/self.num_prob)] * self.num_prob
        else:
            num1 = 0
            num2 = 0
            for i in open(prob_txt):
                if float(i) != 0.0:
                    self.num_prob.append(float(i))
                    self.num_try1.append(num1)
                    self.num_try2.append(num2)
                    num1 += 1
                num2 += 1
            a = int(self.prob[len(self.prob) - 1]) + num2
            self.prob.append(a)
    # calcuklation free energy

    def CAL2(self):
        # kai = np.array(self.MED)
        num1 = 0
        for i in self.MED:
            num1 += len(i[0])
        print(num1, len(self.num_prob))
        if self.num_CAL_ENERGY == 1:
            print("We have already calcuration")
            self.num_CAL_ENERGY += 1
        elif self.num_CAL_ENERGY == 2:
            print("restart FREE_ENERGY_LANDSCAPE")
        elif num1 != len(self.num_prob):
            print("tne prob number different from TRR file")
        else:
            self.num_CAL_ENERGY += 1
            num2 = 0
            for i in range(len(self.MED)):
                for j in range(len(self.MED[i][0])):
                    self.kai.append(math.log(self.num_prob[num2])*(-1.986)*0.3)
                    # log(prob)*R*T
                    num2 += 1
            del self.MED

    def CAL3(self, kizami):
        self.kizami = int(kizami)
        # self.clt = [[0.0 for i in range(kizami)]for j in range(kizami)]
        # print(self.clt.type())
        self.clt = [[float(0)]*self.kizami]*self.kizami
        self.clt = np.array(self.clt).tolist()
        # print(self.clt.type())
        a = [[[0]]*self.kizami]*self.kizami
        self.clustering_fram = np.array(a).tolist()
        """
        for i in range(kizami):
            x3 = (x_haba * i) + x_min
            x.append(x3)
            y3 = (y_haba + i) + y_min
            y.append(y3)
        """
        for num_frame in self.early_calculate_frame():
            pass
        """
        for i in self.clt:
            print(i)
        for j6 in range(len(self.num_prob)):
            num_x = kizami - (int((self.y[j6] - y_min)/y_haba)+1)
            num_y = (int((self.x[j6] - x_min)/x_haba) - 1)
            # prob_value*-1.986*temperature/1000
            clt[num_x][num_y] += math.log(float(self.kai[j6])) * (-1.986) * 0.3
            self.clustering_fram[num_x][num_y].append(j6)
            print(j6, "/", len(self.num_prob))
        """
        # clt = np.array(clt) - float(np.amin(np.array(clt)))
        # self.clt = clt.tolist()
        for num_cal in self.early_calculate_clt():
            pass
        """
        print("?????????????????????????????????????????????????")
        for i in self.clt:
            print(i)
        for i, j in itertools.product(range(self.kizami), range(self.kizami)):
            print(i, j)
            if clt[i][j] - float(np.amin(np.array(clt))) < 15:
                clt[i][j] -= float(np.amin(np.array(clt)))
            else:
                clt[i][j] = -410
        """

    def early_calculate_frame(self):
        x_min = float(min(self.x))
        x_max = float(max(self.x))
        y_min = float(min(self.y))
        y_max = float(max(self.y))
        x_jiku = (x_max - x_min)
        y_jiku = (y_max - y_min)
        x_haba = x_jiku/self.kizami
        y_haba = y_jiku/self.kizami
        for j6 in range(len(self.num_prob)):
            num_x = self.kizami - (int((self.y[j6] - y_min)/y_haba)+1)
            num_y = (int((self.x[j6] - x_min)/x_haba) - 1)
            # prob_value*-1.986*temperature/1000
            self.clt[num_x][num_y] = self.clt[num_x][num_y] + self.kai[j6]
            self.clustering_fram[num_x][num_y].append(j6)
            # print(self.kai[j6])
            yield

    def early_calculate_clt(self):
        kai = []
        for i, j in itertools.product(range(self.kizami), range(self.kizami)):
            if self.clt[i][j] > 0:
                kai.append(float(self.clt[i][j]))
            yield
        a = max(kai)
        print(a)
        del kai
        for i, j in itertools.product(range(self.kizami), range(self.kizami)):
            if self.clt[i][j] == 0:
                self.clt[i][j] = -2
            else:
                self.clt[i][j] = abs(self.clt[i][j] - a)
            # print(self.clt[i][j])
            yield

    def draw_map(self, file_path="heat_sample.png", title="heat map",
                 x_label="x", y_label="y"):
        MAX = np.amax(np.array(self.clt))
        # MIN = np.amin(np.array(self.clt))
        df = pd.DataFrame(self.clt,
                          index=range(1, int(self.kizami) + 1),
                          columns=range(1, int(self.kizami) + 1))
        clt1 = []
        x_min = float(min(self.x))
        x_max = float(max(self.x))
        y_min = float(min(self.y))
        y_max = float(max(self.y))
        x_jiku = (x_max - x_min)
        y_jiku = (y_max - y_min)
        x_haba = x_jiku/self.kizami
        y_haba = y_jiku/self.kizami
        for i in self.clt:
            for j in i:
                clt1.append(j)
        df_mask = (df < -1)
        cen = statistics.median(clt1) + 1
        ax = sns.heatmap(df, cmap="jet", vmin=0,
                         vmax=MAX, mask=df_mask, center=cen,
                         cbar_kws={'label': str(title)})
        ax.grid()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        # plt.figure(figsize=(20,10))
        plt.xticks([0, self.kizami/4, self.kizami/2, self.kizami*3/4,
                   self.kizami],
                   [round(x_min, 2), round(x_min+(x_haba * self.kizami/4), 2),
                   round((x_max+x_min)/2, 2),
                   round(x_max-(x_haba * self.kizami/4), 2), round(x_max, 2)])
        plt.yticks([0, self.kizami/4, self.kizami/2, self.kizami*3/4,
                    self.kizami],
                   [round(y_max, 2), round(y_max-(y_haba * self.kizami/4), 2),
                   round((y_max+y_min)/2, 2),
                   round(y_min+(y_haba * self.kizami/4), 2), round(y_min, 2)])
        plt.savefig("{0}".format(file_path), bbox_inches="tight")
        plt.close()

    def pick_up(self, trr_path, pdb_path, save_path="."):
        dic = {}
        """
        frm = []
        for i in trr_path_list:
            u = MDAnalysis.Universe(i)
            frm += u.trajectory
        """
        for i, j in itertools.product(range(self.kizami), range(self.kizami)):
            if self.clt[i][j] < 5 and self.clt[i][j] >= 0:
                KKK = []
                num1 = 0
                for kkk in self.clustering_fram[i][j]:
                    if kkk != 0:
                        KKK.append(self.num_try2[self.num_try1.index(kkk)])
                    num1 += 1
                dic[i, j] = KKK
        del self.kizami, self.clt, self.clustering_fram
        for i, f in dic.items():
            f = [int(j) for j in f]
            print(f)
            for h in trr_path:
                if len(f) > 0:
                    if f[0] < self.prob[trr_path.index(h) + 1]:
                        c = [int(jj) for jj in f if jj < self.prob[trr_path.index(h) + 1]]
                        print(trr_path)
                        tool.MD_to_pdb_chain(trr_path,
                                             pdb_path, c)
                        for ii in range(len(c)):
                            f.pop(0)
                """
                if f[0] > self.prob[j + 1]:
                    f = [int(k - 900000) for k in range(f)]
                if f[0] > 900000 and j < 3:
                    f = [int(k - 900000) for k in range(f)]
                elif f[0] > 1000000 and j > 2:
                    f = [int(k - 1000000) for k in range(f)]
                for j1 in range(len(f)):
                    if f[j1] < 0:
                        f.pop(j1)
                """
            os.rename("pdb_file", str(i))
            shutil.move(str(i), save_path)
        del self
