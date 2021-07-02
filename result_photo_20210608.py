
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
import pandas as pd
import statistics
import math

path_pca_66 = ["/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/aff.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/eaf.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/taf.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/aff2.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/eaf2.txt", "/lustre7/home/lustre3/satoshi/sa/ALL/txt/66/66-youso/taf2.txt"]
path_pca_2024 =  ["/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/aff.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/eaf.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/taf.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/aff2.txt","/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/eaf2.txt", "/lustre7/home/lustre3/satoshi/sa/ALL/txt/2024/2024_youso/taf2.txt"]
probpath=["/lustre7/home/lustre3/satoshi/sa/clt/txt/prob_aff.txt","/lustre7/home/lustre3/satoshi/sa/clt/txt/prob_eaf.txt","/lustre7/home/lustre3/satoshi/sa/clt/txt/prob_taf.txt","/lustre7/home/lustre3/satoshi/sa/clt_kai/dimention66/txt/aff4/prob.dat","/lustre7/home/lustre3/satoshi/sa/clt_kai/dimention66/txt/eaf1/prob.dat","/lustre7/home/lustre3/satoshi/sa/clt_kai/dimention66/txt/taf7/prob.dat"]
name=["MED26NTD(WT)-aff4","MED26NTD(WT)-eaf1","MED26NTD(WT)-taf7","MED26NTD(mutant)-aff4","MED26NTD(mutant)-eaf1","MED26NTD(mutant)-taf7",]
MM = [chr(i) for i in range(97, 97+26)]
num_MM = 0

TTT = path_pca_2024

def ALL(kizami):
    plt.subplot(3,3,1)
    kizami = int(kizami)
    x =[]
    x1 = []
    y1 = []
    for i in range(6):
        x = []
        for j4 in  open("{0}".format(TTT[i])):
            x.append(float(j4))
        for j5 in range(len(x)):
            if j5 < (len(x)/2):
                x1.append(x[j5])
            else:
                y1.append(x[j5])
    x_max = max(x1)
    x_min = min(x1)
    y_max = max(y1)
    y_min = min(y1)
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for j in range(6):
        for i in open("{0}".format(probpath[j])):
            prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    x_haba = x_jiku/kizami
    y_haba = y_jiku/kizami
    clt = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    clt5 = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    x = []
    y = []
    #print(len(x1))
    #print(len(prob_kai))
    for i in range(kizami):
        x3 = (x_haba * i) + x_min
        x.append(x3)
        y3 = (y_haba + i) + y_min
        y.append(y3)
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba)+1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
        #print("end")
        #WRITE(ligand,"pca",kizami,clt,savepath[0])
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
    m = float(min(clt2))
    M = float(max(clt2))
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                    clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
                else :
                    clt5[j6][j7] += (-1)
    plt.rcParams["font.size"] = 20
    df = pd.DataFrame(clt5,
                      index = [y_max - (y_haba * i) for i in range(kizami)],
                      columns = [x_min + (x_haba * j7) for j7 in range(kizami)])
    RE_date = [y_max, y_min, y_haba, x_max, x_min, x_haba]
    df_mask = (df <= 0)
    cen = statistics.median(clt2) + 1
    ax = sns.heatmap(df , cmap="jet",vmin = 0 ,vmax = 10,center = cen, mask = df_mask, cbar_kws={'label': 'PMF'})
    ax.figure.axes[-1].yaxis.label.set_size(FONT)
    ax_pos = ax.get_position()
    # soutaiteki ni iti wo syutoku
    #plt.text(ax_pos.x1 - 2, ax_pos.y1 + 5, r"$\bf{" + str(MM[num_MM]) + "}$", fontsize = 20)
    #num_MM += 1
    #ax.grid()
    #plt.subplots_adjust(wspace=0.5, hspace=0.5)
    #plt.subplot(3,2,j+1)
    cax = ax.collections[0].colorbar.ax
    cax.tick_params(which='major', labelsize=18)
    cax.xaxis.label.set_fontsize(20)
    plt.xlabel("PC1", fontsize=FONT)
    plt.ylabel("PC2", fontsize=FONT)
    plt.title("MED26NTD(WT)\nMED26NTD(mutant)",fontsize=FONT)
    #plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[x_min,x_min+(x_haba * kizami/4),0,x_max-(x_haba * kizami/4),x_min])
    #plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[y_max,y_max-(y_haba * kizami/4),0,y_min+(y_haba * kizami/4),y_min])
    plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(x_min,2),round(x_min+(x_haba * kizami/4),2),round((x_max+x_min)/2,2),round(x_max-(x_haba * kizami/4),2),round(x_max,2)])
    plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(y_max,2),round(y_max-(y_haba * kizami/4),2),round((y_max+y_min)/2,2),round(y_min+(y_haba * kizami/4),2),round(y_min,2)])
    return RE_date


def mutant(kizami, RE_date):
    # RE_date = [y_max, y_min, y_haba, x_max, x_min, x_haba]
    plt.subplot(3,3,3)
    kizami = int(kizami)
    x =[]
    x1 = []
    y1 = []
    for i in range(3,6):
        x = []
        for j4 in  open("{0}".format(TTT[i])):
            x.append(float(j4))
        for j5 in range(len(x)):
            if j5 < (len(x)/2):
                x1.append(x[j5])
            else:
                y1.append(x[j5])
    x_max = RE_date[3]
    x_min = RE_date[4]
    y_max = RE_date[0]
    y_min = RE_date[1]
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for j in range(3,6):
        for i in open("{0}".format(probpath[j])):
            prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    x_haba = x_jiku/kizami
    y_haba = y_jiku/kizami
    clt = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    clt5 = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    x = []
    y = []
    #print(len(x1))
    #print(len(prob_kai))
    for i in range(kizami):
        x3 = (x_haba * i) + x_min
        x.append(x3)
        y3 = (y_haba + i) + y_min
        y.append(y3)
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba)+1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
        #print("end")
        #WRITE(ligand,"pca",kizami,clt,savepath[0])
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
    m = float(min(clt2))
    M = float(max(clt2))
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                    clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
                else :
                    clt5[j6][j7] += (-1)
    plt.rcParams["font.size"] = 20
    df = pd.DataFrame(clt5,
                      index = [y_max - (y_haba * i) for i in range(kizami)],
                      columns = [x_min + (x_haba * j7) for j7 in range(kizami)])
    df_mask = (df <= 0)
    cen = statistics.median(clt2) + 1
    ax = sns.heatmap(df , cmap="jet",vmin = 0 ,vmax = 10,center = cen, mask = df_mask, cbar_kws={'label': 'PMF'})
    ax.figure.axes[-1].yaxis.label.set_size(FONT)
    ax_pos = ax.get_position()
    # soutaiteki ni iti wo syutoku
    #plt.text(ax_pos.x1 - 2, ax_pos.y1 + 5, r"$\bf{" + str(MM[num_MM]) + "}$", fontsize = 20)
    #num_MM += 1
    #ax.grid()
    #plt.subplots_adjust(wspace=0.5, hspace=0.5)
    #plt.subplot(3,2,j+1)
    cax = ax.collections[0].colorbar.ax
    cax.tick_params(which='major', labelsize=18)
    cax.xaxis.label.set_fontsize(20)
    plt.xlabel("PC1", fontsize=FONT)
    plt.ylabel("PC2", fontsize=FONT)
    plt.title("MED26NTD(mutant)",fontsize=FONT)
    #plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[x_min,x_min+(x_haba * kizami/4),0,x_max-(x_haba * kizami/4),x_min])
    #plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[y_max,y_max-(y_haba * kizami/4),0,y_min+(y_haba * kizami/4),y_min])
    plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(x_min,2),round(x_min+(x_haba * kizami/4),2),round((x_max+x_min)/2,2),round(x_max-(x_haba * kizami/4),2),round(x_max,2)])
    plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(y_max,2),round(y_max-(y_haba * kizami/4),2),round((y_max+y_min)/2,2),round(y_min+(y_haba * kizami/4),2),round(y_min,2)])

def WT(kizami, RE_date):
    plt.subplot(3,3,2)
    kizami = int(kizami)
    x =[]
    x1 = []
    y1 = []
    for i in range(3):
        x = []
        for j4 in  open("{0}".format(TTT[i])):
            x.append(float(j4))
        for j5 in range(len(x)):
            if j5 < (len(x)/2):
                x1.append(x[j5])
            else:
                y1.append(x[j5])
    x_max = RE_date[3]
    x_min = RE_date[4]
    y_max = RE_date[0]
    y_min = RE_date[1]
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for j in range(3):
        for i in open("{0}".format(probpath[j])):
            prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    x_haba = x_jiku/kizami
    y_haba = y_jiku/kizami
    clt = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    clt5 = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
    x = []
    y = []
    #print(len(x1))
    #print(len(prob_kai))
    for i in range(kizami):
        x3 = (x_haba * i) + x_min
        x.append(x3)
        y3 = (y_haba + i) + y_min
        y.append(y3)
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba)+1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
        #print("end")
        #WRITE(ligand,"pca",kizami,clt,savepath[0])
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
    m = float(min(clt2))
    M = float(max(clt2))
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                    clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
                else :
                    clt5[j6][j7] += (-1)
    plt.rcParams["font.size"] = 20
    df = pd.DataFrame(clt5,
                      index = [y_max - (y_haba * i) for i in range(kizami)],
                      columns = [x_min + (x_haba * j7) for j7 in range(kizami)])
    df_mask = (df <= 0)
    cen = statistics.median(clt2) + 1
    ax = sns.heatmap(df , cmap="jet",vmin = 0 ,vmax = 10,center = cen, mask = df_mask, cbar_kws={'label': 'PMF'})
    ax.figure.axes[-1].yaxis.label.set_size(FONT)
    ax_pos = ax.get_position()
    # soutaiteki ni iti wo syutoku
    #plt.text(ax_pos.x1 - 2, ax_pos.y1 + 5, r"$\bf{" + str(MM[num_MM]) + "}$", fontsize = 20)
    #num_MM += 1
    #ax.grid()
    #plt.subplots_adjust(wspace=0.5, hspace=0.5)
    #plt.subplot(3,2,j+1)
    cax = ax.collections[0].colorbar.ax
    cax.tick_params(which='major', labelsize=18)
    cax.xaxis.label.set_fontsize(20)
    plt.xlabel("PC1", fontsize=FONT)
    plt.ylabel("PC2", fontsize=FONT)
    plt.title("MED26NTD(WT)",fontsize=FONT)
    #plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[x_min,x_min+(x_haba * kizami/4),0,x_max-(x_haba * kizami/4),x_min])
    #plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[y_max,y_max-(y_haba * kizami/4),0,y_min+(y_haba * kizami/4),y_min])
    plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(x_min,2),round(x_min+(x_haba * kizami/4),2),round((x_max+x_min)/2,2),round(x_max-(x_haba * kizami/4),2),round(x_max,2)])
    plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(y_max,2),round(y_max-(y_haba * kizami/4),2),round((y_max+y_min)/2,2),round(y_min+(y_haba * kizami/4),2),round(y_min,2)])


kizami_hai = np.arange(1,1001)
FONT = 24

for kizami in kizami_hai:
    num_MM = 0
    fig = plt.figure(figsize=(45,45), tight_layout=True)
    plt.subplots_adjust(wspace=1, hspace=1)
    RE_date = ALL(kizami)
    # RE_date = [y_max, y_min, y_haba, x_max, x_min, x_haba]
    WT(kizami, RE_date)
    mutant(kizami, RE_date)
    MM = [chr(i) for i in range(97, 97+26)]
    for j in range(len(TTT)):
        plt.subplot(3,3,j+1+3)
        kizami = int(kizami)
        x =[]
        x1 = []
        y1 = []
        x = []
        for j4 in  open("{0}".format(TTT[j])):
            x.append(float(j4))
        for j5 in range(len(x)):
            if j5 < (len(x)/2):
                x1.append(x[j5])
            else:
                y1.append(x[j5])
        x_max = RE_date[3]
        x_min = RE_date[4]
        y_max = RE_date[0]
        y_min = RE_date[1]
        x_jiku = (x_max - x_min)
        y_jiku = (y_max - y_min)
        x_haba = x_jiku / kizami
        y_haba = y_jiku / kizami
        prob = []
        prob_kai = []
        for i in open("{0}".format(probpath[j])):
            prob.append(float(i))
        for i in range(len(prob)):
            if prob[i] > 0:
                prob_kai.append(prob[i])
        x_haba = x_jiku/kizami
        y_haba = y_jiku/kizami
        clt = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
        clt5 = [[0.0 for i in range(kizami)] for j7 in range(kizami)]
        x = []
        y = []
        #print(len(x1))
        #print(len(prob_kai))
        for i in range(kizami):
            x3 = (x_haba * i) + x_min
            x.append(x3)
            y3 = (y_haba + i) + y_min
            y.append(y3)
        for j6 in range(len(x1)):
            clt[kizami - (int((y1[j6] - y_min)/y_haba)+1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
        #print("end")
        #WRITE(ligand,"pca",kizami,clt,savepath[0])
        clt2 = []
        for j6 in range(kizami):
            for j7 in range(kizami):
                if clt[j6][j7] > 0.0:
                    if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                        clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
        m = float(min(clt2))
        M = float(max(clt2))
        clt2 = []
        for j6 in range(kizami):
            for j7 in range(kizami):
                if clt[j6][j7] > 0.0:
                    if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                        clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                        clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
                    else :
                        clt5[j6][j7] += (-1)
        plt.rcParams["font.size"] = 20
        df = pd.DataFrame(clt5,
                          index = [y_max - (y_haba * i) for i in range(kizami)],
                          columns = [x_min + (x_haba * j7) for j7 in range(kizami)])
        df_mask = (df <= 0)
        cen = statistics.median(clt2) + 1
        ax = sns.heatmap(df , cmap="jet",vmin = 0 ,vmax = 10,center = cen, mask = df_mask, cbar_kws={'label': 'PMF'})
        ax.figure.axes[-1].yaxis.label.set_size(FONT)
        ax_pos = ax.get_position()
        # soutaiteki ni iti wo syutoku
        #plt.text(ax_pos.x1 , ax_pos.y1 , r"$\bf{" + str(MM[num_MM]) + "}$", fontsize = 20)
        #print(num_MM)
        num_MM += 1
        #ax.grid()
        #plt.subplots_adjust(wspace=0.5, hspace=0.5)
        #plt.subplot(3,2,j+1)
        cax = ax.collections[0].colorbar.ax
        cax.tick_params(which='major', labelsize=18)
        cax.xaxis.label.set_fontsize(20)
        plt.xlabel("PC1", fontsize=FONT)
        plt.ylabel("PC2", fontsize=FONT)
        plt.title(name[j],fontsize=FONT)
        #plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[x_min,x_min+(x_haba * kizami/4),0,x_max-(x_haba * kizami/4),x_min])
        #plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[y_max,y_max-(y_haba * kizami/4),0,y_min+(y_haba * kizami/4),y_min])
        plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(x_min,2),round(x_min+(x_haba * kizami/4),2),round((x_max+x_min)/2,2),round(x_max-(x_haba * kizami/4),2),round(x_max,2)])
        plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(y_max,2),round(y_max-(y_haba * kizami/4),2),round((y_max+y_min)/2,2),round(y_min+(y_haba * kizami/4),2),round(y_min,2)])
    plt.savefig("/lustre7/home/lustre3/satoshi/sa/ALL/paper_20210608/2024/kizami{0}.png".format(kizami),bbox_inches="tight")
    plt.close()
    print(kizami)
