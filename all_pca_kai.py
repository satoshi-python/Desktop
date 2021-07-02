
import MDAnalysis
import satoshi_pca as SAS
path = "/lustre7/home/lustre3/satoshi/MED"
TRR = ["/aff4/test_all.trr", "/eaf1/test_all.trr", "/taf7/test_all.trr",
       "/aff4_kai/run_all.trr", "/eaf1_kai/run_all.trr",
       "/taf7_kai/run_all.trr"]
PDB = ["/aff4/HEN.pdb", "/eaf1/HEN.pdb", "/taf7/HEN.pdb",
       "/aff4_kai/aff4kai.pdb", "/eaf1_kai/eaf1kai.pdb",
       "/taf7_kai/taf7kai.pdb"]
PROB = ["/aff4/prob.txt", "/eaf1/prob.txt", "/taf7/prob.txt",
        "/aff4_kai/prob.dat", "/eaf1_kai/prob.dat", "/taf7_kai/prob.dat"]


def PDB_cal(num1):
    num_pdb = []
    RESIDUE = ["N", "C"]
    for i in open(path + PDB[num1]):
        f = i.split()
        if f[2] in RESIDUE:
            num_pdb.append(int(f[1]))
    print(len(num_pdb))
    return num_pdb


def PROB_cal(num1):
    num_prob = []
    num2 = 0
    for i in open(path + PROB[num1], "r"):
        if float(i) != 0:
            num_prob.append(num2)
        num2 += 1
    return num_prob


def TRR_cal():
    kai_zahyou = []
    for trr in range(6):
        num_pdb = PDB_cal(trr)
        num_prob = PROB_cal(trr)
        u = MDAnalysis.Universe(path + TRR[trr])
        frm = u.trajectory
        frm_itr = iter(frm)
        del frm, u
        print(len(num_prob))
        """
        for i in num_prob:
            kai = []
            x = float(frm[i][0][0])
            y = float(frm[i][0][1])
            z = float(frm[i][0][2])
            for j in num_pdb:
                kai.append(str(float(frm[i][j][0]) - x))
                kai.append(str(float(frm[i][j][1]) - y))
                kai.append(str(float(frm[i][j][2]) - z))
            kai_zahyou.append(kai)
            print("kai", len(kai), " kai_zahyou", len(kai_zahyou), "/",
                  len(num_prob), " num_pdb", len(num_pdb))
        num2 = 0
        while True:
            try:
                kai = []
                FRM = next(frm_itr)
                if num2 in num_prob:
                    x = float(FRM[0][0])
                    y = float(FRM[0][1])
                    z = float(FRM[0][2])
                    for j in num_pdb:
                        kai.append(str(float(FRM[j][0]) - x))
                        kai.append(str(float(FRM[j][1]) - y))
                        kai.append(str(float(FRM[j][2]) - z))
                    kai_zahyou.append(kai)
                    print("kai", len(kai), " kai_zahyou", len(kai_zahyou), "/",
                          len(num_prob), " num_pdb", len(num_pdb))
                    del x, y, z, FRM
            except StopIteration:
                break
            num2 += 1
            """
        for i in safe_mem(frm_itr, num_prob, num_pdb):
            kai_zahyou.append(i)
            print("kai_zahyou", len(kai_zahyou), "/",
                  len(num_prob), " num_pdb", len(num_pdb))
        del frm_itr, num_prob, num_pdb
    return kai_zahyou


def safe_mem(frm_itr, num_prob, num_pdb):
    num2 = 0
    while True:
        try:
            kai = []
            FRM = next(frm_itr)
            if num2 in num_prob:
                x = float(FRM[0][0])
                y = float(FRM[0][1])
                z = float(FRM[0][2])
                for j in num_pdb:
                    kai.append(str(float(FRM[j][0]) - x))
                    kai.append(str(float(FRM[j][1]) - y))
                    kai.append(str(float(FRM[j][2]) - z))
                yield kai
                # del x, y, z, FRM, kai
        except StopIteration:
            del kai
            break
        num2 += 1


def PPP():
    kai = [0]
    num1 = 0
    for i in range(6):
        for i in open(path + PROB[i], "r"):
            if float(i) != 0:
                num1 += 1
        kai.append(num1)
    return kai


if __name__ == '__main__':
    kai = SAS.pca(TRR_cal())
    kai = kai.tolist()
    path1 = "/lustre7/home/lustre3/satoshi/ALL_PCA/txt/"
    ttt = "_20201207_2.txt"
    ligands = ["aff", "eaf", "taf", "affkai", "eafkai", "tafkai"]
    num = PPP()
    for i in range(6):
        f = open(path1+ligands[i]+ttt, "w")
        for j in range(int(num[i]), int(num[i+1])):
            f.write(str(kai[j][0]))
            f.write(" ")
            f.write(str(kai[j][1]))
            f.write("\n")
        f.close()
