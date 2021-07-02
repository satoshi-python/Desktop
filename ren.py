
import tool
import sys

a = sys.argv
path = "/lustre7/home/lustre3/satoshi/ALL_PCA/txt"
txt_path = [path+"/aff_20201207_2.txt", path+"/eaf_20201207_2.txt",
            path+"/taf_20201207_2.txt", path+"/affkai_20201207_2.txt",
            path+"/eafkai_20201207_2.txt", path+"/tafkai_20201207_2.txt"]
path = "/lustre7/home/lustre3/satoshi/MED"
prob_path = [path+"/aff4/prob.txt", path+"/eaf1/prob.txt",
             path+"/taf7/prob.txt", path+"/aff4_kai/prob.dat",
             path+"/eaf1_kai/prob.dat", path+"/taf7_kai/prob.dat"]

trr_path = [path+"/aff4/test_all.trr", path+"/eaf1/test_all.trr",
            path+"/taf7/test_all.trr", path+"/aff4_kai/run_all.trr",
            path+"/eaf1_kai/run_all.trr", path+"/taf7_kai/run_all.trr"]
pdb_path = [path+"/aff4/HEN.pdb", path+"/eaf1/HEN.pdb", path+"/taf7/HEN.pdb",
            path+"/aff4_kai/aff4kai.pdb", path+"/eaf1_kai/eaf1kai.pdb",
            path+"/taf7_kai/taf7kai.pdb"]
kizami = int(a[1])
print(kizami)
tool.FREE_FROM_PCA_EACH(txt_path, prob_path, kizami, trr_path, pdb_path)
for i in range(6):
    num_prob = 0
    num_trr = 0
    for prob in open(prob_path[i], "r"):
        if float(prob) != 0:
            num_prob += 1
    for trr in open(txt_path[i], "r"):
        num_trr += 1
    print(num_prob, num_trr)
