
import MDAnalysis
import sys
import itertools
import tool
from argparse import ArgumentParser
"""
a = sys.argv
a.pop(0)
kai1 = [i for i in a if ".trr" in i]
kai2 = [i for i in a if ".pdb" in i]
kai3 = [i for i in a if "prob" in i]
kai4 = [i for i in a if ".trr" not in i and ".pdb" not in i
        and ".txt" not in i and ".dat" not in i]
a = []
a.append([0, 0])
a.append(kai1)
a.append(kai2)
a.append(kai3)
a.append(kai4)
"""
# a = [i.split() for i in a]
# a[1]:trr_path_list
# a[2]:pdb_path_list
# a[3]:prob_path_list
# a[4]:data processing
"""
1:only CA coordinates
2:chain A and chain B | only CA
3:select residue coordinates like["3","C","N","O"]
4:chain A and chain B | selesct ATOM
"""
a=[]
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("-pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("-prob","--probtxt",type=str,help="path of prob")
    argparser.add_argument("-cal","--caluculation",type=str,help="way of data processing")
    return argparser.parse_args()

def PDB_cal1(num1):
    num_pdb = []
    for i in open(a[1][num1], "r"):
        f = i.split()
        if f[2] == "CA":
            num_pdb.append(int(f[1]))
    return num_pdb

def PDB_cal2(num1):
    num_pdb = []
    kai_a = []
    kai_b = []
    for i in open(a[1][num1], "r"):
        f = i.split()
        if f[4] == "A" and f[2] == "CA":
            kai_a.append(int(f[1]))
        if f[4] == "B" and f[2] == "CA":
            kai_b.append(int(f[1]))
    num_pdb.append(kai_a)
    num_pdb.append(kai_b)
    return num_pdb

def PDB_cal3(num1):
    num_pdb = []
    for i in open(a[1][num1], "r"):
        f = i.split()
        if f[2] in a[3]:
            num_pdb.append(int(f[1]))
    return num_pdb


def PDB_cal4(num1):
    num_pdb = []
    kai_a = []
    kai_b = []
    for i in open(a[1][num1], "r"):
        f = i.split()
        if f[4] == "A" and f[2] in a[3]:
            kai_a.append(int(f[1]))
        if f[4] == "B" and f[2] in a[3]:
            kai_b.append(int(f[1]))
    num_pdb.append(kai_a)
    num_pdb.append(kai_b)
    return num_pdb


def PROB_cal(num1):
    num_prob = []
    num2 = 0
    for i in open(a[2][num1], "r"):
        if float(i) != 0:
            num_prob.append(num2)
        num2 += 1
    return num_prob


def main():
    global a
    args = get_option()
    a.append(str(args.trajectory).split(","))
    a.append(str(args.protein).split(","))
    a.append(str(args.probtxt).split(","))
    a.append(str(args.caluculation).split(","))
    if len(a[0]) == len(a[1]) and len(a[1]) == len(a[2]):
        print("go")
        for i in a[0]:
            num1 = a[0].index(i)
            if len(a[3]) > 1:
                if int(a[3][0]) == 4:
                    num_pdb = PDB_cal4(num1)
                elif int(a[3][0]) == 3:
                    num_pdb = PDB_cal3(num1)
            else:
                if int(a[3][0]) == 1:
                    num_pdb = PDB_cal1(num1)
                elif int(a[3][0]) == 2:
                    num_pdb = PDB_cal2(num1)
            num_prob = PROB_cal(num1)
            u = MDAnalysis.Universe(i)
            frm = u.trajectory
            del u
            if int(a[3][0]) == 2 or int(a[3][0]) == 4:
                for i in safe_mem_distance(frm, num_prob, num_pdb):
                    k = i.replace(",", "")
                    k = k.replace("[", "")
                    k = k.replace("]", "")
                    print(k)
            else:
                for i in safe_mem_coordinates(frm, num_prob, num_pdb):
                    k = i.replace(",", "")
                    k = k.replace("[", "")
                    k = k.replace("]", "")
                    print(k)
        else:
            print("not match kind of trr or pdb, prob")


def safe_mem_distance(frm, num_prob, num_pdb):
    num2 = 0
    for frm_num in frm:
        try:
            kai = []
            if num2 in num_prob:
                for i in safe_safe_distance(num_pdb, frm_num):
                    kai.append(i)
                yield str(kai)
                # del x, y, z, FRM, kai
        except StopIteration:
            del kai
            break
        num2 += 1


def safe_safe_distance(num_pdb, FRM):
    for j1, j2 in itertools.product(num_pdb[0], num_pdb[1]):
        yield tool.contact(float(FRM[j1][0]),
                           float(FRM[j1][1]),
                           float(FRM[j1][2]),
                           float(FRM[j2][0]),
                           float(FRM[j2][1]),
                           float(FRM[j2][2]))


def safe_mem_coordinates(frm, num_prob, num_pdb):
    frm_itr = iter(frm)
    num2 = 0
    while True:
        try:
            kai = []
            FRM = next(frm_itr)
            if num2 in num_prob:
                for j1 in num_pdb:
                    for j2 in range(3):
                        kai.append(float(FRM[j1][j2]))
                yield str(kai)
                # del x, y, z, FRM, kai
        except StopIteration:
            del kai
            break
        num2 += 1


if __name__ == '__main__':
    main()
