
import satoshi_cluster
import MDAnalysis
import tool
import random
import itertools
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("--trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("--pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("--prob","--probtxt",type=str,help="path of prob")
    argparser.add_argument("--clu","--cluster",type=str,help="number of cluster")
    return argparser.parse_args()

def load_pdb(pdb_path):
    num_pdb = []
    A = []
    B = []
    for i in open(pdb_path,"r"):
        f = i.split()
        if f[2] == "CA" and f[4] == "A":
            A.append(int(f[1]))
        elif f[2] == "CA" and f[4] == "B":
            B.append(int(f[1]))
    num_pdb.append(A)
    num_pdb.append(B)
    return num_pdb

def load_prob(prob_path):
    num_prob = []
    num2 = 0
    for i in open(prob_path, "r"):
        if float(i) != 0:
            num_prob.append(num2)
        num2 += 1
    return num_prob

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

def PROB(path):
    num_prob = []
    num2 = 0
    num1 = 0
    a1 = []
    a2 = []
    for i in open(path, "r"):
        a2.append(num1)
        if float(i) != 0:
            a1.append(num2)
            num2 += 1
        else:
            a1.append("*")
        num1 += 1
    num_prob.append(a1)
    num_prob.append(a2)
    return num_prob

def main():
    args = get_option()
    pdb = str(args.pdb)
    prob = str(args.prob)
    clu = int(args.clu)
    num_pdb = load_pdb(pdb)
    num_prob = load_prob(prob)
    u = MDAnalysis.Universe(str(args.trr))
    frm = u.trajectory
    date = []
    print("ok")
    for i in safe_mem_distance(frm, num_prob, num_pdb):
        k = i.replace(",", "")
        k = k.replace("[", "")
        k = k.replace("]", "")
        date.append(k.split())
    print("calculation")
    a = satoshi_cluster.k_means(date,clu,random.randint(1, 10000))
    for i in a:
        print(i)
    frm_num = []
    num_prob = PROB(prob)
    for i in range (100):
        tool.MD_to_pdb_chain(str(args.trr),pdb,[num_prob[1][num_prob[0].index(random.choice([n for n, v in enumerate(a) if v == i]))]])

if __name__ == '__main__':
    main()
