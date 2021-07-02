
import tool
from argparse import ArgumentParser
import random

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("--trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("--pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("--txt","--clutxt",type=str,help="path of prob")
    argparser.add_argument("--num","--cluster",type=str,help="number of cluster")
    argparser.add_argument("--TAKE","--number",type=str,help="number of cluster")

    return argparser.parse_args()


def load_prob():
    MED = []
    a = []
    for i in open("/home/biostr1/SOTU18/got/ALL/cluster/txt/sum_prob.txt"):
        f = i.split()
        if "t" not in f[0]:
            a.append(f)
        else:
            if len(a) > 1:
                MED.append(a)
            a = []
    return MED

def load_txt(path):
    ref = []
    for i in open(path,"r"):
        ref.append(int(i))
    return ref

def take_out(MED, txt):
    ref = []
    for i in MED:
        if int(i[1]) in txt:
            ref.append(int(i[0]))
    print(ref)
    return ref


def main():
    prob = load_prob()
    args = get_option()
    print(args)
    txt_list = load_txt(str(args.txt))
    a = take_out(prob[int(args.num)],txt_list)
    tool.TRR_TO_PDB(str(args.trr), str(args.pdb), random.sample(a,int(args.TAKE)))

if __name__ == '__main__':
    main()
