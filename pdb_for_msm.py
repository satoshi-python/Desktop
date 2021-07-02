
import tool
import itertools
import random
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("--trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("--pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("--prob","--probtxt",type=str,help="path of prob")
    argparser.add_argument("--txt","--txt_path",type=str,help="path of txt")
    return argparser.parse_args()

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

def load_txt(path):
    a = []
    num1 = 0
    for i in open(path,"r"):
        if i[:-1].isdecimal()==True:
            a.append(int(i))
        num1 += 1
    print(num1,len(a))
    return a

def main():
    args = get_option()
    pdb = str(args.pdb)
    prob = str(args.prob)
    num_prob = PROB(prob)
    a= load_txt(str(args.txt))
    #for i in range (100):
        #print([n for n, v in enumerate(a) if v == i])
    for i in range (100):
        tool.MD_to_pdb_chain(str(args.trr),pdb,[num_prob[1][num_prob[0].index(random.choice([n for n, v in enumerate(a) if v == i]))]])

if __name__ == '__main__':
    main()
