
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("--pdb","--protein",type=str,help="path of pdb")
    return argparser.parse_args()

def load_pdb(path):
    num = 0
    CL = ["ASP","GLU"]
    NA = ["ARG","LYS","HIS"]
    for i in open(path,"r"):
        f = i.split()
        if f[2] == "CA":
            print(f[3])
            if f[3] in NA:
                num += 1
            elif f[3] in CL:
                num -= 1
        print(num)
    print(num)

def main():
    args = get_option()
    pdb = load_pdb(args.pdb)

if __name__ == '__main__':
    main()
