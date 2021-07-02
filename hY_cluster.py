import tool
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("-pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("-prob","--probtxt",type=str,help="path of prob")
    argparser.add_argument("-name","--name_kei",type=str,help="name of file")
    return argparser.parse_args()

def make_list(path):
    ll = []
    for i in open(path,"r"):
        ll.append(int(i))
    return ll

def main():
    args = get_option()
    L = make_list(str(args.name_kei))
    J,K = tool.HY_CONTACT(str(args.protein),str(args.probtxt),str(args.trajectory),L)
    tool.heat_map("hy.png","title","MED26","ligand",J)
    tool.heat_map("hy_non_prob.png","title","MED26","ligand",K)
    tool.yoko_tate(J,"hy_YT.pdf")
    H = tool.MD_LIST_CONTACT(str(args.protein),str(args.trajectory),str(args.probtxt),str(args.name_kei))
    print(len(H))
    print(len(H[0]))
    tool.yoko_tate(H,"con_YT.pdf")

if __name__ == '__main__':
    main()
