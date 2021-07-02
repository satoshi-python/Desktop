
import MDAnalysis

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("-pdb","--protein",type=str,help="path of pdb")
    return argparser.parse_args()

def main():
    args = get_option()
    u = MDAnalysis.Universe(str(args.protein),str(args.trajectory))
    frm = u.trajectory
    print("number of frame:",len(frm))
    print("number of atom:", len(frm[0]))
    print("cordinate of atom[first atom]:",frm[0][0])

if __name__ == '__main__':
    main()
