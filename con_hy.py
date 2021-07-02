
import tool
import threading
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr","--trajectory",type=str,help="path of trr")
    argparser.add_argument("-pdb","--protein",type=str,help="path of pdb")
    argparser.add_argument("-prob","--probtxt",type=str,help="path of pdb")
    argparser.add_argument("-name","--name_kei",type=str,help="name of file")
    return argparser.parse_args()


class HY_Calculation(threading.Thread):
    def __init__(self,thread_name, niji, doner, ac, doner_cal, trr, prob):
        self.name = thread_name
        self.ac = ac
        self.do = doner
        self.niji = niji
        self.do_cal = doner_cal
        self.syoki = [0.0] * 92
        self.prob = prob
        self.trr = trr
        del u, niji, ac, doner, LIST

    def ca2(self):
        for self.num1, self.num2 in itertools.product(self.num_a,
                                                      self.num_b):
            # define doner or accepter
            if self.num1 in self.ac[0] and self.num2 in self.ac[1]:
                # distace
                if contact_list(self.frm[self.num1],self.frm[self.num2]) <= 3.3:
                    # chain_a is donner
                    if self.num1 in self.do[0]:
                        # angle
                        num_angle = []
                        for ii in range(1, int(self.do_cal[0][self.num_a.index(self.num1)]) +1):
                            num_angle.append(angle(self.frm[self.num1][0], self.frm[self.num1][1],self.frm[self.num1][2],self.frm[self.num2][0],self.frm[self.num2][1],self.frm[self.num2][2],self.frm[self.num1 + ii][0],self.frm[self.num1 + ii][1],self.frm[self.num1 + ii][2]))
                        if len([kkknn for kkknn in num_angle if kkknn >= 120]) > 0:
                            self.syoki[self.niji[0].index(self.num_a)] += self.prob[self.num_prob]
                            # print("ok  ", self.num1, self.num2)
                            break
                    # chain_b is donner
                    if self.num2 in self.do[1]:
                        # anngle
                        num_angle = []
                        for ii in range(1, int(self.do_cal[1][self.num_b.index(self.num2)]) +1):
                            num_angle.append(angle(self.frm[self.num1][0],self.frm[self.num1][1],self.frm[self.num1][2],self.frm[self.num2][0],self.frm[self.num2][1],self.frm[self.num2][2],self.frm[self.num2 + ii][0],self.frm[self.num2 + ii][1],self.frm[self.num2 + ii][2]))
                        if len([kkknn for kkknn in num_angle if kkknn >= 120]) > 0:
                            self.syoki[self.niji[0].index(self.num_a)] += self.prob[self.num_prob]
                            # print("ok  ", self.num1, self.num2)
                            break
                elif contact_list(self.frm[self.num1][0], self.frm[self.num1][1], self.frm[self.num1][2],self.frm[self.num2][0],self.frm[self.num2][1],self.frm[self.num2][2]) > 15:
                    break


    def ca1(self):
        for self.num_a, self.num_b in itertools.product(self.niji[0],
                                                        self.niji[1]):
            yield

    def calculation(self):
        for i in range(self.frm):
            if self.prob[i] != 0.0:
                for i in self.ca1():
                        pass

def main():
    args = get_option()
    PDB = tool.sepa(str(args.protein))
    # print(PDB[0])
    niji = tool.NOT_H(str(args.protein))
    DONER = []
    ACCEPTER = []
    DONER_cal = []
    for i in range(2):
        DONER.append(donner(PDB[i]))
        DONER_cal.append(donner_cal(PDB[i]))
        ACCEPTER.append(accepter(PDB[i]))
    # print(len(DONER[1]))
    # print(len(DONER_cal[1]))
    u = MDAnalysis.Universe(str(args.trajectory))
    frm = u.trajectory
    PROB = []
    for i in open(str(args.probtxt)):
        PROB.append(float(i))
    for i in range(1,(len(frm)/1000) + 1):
        cal = HY_Calculation(i,niji, DONER, ACCEPTER, DONER_cal, trr[(i-1) * 1000:(i-1) * 1000 + 1000], prob[(i-1) * 1000:(i-1) * 1000 + 1000])
    cal.calculation()
    return cal.syoki,cal.syoki2


if __name__ == '__main__':
    main()
