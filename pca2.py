
import satoshi_pca as SAS
import sys
a = sys.argv
# a[1]:result of pca1
# a[2]:prob_path_list
# a[3]:save_path_list
a.pop(0)
kai1 = a[0]
kai2 = [i for i in a if "prob" in i]
kai3 = [i for i in a if i not in kai2 and i not in kai1]
print(kai3)
a = []
a.append([0, 0])
a.append(kai1)
a.append(kai2)
a.append(kai3)


def PPP():
    global a
    kai = [0]
    num1 = 0
    for i in a[2]:
        for j in open(i, "r"):
            if float(j) != 0:
                num1 += 1
        kai.append(num1)
    return kai


def main():
    global a
    kai = []
    """
    for i in open(a[1], "r"):
        f = i.split()
        kai.append(f)
        if len(kai) == 10:
            break
    """
    with open(a[1]) as f:
        while True:
            a = f.readline().split()
            if not a:
                break
            kai.append([float(i.strip()) for i in a])
    print(len(kai))
    kai = SAS.pca(kai)
    kai = kai.tolist()
    num = PPP()
    for i in range(len(a[2])):
        f = open(a[3][i], "w")
        for j in range(int(num[i]), int(num[i+1])):
            if j < len(kai):
                f.write(str(kai[j][0]))
                f.write(" ")
                f.write(str(kai[j][1]))
                f.write("\n")
        f.close()


if __name__ == '__main__':
    main()
