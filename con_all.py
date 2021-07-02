
txt_path=["txt/aff","txt/eaf","txt/taf","kaitxt/affkai","kaitxt/eafkai","kaitxt/tafkai"]
for f in txt_path:
    if "kai" in f:
        n1 = 201
    else:
        n1 = 181
    date = [0.0] * 92
    for n2 in range(1, n1):
        for i in open("{0}_{1}.txt".format(f,n2)):
            fi = i.split()
            date = [float(j1) + float(j2) for (j1, j2) in zip(date, fi)]
    f2 = open("{0}_ALL.txt".format(f),"w")
    for ii in date:
        a =  str(ii) + "\n"
        f2.write(a)
    f2.close()
