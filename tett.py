
test = [str(i)for i in range(1, 100)]
kai = [0]
for i in range(10):
    kai1 = []
    for j in range(len(test)):
        if i == int(test[j][0:1]):
            kai1.append(test[j])
    kai = kai + kai1
test = []
for i in kai:
    a = str(str(i).zfill(4)) + ".trr"
    test.append(a)
print(test)
