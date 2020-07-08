import sys

edgef = sys.argv[1]
nscaf = int(sys.argv[2])
threshold = float(sys.argv[3])

edge = []
with open(edgef) as ifs:
    for l in ifs.readlines():
        try:
            edge.append((int(l.split()[0]), int(l.split()[1]), float(l.split()[2])))
        except:
            continue

edge = sorted(edge, key=lambda x: x[0])

def find(node):
    for e in edge:
        if e[0] == node or e[1] == node:
            if e[0] % 2 == 0 and e[1] - e[0] == 1:
                continue
            else:
                return e

def findpath(node):
    path = [node]
    n1, n2, value = find(node)
    total_edge = 0
    valid_edge = 0
    # path.append(value)
    if n1 == node:
        n = n2
    else:
        n = n1
    while n not in path and value > threshold:
        total_edge += 1
        if abs(n1 - n2) == 1:
            valid_edge += 1
        path.append(value)
        path.append(n)
        if n % 2 == 0:
            ll = n + 1
        else:
            ll = n - 1
        path.append(ll)
        n1, n2, value = find(ll)
        if n1 == ll:
            n = n2
        else:
            n = n1
        # path.append(value)
    if node % 2 == 0:
        n = node + 1
    else:
        n = node - 1
    path = [n] + path
    n1, n2, value = find(n)
    if n1 == n:
        n = n2
    else:
        n = n1
    while n not in path and value > threshold:
        total_edge += 1
        if abs(n1 - n2) == 1:
            valid_edge += 1
        # path = [n] + path
        path = [n, value] + path
        if n % 2 == 0:
            ll = n + 1
        else:
            ll = n - 1
        path = [ll] + path
        n1, n2, value = find(ll)
        if n1 == ll:
            n = n2
        else:
            n = n1
    return path, total_edge, valid_edge

def findnode(i, graph):
    for e in graph:
        if e[0] == i or e[1] == i:
            if e[0] % 2 == 0 and e[1] - e[0] == 1:
                continue
            else:
                if e[0] == i:
                    return e[1]
                else:
                    return e[0]

def switch_edge(i, exp, graph):
    orig1 = -1
    orig2 = -1
    for idx, e in enumerate(graph):
        if e[0] == i or e[1] == i:
            orig1 = idx
        if e[0] == exp or e[1] == exp:
            orig2 = idx
    tmp1, tmp2 = graph[orig1]
    tmp3, tmp4 = graph[orig2]
    # print "i, exp: ", i, exp
    # print "orig1, orig2: ", orig1, orig2
    # print "edge1, edge2: ", graph[orig1], graph[orig2]
    # if tmp2 == i or tmp4 == exp:
    #     raise IndexError("swith index error")
    if tmp1 == i and tmp3 == exp:
        graph[orig1] = [tmp1, tmp3]
        graph[orig2] = [tmp2, tmp4]
    if tmp1 == i and tmp4 == exp:
        graph[orig1] = [tmp1, tmp4]
        graph[orig2] = [tmp2, tmp3]
    if tmp2 == i and tmp3 == exp:
        graph[orig1] = [tmp2, tmp3]
        graph[orig2] = [tmp1, tmp4]
    if tmp2 == i and tmp4 == exp:
        graph[orig1] = [tmp2, tmp4]
        graph[orig2] = [tmp1, tmp3]

def valid(graph):
    ret = True
    for i in range(1, 2*nscaf, 2):
        expect = i + 1
        if (i + 1) % 40 == 0:
            expect = i + 1 - 40
        if (findnode(i, graph)) != expect:
            ret = False
            break
    return ret

def get_switch(graph):
    switch = 0
    while not valid(graph):
        for i in range(1, 2*nscaf, 2):
            node = findnode(i, graph)
            expect = i + 1
            if (i + 1) % 40 == 0:
                expect = i + 1 - 40
            if node != expect:
                switch_edge(i, expect, graph)
                switch += 1
                break
    return switch

checked = []
paths = []
total_edge = 0
valid_edge = 0
for i in range(2*nscaf):
    if i not in checked:
        path, tedge, vedge = findpath(i)
        total_edge += tedge
        valid_edge += vedge
        paths.append(path)
        checked = path + checked

G = [e[:2] for e in edge]
print(path)
switch = get_switch(G)

def print_path(paths):
    for no, path in enumerate(paths):
        outstring = ""
        for idx, name in enumerate(path):
            if idx % 3 == 0:
                outstring += str(name)
            if idx % 3 == 1:
                outstring += ", " + str(name)
            if idx % 3 == 2:
                if idx != len(path) - 1:
                    outstring += ", " + str(name) + " -> "
                else:
                    outstring += ", " + str(name)
        print "------------------------------ path %d --------------------------------" % no
        print outstring
        print "-----------------------------------------------------------------------"

print "path number: ", len(paths)
print_path(paths)
print "switch (edit distance): ", switch
print "switch error rate: {}%".format(float(switch) * 100 / nscaf)
print "total edge count: ", total_edge
print "valid edge count: ", valid_edge
error_rate = float(total_edge - valid_edge) / float(total_edge) if total_edge > 0 else 0
print "error rate: {}%".format(100 * error_rate)

