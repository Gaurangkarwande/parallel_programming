import os
from collections import defaultdict

def main():
    data_file = './data/Graph_directed.txt'
    num_edges = 0

    with open(data_file) as f:
        lines = f.readlines()
        for line in lines:
            source, dest, weight = line.split()
            source = int(source)
            dest = int(dest)
            weight = float(weight)
            A_list[source].append((dest, weight))
            num_edges += 1

    V = []
    I = []
    prev = 0
    for i in range(min(A_list.keys()), max(A_list.keys())+1):
        V.append(i)
        I.append(prev)
        prev += len(A_list[i])
    I.append(prev)

    E = []
    W = []
    for vertex in V:
        for edge in A_list[vertex]:
            E.append(edge[0])
            W.append(edge[1])

    v_file = './data/Graph_directed_V.txt'
    i_file = './data/Graph_directed_I.txt'
    e_file = './data/Graph_directed_E.txt'
    w_file = './data/Graph_directed_W.txt'

    with open(v_file, "w") as f:
        for vertex in V:
            f.write(str(vertex) + "\n")

    with open(i_file, "w") as f:
        for i in I:
            f.write(str(i) + "\n")

    with open(e_file, "w") as f:
        for e in E:
            f.write(str(e) + "\n")

    with open(w_file, "w") as f:
        for w in W:
            f.write(str(w) + "\n")

if __name__ = "__main__":
    main()