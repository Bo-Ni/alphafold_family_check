from tqdm import tqdm
from matplotlib import pyplot as plt
from urllib import request


data_h = open("temp_files/savefilec.txt")
data = data_h.readlines()
data_h.close()
data_a = [i.rstrip().split(";") for i in data]
data_a = [[i[0].split("_"), [i1.split("_") for i1 in i[1].split(",")]] for i in data_a]
for i in range(len(data_a)):
    if not data_a[i][1][0][0]:
        data_a[i][1] = []


data_h = open("../input_files/Pfam-A.clans.tsv")
data = data_h.readlines()
data_h.close()
data_p = [i.rstrip().split("\t") for i in data]
data_p = [[i[0], i[3]] for i in data_p]



result = {}
all_families = []
for i0 in range(len(data_a)):
    temp = []
    for i in range(len(data_a[i0][1])):
        for i1 in range(len(data_p)):
            if data_p[i1][0] == data_a[i0][1][i][0]:
                temp.append(data_p[i1][1])
    for i2 in temp:
        if i2 not in all_families:
            all_families.append(i2)
    if ";".join(temp) not in result:
        result[";".join(temp)] = [data_a[i0][0][0]]
    else:
        result[";".join(temp)].append(data_a[i0][0][0])

list_names = list(result.keys())
list_names.sort(key=lambda i: len(result[i]), reverse=True)

# for i in list_names:
#     print(len(result[i]), i, result[i])

# print(len(result))
# print(len(all_families))

data_h = open("../input_files/AF_human - Arkusz1.tsv")
data = data_h.readlines()
data_h.close()
data_k = [i.rstrip().split("\t") for i in data[1:]]

data_k = [[i[0], [i1.strip() for i1 in i[3].split(",")], i[-1]] for i in data_k]


result_cut = {}
for i0 in result:
    new_temp = []
    for i in range(len(result[i0])):
        for i1 in data_k:
            if result[i0][i] == i1[0]:
                # print(i1[1][0], len(i1[1][0]))
                if len(i1[1][0]) == 0:
                    new_temp.append([result[i0][i], i1[2]])
                    # result[i0][i] = f"{result[i0][i]}*"  # Znakowanie
    if new_temp:
        result_cut[i0] = new_temp

# list_names = list(result.keys())
# list_names.sort(key=lambda i: len(result[i]), reverse=True)
# for i in list_names:
#     print(len(result[i]), i, result[i])
#
#
# list_names = list(result_cut.keys())
# list_names.sort(key=lambda i: len(result_cut[i]), reverse=True)
# for i in list_names:
#     print(len(result_cut[i]), i, result_cut[i])




# TODO - alternatywne klastrowanie
# TODO - 90% check


list_names = list(result_cut.keys())
list_names.sort(key=lambda i: len(result_cut[i]), reverse=True)
for i in list_names:
    if not i:
        print("Domains (Pfam families):", "NO FAMILY", " Number of hits: ", len(result_cut[i]))
    else:
        print("Domains (Pfam families):", i, " Number of hits: ", len(result_cut[i]))
    for i1 in result_cut[i]:
        print("\t" + " -- knot: ".join(i1))






