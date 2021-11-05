from tqdm import tqdm
from matplotlib import pyplot as plt
from urllib import request

def draw_blobs(data):
    plot_title = data[0][0]
    main_range = [int(data[0][1]), int(data[0][2])]
    full_range = [0, data[0][-1]]
    mini_ranges = [[i[0], int(i[2]), int(i[3])] for i in data[1]]

    plt.figure(figsize=((8, 4)))
    plt.hlines(0.5, full_range[0], full_range[1], color="k", linewidth=3)
    plt.hlines(0.5, main_range[0], main_range[1], color="g", linewidth=10)
    plt.text(x=((full_range[1] - full_range[0]) // 2), y=0.75, s=plot_title, fontsize=15)
    plt.text(main_range[0], 0.6, str(main_range[0]))
    plt.text(main_range[1], 0.6, str(main_range[1]))
    for i in mini_ranges:
        plt.hlines(0.4, i[1], i[2], linewidth=5, color="b")
        plt.plot([i[1], i[2]], [0.4, 0.4], "o", c="r", markersize=4)
        plt.text(i[1] + ((i[2] - i[1]) // 2), 0.3, i[0], rotation_mode="anchor", rotation="-45", fontsize=11)
        plt.text(i[1], 0.3, str(i[1]), rotation_mode="anchor", rotation="-45")
        # plt.text(i[2], 0.3, str(i[1]), rotation_mode="anchor", rotation="-45")
    plt.box(False)
    plt.yticks([])
    plt.tight_layout()
    plt.ylim(-1, 1)
    plt.savefig(f"plots/{plot_title}.png")
    plt.show()


data_h = open("temp_files/savefilec.txt")
data = data_h.readlines()
data_h.close()
data_a = [i.rstrip().split(";") for i in data]
data_a = [[i[0].split("_"), [i1.split("_") for i1 in i[1].split(",")]] for i in data_a]
for i in range(len(data_a)):
    if not data_a[i][1][0][0]:
        data_a[i][1] = []


data_h = open("input_files/Pfam-A.clans.tsv")
data = data_h.readlines()
data_h.close()
data_p = [i.rstrip().split("\t") for i in data]
data_p = [[i[0], i[3]] for i in data_p]

# for i0 in tqdm(range(len(data_a))):
#     url = f"https://www.uniprot.org/uniprot/{data_a[i0][0][0]}.txt"
#     temp_url = request.urlopen(url)
#     temp_url = temp_url.readlines()
#     temp_url = int(temp_url[0].decode("utf-8").rstrip().split(";")[-1].strip().rstrip(" AA."))
#     data_a[i0][0].append(temp_url)
#     for i in range(len(data_a[i0][1])):
#         for i1 in range(len(data_p)):
#             if data_a[i0][1][i][0] == data_p[i1][0]:
#                 t = [data_p[i1][1]]
#                 t.extend(data_a[i0][1][i])
#                 data_a[i0][1][i] = t
#     draw_blobs(data_a[i0])

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

data_h = open("input_files/AF_human - Arkusz1.tsv")
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






