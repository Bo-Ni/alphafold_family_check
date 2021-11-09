from urllib import request
from tqdm import tqdm
from matplotlib import pyplot as plt


def draw_blobs(data):
    plot_title = data[0][0]
    main_range = [int(data[0][1]), int(data[0][2])]
    full_range = [0, data[0][-1]]
    try:
        print("OK", data)
        mini_ranges = [[i[0], int(i[2]), int(i[3])] for i in data[1]]
    except:
        print(data)

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


data_h = open("temp_files/savefile_human.txt")
data = data_h.readlines()
data_h.close()
data_a = [i.rstrip().split(";") for i in data]
data_a = [[i[0].split("_"), [i1.split("_") for i1 in i[1].split(",")]] for i in data_a]
for i in range(len(data_a)):
    if not data_a[i][1][0][0]:
        data_a[i][1] = []

data_h = open("./input_files/Pfam-A.clans.tsv")
data = data_h.readlines()
data_h.close()
data_p = [i.rstrip().split("\t") for i in data]
data_p = [[i[0], i[3]] for i in data_p]

for i0 in tqdm(range(len(data_a))):
    url = f"https://www.uniprot.org/uniprot/{data_a[i0][0][0]}.txt"
    temp_url = request.urlopen(url)
    temp_url = temp_url.readlines()
    temp_url = int(temp_url[0].decode("utf-8").rstrip().split(";")[-1].strip().rstrip(" AA."))
    data_a[i0][0].append(temp_url)
    for i in range(len(data_a[i0][1])):
        for i1 in range(len(data_p)):
            if data_a[i0][1][i][0] == data_p[i1][0]:
                t = [data_p[i1][1]]
                t.extend(data_a[i0][1][i])
                data_a[i0][1][i] = t
    draw_blobs(data_a[i0])