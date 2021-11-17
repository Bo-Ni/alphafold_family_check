from urllib import request
from tqdm import tqdm
from matplotlib import pyplot as plt
import os


def DecodeSavefile(input_file):
    file_h = open(input_file)
    file = file_h.readlines()
    file_h.close()
    file = [i.rstrip().split(";") for i in file]

    file = [[
        i[0].split("_"),
        i[1],
        i[2],
        [i1.split("_") for i1 in i[3].split(",")],
        i[4],
        i[5],
        [i1 for i1 in i[6].split(",")]
    ] for i in file]

    return file


# def DrawPlots():  # TODO MIGHT NEED AN UPDATE
#     # data_h = open("temp_files/savefile_human.txt")
#     # data = data_h.readlines()
#     # data_h.close()
#     # data_a = [i.rstrip().split(";") for i in data]
#     # data_a = [[i[0].split("_"), [i1.split("_") for i1 in i[1].split(",")]] for i in data_a]
#     # for i in range(len(data_a)):
#     #     if not data_a[i][1][0][0]:
#     #         data_a[i][1] = []
#     #
#     # data_h = open("./input_files/Pfam-A.clans.tsv")
#     # data = data_h.readlines()
#     # data_h.close()
#     # data_p = [i.rstrip().split("\t") for i in data]
#     # data_p = [[i[0], i[3]] for i in data_p]
#     #
#     # for i0 in tqdm(range(len(data_a))):
#     #     url = f"https://www.uniprot.org/uniprot/{data_a[i0][0][0]}.txt"
#     #     temp_url = request.urlopen(url)
#     #     temp_url = temp_url.readlines()
#     #     temp_url = int(temp_url[0].decode("utf-8").rstrip().split(";")[-1].strip().rstrip(" AA."))
#     #     data_a[i0][0].append(temp_url)
#     #     for i in range(len(data_a[i0][1])):
#     #         for i1 in range(len(data_p)):
#     #             if data_a[i0][1][i][0] == data_p[i1][0]:
#     #                 t = [data_p[i1][1]]
#     #                 t.extend(data_a[i0][1][i])
#     #                 data_a[i0][1][i] = t
#     #     draw_blobs(data_a[i0])
#
#     def DrawBlobs(data):
#         plot_title = data[0][0]
#         main_range = [int(data[0][1]), int(data[0][2])]
#         full_range = [0, data[0][-1]]
#         try:
#             print("OK", data)
#             mini_ranges = [[i[0], int(i[2]), int(i[3])] for i in data[1]]
#         except:
#             print(data)
#
#         plt.figure(figsize=((8, 4)))
#         plt.hlines(0.5, full_range[0], full_range[1], color="k", linewidth=3)
#         plt.hlines(0.5, main_range[0], main_range[1], color="g", linewidth=10)
#         plt.text(x=((full_range[1] - full_range[0]) // 2), y=0.75, s=plot_title, fontsize=15)
#         plt.text(main_range[0], 0.6, str(main_range[0]))
#         plt.text(main_range[1], 0.6, str(main_range[1]))
#         for i in mini_ranges:
#             plt.hlines(0.4, i[1], i[2], linewidth=5, color="b")
#             plt.plot([i[1], i[2]], [0.4, 0.4], "o", c="r", markersize=4)
#             plt.text(i[1] + ((i[2] - i[1]) // 2), 0.3, i[0], rotation_mode="anchor", rotation="-45", fontsize=11)
#             plt.text(i[1], 0.3, str(i[1]), rotation_mode="anchor", rotation="-45")
#             # plt.text(i[2], 0.3, str(i[1]), rotation_mode="anchor", rotation="-45")
#         plt.box(False)
#         plt.yticks([])
#         plt.tight_layout()
#         plt.ylim(-1, 1)
#         plt.savefig(f"plots/{plot_title}.png")
#         plt.show()


def SplitResultToCategories(data, output="output_info.txt"):
    """
    :param data: input data as list returned from DecodeSavefile function
    :param output: Path to output filename, where results will be saved
    """
    category1 = []  # Known structure for protein, known domains
    category2 = []  # Known structure for homologs, known domains
    category3 = []  # No structures for homologs, known domains for proteins of homologs
    category4 = []  # Known structure for protein or homologs, domain unknown for protein or known for homologs
    category5 = []  # No structures for homologs, no domains for homologs

    for i in tqdm(range(len(data))):
        if data[i][1] == "prot" and data[i][4] == "prot":
            category1.append(data[i])
        elif data[i][1] == "prot" and data[i][4] != "nodata" and data[i][6][0] != '':
            category2.append(data[i])
        elif data[i][3][0][0] != '' and data[i][6][0] == '':
            category3.append(data[i])
        elif data[i][1] != "prot" and data[i][6][0] != '':
            category4.append(data[i])
        else:
            category5.append(data[i])

    category1.sort(key=lambda i: int(i[5]), reverse=True)
    category2.sort(key=lambda i: len(i[6]), reverse=True)

    category3.sort(key=lambda i: int(i[2]), reverse=True)
    category3.sort(key=lambda i: i[1] if i[1] != "prot" else "")

    category4.sort(key=lambda i: len(i[6]), reverse=True)
    category4.sort(key=lambda i: i[4] if i[4] != "prot" else "")

    category5.sort(key=lambda i: i[0][0])


    for i in category5:
        print(i)
    print("-------")
            
    newfile = open(output, "w")

    newfile.write("Category 1: Known structure for protein, known domains\n")
    for i in category1:
        newfile.write(" ".join(["\t", "-".join(i[0]), "\n\t\t -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Structures:", ", ".join(i[6]), "\n"]))

    newfile.write("Category 2: Known structure for homologs, known domains\n")
    for i in category2:
        newfile.write(" ".join(["\t", "-".join(i[0]), "\n\t\t -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Homologs -- {len([i1 for i1 in i[6] if i1])} Structures:" , ", ".join(i[6]), "\n"]))

    newfile.write("Category 3: No structures for homologs, known domains for proteins of homologs\n")
    for i in category3:
        if i[1] == "prot":
            newfile.write(" ".join(["\t", "-".join(i[0]), "\n\t\t -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Homologs -- {len([i1 for i1 in i[6] if i1])} Structures:" , ", ".join(i[6]), "\n"]))
        else:
            newfile.write(" ".join(["\t", "-".join(i[0]), f"\n\t\t -- {i[2]} Homologs -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Homologs -- {len([i1 for i1 in i[6] if i1])} Structures:" , ", ".join(i[6]), "\n"]))

    newfile.write("Category 4: Known structure for protein or homologs, domain unknown for protein or known for homologs\n")
    for i in category4:
        if i[4] == "prot":
            newfile.write(" ".join(["\t", "-".join(i[0]), f"\n\t\t -- {i[2]} Homologs -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Structures:", ", ".join(i[6]), "\n"]))
        else:
            newfile.write(" ".join(["\t", "-".join(i[0]), f"\n\t\t -- {i[2]} Homologs -- Domains:", ", ".join(["-".join(i1) for i1 in i[3]]), f"\n\t\t -- {i[5]} Homologs -- {len([i1 for i1 in i[6] if i1])} Structures:", ", ".join(i[6]), "\n"]))

    newfile.write("Category 5: No structures for homologs, no domains for homologs\n")
    for i in category5:
        newfile.write(" ".join(["\t", "-".join(i[0]), "\n"]))


def decision_func(data):
    all = data[1]
    if data[2][0][1] / all > 0.6:
        data[2] = data[2][0]

    return data


if __name__ == "__main__":
    for i in os.listdir("temp_files"):
        parsed_data = DecodeSavefile(f"temp_files/{i}")
        SplitResultToCategories(parsed_data, output=f"info/{i.split('_')[0]}_data_summary.txt")

