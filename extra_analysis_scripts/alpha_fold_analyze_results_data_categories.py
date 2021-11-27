from urllib import request
from tqdm import tqdm
from matplotlib import pyplot as plt
from alpha_fold_support_scripts import DecodeSavefile
import os


def SplitResultToCategories(data, output="output_info.txt"):
    """
    :param data: input data as list returned from alpha_fold_domain_check.DecodeSavefile function
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


if __name__ == "__main__":
    for i in os.listdir("../temp_files/organism_data"):
        parsed_data = DecodeSavefile(f"../temp_files/organism_data/{i}")
        SplitResultToCategories(parsed_data, output=f"../info/organism_summaries/{i.split('_')[0]}_data_summary.txt")

