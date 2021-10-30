from tqdm import tqdm
from urllib import request
# data_h = open("AF_human - Arkusz1.tsv")
# data = data_h.readlines()
# data_h.close()
# UNI_IDs = [i.split("\t")[0] for i in data[1:]]

data_h = open("Knots on Alphafold proteins - Sulkowska LabTable.csv")
data = data_h.readlines()
data_h.close()
UNI_IDs = [i.split(",")[1].strip('"').rstrip("F1").rstrip("-") for i in data[1:]]
print(len(UNI_IDs))


links = [f"https://www.uniprot.org/uniprot/{i}.txt" for i in UNI_IDs]
result = {"No_family": []}
result_separated = {"No_family": []}
for i0 in tqdm(range(len(links))):
    try:
        temp = request.urlopen(links[i0]).readlines()
    except:
        print("ERROR: ", UNI_IDs[i0])
        continue
    temp = [i.decode("utf-8") for i in temp]
    temp = [i for i in temp if "Pfam" in i]
    if len(temp) == 0:
        result["No_family"].append(UNI_IDs[i0])
        result_separated["No_family"].append(UNI_IDs[i0])
    else:
        temp = [i.split(";")[1].strip() for i in temp]
        for i1 in temp:
            if i1 not in result_separated:
                result_separated[i1] = [UNI_IDs[i0]]
            else:
                result_separated[i1].append(UNI_IDs[i0])
        if "_".join(temp) in result:
            result["_".join(temp)].append(UNI_IDs[i0])
        else:
            result["_".join(temp)] = [UNI_IDs[i0]]


# families_list = list(result_separated.keys())
# families_list.sort(key=lambda i: len(result_separated[i]), reverse=True)
# for i in families_list:
#     print(len(result_separated[i]), i, result_separated[i])

families_list = list(result.keys())
families_list.sort(key=lambda i: len(result[i]), reverse=True)
for i in families_list:
    print(len(result[i]), i, result[i])
