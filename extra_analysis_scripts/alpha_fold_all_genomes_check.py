from alpha_fold_main import AlphaFoldWorkflow
from tqdm import tqdm

## ----- ALL DATA ------

data_h = open("../input_files/AF_non-human - All.tsv")
data = data_h.readlines()
data_h.close()
UNI_IDso = [i.split("\t") for i in data[1:]]

ref_list = []
UNI_IDs = []
for i in UNI_IDso:
    if i[19]:
        try:
            name = i[17]
            knots = i[19].split("{")[1].split("}")[0].split(")")
            for i1 in knots:
                knot_type = i1.split(": (")[0].strip("'").lstrip(", '")
                temp_range = i1.split(": (")[1].split(", ")
                UNI_IDs.append([i[1], [name, int(temp_range[0]), int(temp_range[1])]])
                ref_list.append([name, int(temp_range[0]), int(temp_range[1]), knot_type])

        except:
            continue

genomes = {}
for i in UNI_IDs:
    if i[0] in genomes:
        genomes[i[0]].append(i[1])
    else:
        genomes[i[0]] = [i[1]]

dict_names = list(genomes.keys())
dict_names.sort(key=lambda i: len(genomes[i]))

for i in tqdm(dict_names):
    print(i)
    result = AlphaFoldWorkflow(genomes[i], savefile=f"temp_files/{i}_alpha_fold.txt")

