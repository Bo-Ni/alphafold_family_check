import os
from alpha_fold_analyze_results import DecodeSavefile


def get_human_knot_info():
    data_h = open("../input_files/Knots on Alphafold proteins - Sulkowska LabTable.csv")
    data = data_h.readlines()
    data_h.close()
    knot_data = [[i.split(",")[1].strip('"').rstrip("F1").rstrip("-"),
                i.split(",")[3].strip('"')] for i in data[1:]]
    return knot_data


def get_non_human_knot_info():
    data_h = open("../input_files/AF_non-human - All.tsv")
    data = data_h.readlines()
    data_h.close()
    UNI_IDso = [i.split("\t") for i in data[1:]]

    ref_list = []
    for i in UNI_IDso:
        if i[19]:
            try:
                name = i[17]
                knots = i[19].split("{")[1].split("}")[0].split(")")
                for i1 in knots:
                    knot_type = i1.split(": (")[0].strip("'").lstrip(", '")
                    ref_list.append([name, knot_type])
            except:
                continue
    return ref_list


knot_list = get_human_knot_info()
knot_list.extend(get_non_human_knot_info())
knot_list = [[i[0], i[1].replace("_", "")] for i in knot_list]

result_dict = {}
organism_files = os.listdir("../temp_files")
for organism in organism_files:
    organism_data = DecodeSavefile(f"../temp_files/{organism}")
    for single_data in organism_data:
        knot_info = [i0[1] for i0 in knot_list if i0[0] == single_data[0][0]]
        knot_info = ", ".join([i for i in knot_info if i])
        if single_data[1] == "prot":  # ----- Protein data
            family_id = "--".join([i0[0] for i0 in single_data[3]])
            if family_id in result_dict:
                result_dict[family_id].append([organism.split('_')[0], '_'.join(single_data[0]), knot_info])
            else:
                result_dict[family_id] = [[organism.split('_')[0], '_'.join(single_data[0]), knot_info]]
        elif "hom" in single_data[1]:  # ----- Homolog data
            pass
        elif single_data[1] == "nodata":  # ----- Other
            if "No data" in result_dict:
                result_dict["No data"].append([organism.split('_')[0], '_'.join(single_data[0]), knot_info])
            else:
                result_dict["No data"] = [[organism.split('_')[0], '_'.join(single_data[0]), knot_info]]

family_names = list(result_dict.keys())
family_names.sort(key=lambda i: len(result_dict[i]), reverse=True)


newfile = open("../info/summary2.txt", "w")
for i in family_names:
    temp_org_check = []
    temp_top_check = []
    for i1 in result_dict[i]:
        if i1[0] not in temp_org_check:
            temp_org_check.append(i1[0])
        if i1[2] not in temp_top_check:
            temp_top_check.append(i1[2])
    newfile.write(f"Group -- N. of hits: {len(result_dict[i])} -- N. of organisms: {len(temp_org_check)} -- N. of topologies: {len(temp_top_check)} -- Family: " + i + "\n")
    for i1 in result_dict[i]:
        newfile.write("\t" + "---".join(i1) + "\n")
newfile.close()
