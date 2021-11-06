from alpha_fold_domain_check import alpha_fold_knot_family_check


## ----- HUMAN -----

# data_h = open("Knots on Alphafold proteins - Sulkowska LabTable.csv")
# data = data_h.readlines()
# data_h.close()
# UNI_IDs = [[i.split(",")[1].strip('"').rstrip("F1").rstrip("-"),
#             int(i.split(",")[5].strip('"').split("-")[0]),
#             int(i.split(",")[5].strip('"').split("-")[1])] for i in data[1:]]
#
# result = alpha_fold_knot_family_check(UNI_IDs, parsed_tsv="uni_parsed_data_c.tsv", savefile="savefile.txt")
#
# for i in result:
#     print(i)


## ----- ALL DATA ------

data_h = open("input_files/AF_non-human - All (1).tsv")
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
                UNI_IDs.append([name, int(temp_range[0]), int(temp_range[1])])
                ref_list.append([name, int(temp_range[0]), int(temp_range[1]), knot_type])

        except:
            continue

for i in range(len(UNI_IDs)):
    print(UNI_IDs[i], ref_list[i])

result = alpha_fold_knot_family_check(UNI_IDs, savefile="savefile.txt")

for i in result:
    print(i)
