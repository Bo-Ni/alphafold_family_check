from tqdm import tqdm
from urllib import request

def parse_data_from_Uniprot(data, filename):
    n = 50
    UNI_IDs = [i[0] for i in data]
    UNI_IDs = [UNI_IDs[i:i + n] for i in range(0, len(UNI_IDs), n)]
    urls = [
        f"https://www.uniprot.org/uniprot/?query=id%3A{'+OR+id%3A'.join(i)}&format=tab&columns=id,database(Pfam),organism-id&sort=score"
        for i in UNI_IDs]

    parsed_data = []
    for i in tqdm(urls):
        temp_h = request.urlopen(i)
        temp = [i0.decode("utf-8") for i0 in temp_h.readlines()[1:]]
        parsed_data.extend(temp)

    if filename:
        newfile = open(filename, "w")
    else:
        newfile = open("temp_files/uni_parsed_data.tsv", "w")

    for i in parsed_data:
        newfile.write(i)
    newfile.close()

    parsed_data = [i.rstrip().split("\t") for i in parsed_data]
    return parsed_data


def check_all_families(Uni_data):
    result = {"No_family": []}
    result_separated = {"No_family": []}
    for i0 in tqdm(range(len(Uni_data))):
        if len(Uni_data[i0][1]) == 0:
            result["No_family"].append(Uni_data[i0][0])
            result_separated["No_family"].append(Uni_data[i0][0])
        else:
            temp = Uni_data[i0][1].split(";")
            for i1 in temp:
                if i1 not in result_separated:
                    result_separated[i1] = [Uni_data[i0][0]]
                else:
                    result_separated[i1].append(Uni_data[i0][0])
            if "_".join(temp) in result:
                result["_".join(temp)].append(Uni_data[i0][0])
            else:
                result["_".join(temp)] = [Uni_data[i0][0]]

    families_list = list(result_separated.keys())
    families_list.sort(key=lambda i: len(result_separated[i]), reverse=True)
    for i in families_list:
        print(len(result_separated[i]), i, result_separated[i])

    families_list = list(result.keys())
    families_list.sort(key=lambda i: len(result[i]), reverse=True)
    for i in families_list:
        print(len(result[i]), i, result[i])


## UNUSED CODE - FASTER GENERAL FAMILY CHECK
    # if not parsed_tsv or not os.path.exists(parsed_tsv):
    #     Uni_data = parse_data_from_Uniprot(data, parsed_tsv)
    # elif os.path.exists(parsed_tsv):
    #     file_h = open(parsed_tsv)
    #     file = file_h.readlines()
    #     file_h.close()
    #     Uni_data = [i.rstrip().split("\t") for i in file]
    # check_all_families(Uni_data)


# -----------------------------------------------------------------------------------------------------
