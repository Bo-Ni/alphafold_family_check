from tqdm import tqdm
from urllib import request
import xmltodict
import os


def alpha_fold_knot_family_check(data, parsed_tsv=None, savefile=None):
    """
    :param data: list of [Uniprot IDs, start_index, end_index]
    :param parsed_tsv: path to parsed Uniprot table containing: UniID, PfamFamilyID, OrganismID. If file doesn't exist, it will be created.
    :param savefile: path to temporary calculations saved in batches, made to recover after crash. If file doesn't exist, it will be created.
    """

    def parse_data_from_Uniprot(data, filename):
        n=50
        UNI_IDs = [i[0] for i in data]
        UNI_IDs = [UNI_IDs[i:i+n] for i in range(0, len(UNI_IDs), n)]
        urls = [f"https://www.uniprot.org/uniprot/?query=id%3A{'+OR+id%3A'.join(i)}&format=tab&columns=id,database(Pfam),organism-id&sort=score" for i in UNI_IDs]

        parsed_data = []
        for i in tqdm(urls):
            temp_h = request.urlopen(i)
            temp = [i0.decode("utf-8") for i0 in temp_h.readlines()[1:]]
            parsed_data.extend(temp)

        if filename:
            newfile = open(filename, "w")
        else:
            newfile = open("uni_parsed_data.tsv", "w")

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

    def check_knot_families(uni_id, dom_range):
        def single_check(single, dom_range):
            start = int(single["location"]["@start"])
            end = int(single["location"]["@end"])
            pfam = single["@accession"]
            if start <= dom_range[0] and end >= dom_range[1] or dom_range[1] >= start >= dom_range[0]:
                return [pfam, start, end]

        temp = request.urlopen(f"http://pfam.xfam.org/protein/{uni_id}?output=xml")
        data = xmltodict.parse(temp)
        result = []
        try:
            if type(data["pfam"]["entry"]["matches"]["match"]) == type([]):
                for i in data["pfam"]["entry"]["matches"]["match"]:
                    single_result = single_check(i, dom_range)
                    if single_result:
                        result.append(single_result)
            else:
                single_check(data["pfam"]["entry"]["matches"]["match"], dom_range)
        except:
            pass

        return result

    def continue_search(data, savefile_name):
        temp_h = open(savefile_name)
        temp = temp_h.readlines()
        temp_h.close()

        newdata = []
        for i in data:
            t_i = "_".join([str(i0) for i0 in i])
            state = False
            for i1 in temp:
                if t_i in i1:
                    state = True
            if not state:
                newdata.append(i)

        return newdata

    ## UNUSED CODE - FASTER GENERAL FAMILY CHECK
    # if not parsed_tsv or not os.path.exists(parsed_tsv):
    #     Uni_data = parse_data_from_Uniprot(data, parsed_tsv)
    # elif os.path.exists(parsed_tsv):
    #     file_h = open(parsed_tsv)
    #     file = file_h.readlines()
    #     file_h.close()
    #     Uni_data = [i.rstrip().split("\t") for i in file]
    # check_all_families(Uni_data)

    if not savefile or not os.path.exists(savefile):
        savefile_name = "savefile.txt"
    elif os.path.exists(savefile):
        savefile_name = savefile
        data = continue_search(data, savefile)

    final_result = []
    batch_size = 10
    data = [data[i:i + batch_size] for i in range(0, len(data), batch_size)]

    for i0 in tqdm(range(len(data))):
        temp_batch_string = []
        for i in range(len(data[i0])):
            domain_result = check_knot_families(data[i0][i][0], [data[i0][i][1], data[i0][i][2]])
            temp_string = ";".join(["_".join([str(i0) for i0 in data[i0][i]]),
                                    ",".join(["_".join([str(i0) for i0 in i]) for i in domain_result])]) + "\n"
            temp_batch_string.append(temp_string)
            final_result.append([data[i0][i], domain_result])
        temp_file = open(savefile_name, "a")
        temp_file.write("".join(temp_batch_string))
        temp_file.close()


    return final_result
