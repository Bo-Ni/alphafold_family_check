from tqdm import tqdm
from urllib import request
import xmltodict
import os
import multiprocessing
import asyncio


async def AlphaFoldKnotDomainCheck(data):
    """
    :param data: list of [Uniprot IDs, start_index, end_index]
    """

    def CheckPfamFamiliesExist(uni_id, dom_range):
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

    def CheckUniProtFamiliesExist(input_data):
        insert_text = "%20OR%20id:".join(input_data)
        url = f"https://www.uniprot.org/uniprot/?query=id:{insert_text}&format=tab&columns=id,database(Pfam)&sort=score"
        family_exists = request.urlopen(url)
        family_exists = family_exists.readlines()
        family_exists = [i.decode("utf-8").rstrip().split("\t") for i in family_exists[1:]]
        no_family = [i[0] for i in family_exists if len(i) == 1]
        family_exists = [i for i in family_exists if len(i) == 2]
        family_exists = [[i[0], "prot", i[1].rstrip(";").split(";")] for i in family_exists]

        return no_family, family_exists

    def CheckHomologExists(input_data):
        insert_text = "%20OR%20uniprot:".join([i[0] for i in input_data])

        url = f"https://www.uniprot.org/uniref/?query=uniprot:{insert_text}&format=tab&limit=10&columns=id,count,members&sort=score"
        clusters = request.urlopen(url)
        clusters = clusters.readlines()
        clusters = [i.decode("utf-8").rstrip().split("\t") for i in clusters[1:]]
        result = []

        for UniID in input_data:
            temp_clusters = [i for i in clusters if UniID[0] in i[0]]
            chosen_cluster = None
            identity = None
            for i in temp_clusters:
                if "UniRef50" in i[0]:
                    chosen_cluster = i
                    identity = '50_hom'
                    break
                elif "UniRef90" in i[0]:
                    identity = '90_hom'
                    chosen_cluster = i

            if chosen_cluster:
                batch_size = 100
                cluster_size = int(chosen_cluster[1])
                chosen_cluster = chosen_cluster[2].split("; ")
                batches = [chosen_cluster[i:i + batch_size] for i in range(0, len(chosen_cluster), batch_size)]
                all_families_exist = []
                all_no_families = []
                for i1 in batches:
                    temp_no_family, temp_exists = CheckUniProtFamiliesExist(i1)
                    all_no_families.extend(temp_no_family)
                    all_families_exist.extend(temp_exists)
                all_families_exist = [i[-1] for i in all_families_exist]
                all_families_dict = {}
                for i1 in all_families_exist:
                    for i2 in i1:
                        if i2 in all_families_dict:
                            all_families_dict[i2] += 1
                        else:
                            all_families_dict[i2] = 1
                all_families_dict = [[i, all_families_dict[i]] for i in all_families_dict]
                UniID.extend([identity, cluster_size])
                result.append([UniID, all_families_dict])
            else:
                UniID.extend(['no_data', 0])
                result.append([UniID, []])

        return result

    final_result = []
    no_family = []
    for i0 in range(len(data)):
        domain_result = CheckPfamFamiliesExist(data[i0][0], [data[i0][1], data[i0][2]])
        if domain_result:
            data[i0].extend(["prot", 0])
            final_result.append([data[i0], domain_result])
        else:
            no_family.append(data[i0])

    homolog_check = CheckHomologExists(no_family)
    final_result.extend(homolog_check)
    return final_result




if __name__ == "__main__":
    AlphaFoldKnotDomainCheck([["Q9FLD5", 1, 100], ["Q6PL18", 1, 100], ["Q5T9A4", 1, 100]])
