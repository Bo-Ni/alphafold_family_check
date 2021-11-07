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

    final_result = []
    for i0 in range(len(data)):
        domain_result = check_knot_families(data[i0][0], [data[i0][1], data[i0][2]])
        final_result.append([data[i0], domain_result])

    return final_result


if __name__ == "__main__":
    AlphaFoldKnotDomainCheck([["Q9FLD5", 1, 100], ["Q6PL18", 1, 100], ["Q5T9A4", 1, 100]])
