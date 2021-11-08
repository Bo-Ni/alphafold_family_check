from urllib import request
import asyncio

async def AlphaFoldStructureCheck(input_data):
    """
    :param input_data: List of Uniprot IDs ["A0A081RPU7", "A0A087S697",	"A0A087S2L9", ...]
    :return:
    """

    def CheckPDBExists(input_data):
        """
        Checks for PDB IDs of given protein in Uniprot.
        WARNING - Script assumes batches are small enough not to reach the url length limit.
        Anything below 150 Uni IDs should be fine.

        :param input_data: List of Uniprot IDs ["A0A081RPU7", "A0A087S697",	"A0A087S2L9", ...]
        :return: no_pdb: List of Uniprot IDs without structure ['Q9FLD5', 'Q5T9A4', 'Q925I1', ...]
        :return: pdb_exists: List of Uniprot IDs and all structures [['Q12019', ['5FL8', '5JCS']], ['P32795', ['2MV3']], ...]
        """

        insert_text = "%20OR%20id:".join(input_data)
        url = f"https://www.uniprot.org/uniprot/?query=id:{insert_text}&format=tab&columns=id,database(PDB)&sort=score"
        pdb_exists = request.urlopen(url)
        pdb_exists = pdb_exists.readlines()
        pdb_exists = [i.decode("utf-8").rstrip().split("\t") for i in pdb_exists[1:]]
        no_pdb = [i[0] for i in pdb_exists if len(i) == 1]
        pdb_exists = [i for i in pdb_exists if len(i) == 2]
        pdb_exists = [[i[0], "prot", i[1].rstrip(";").split(";")] for i in pdb_exists]

        return no_pdb, pdb_exists

    def CheckHomologExists(input_data):
        insert_text = "%20OR%20uniprot:".join(input_data)
        url = f"https://www.uniprot.org/uniref/?query=uniprot:{insert_text}&format=tab&limit=10&columns=id,count,members&sort=score"
        clusters = request.urlopen(url)
        clusters = clusters.readlines()
        clusters = [i.decode("utf-8").rstrip().split("\t") for i in clusters[1:]]
        result = []

        for UniID in input_data:
            temp_clusters = [i for i in clusters if UniID in i[0]]
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
                all_PDBS_exist = []
                all_no_PDBS = []
                for i1 in batches:
                    temp_no_pdb, temp_exists = CheckPDBExists(i1)
                    all_no_PDBS.extend(temp_no_pdb)
                    all_PDBS_exist.extend(temp_exists)
                result.append([UniID, identity, cluster_size, len(all_PDBS_exist)])
            else:
                result.append([UniID, 'no_data', 0, 0])

        return result

    no_pdb, pdb_exists = CheckPDBExists(input_data)
    homolog_result = CheckHomologExists(no_pdb)

    result = pdb_exists
    result.extend(homolog_result)
    return result


if __name__ == "__main__":
    AlphaFoldStructureCheck(["Q9FLD5", "Q6PL18", "Q5T9A4", "Q925I1", "P32795", "Q12019"])

