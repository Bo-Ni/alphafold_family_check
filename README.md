# Alpha Fold Domain Check

## Current workflow
![alt text](https://github.com/exsto1/alpha_fold_wrapper/blob/7b6eaa8bfa0a14e85901674893ab3cac8830503f/info/alpha_fold.png)

## Usage
```
from alpha_fold_main import AlphaFoldWorkflow
AlphaFoldWorkflow(INPUT_DATA, savefile="PATH/TO/SAVEFILE")
```
INPUT_DATA in this contex is a list with Uniprot IDs, start and end indexes (for example stating range of the knot).
See an example below:

```[["Q9FLD5", 1, 100], ["Q6PL18", 1, 100], ["Q5T9A4", 1, 100]]```

"savefile" flag is especially useful for big datasets, enabling to resume search after unexpected errors, saving data in batches of 10 IDs.
If none provided script will use default value of "savefile.txt" - link this file in the flag if you want to use it as a checkpoint after a crash.

## Savefile format and it's unpacked version
Sample line from the sivefile looks like this:
```P22748_48_285;50hom;179;PF00194_92;prot;10;1ZNC,3F7B,3F7U,3FW3,5IPZ,5JN8,5JN9,5JNA,5JNC,5KU6```
It contains following data separated with ";" character:
- Uniprot ID, start range, end range (as supplied into the workflow) connected with "_" character
- Data set for domain info, one of the following:
    - prot - info from the input protein
    - 50hom - no info from input protein, using UniRef 50%
    - 90hom - no info from input protein and no UniRef 50% dataset, using UniRef 90%
    - nodata - no info from input protein and no UniRef datasets
- Number of homologs found. Artificial 0 for input protein data.
- Depending if dataset comes from protein itself or it's homologs:
    - prot: list of domain info, separated by "," character. Each section contains Pfam ID, start and end indexes connected by "_" character.
    - homologs: list of domain info, separated by "," character. Each section contains Pfam ID and number of homologs it was found for, connected by "_" character.
- Data set for structure info, same as in section 2.
- Number of structures found.
- List of found structure - for protein or collective list for all homologs.

So as listed above, savefile file contains 7 main sections in total. 
So does the parsed list from the alpha_fold_analyze_results.DecodeSavefile function.
All the separators were split, but the structure stays the same.
Example list can be seen below:

``` [['P07451', '25', '259'], '50hom', '77', [['PF00194', '44']], 'prot', '5', ['1Z93', '1Z97', '2HFW', '3UYN', '3UYQ']] ```
