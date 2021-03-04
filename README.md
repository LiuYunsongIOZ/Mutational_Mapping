# Mutational_Mapping(Molecular Biology)
Academic Project (Adaptive Evolution)

This is an academic project started by Dr. Liu Yunsong from Institute of Zoology, Chinese Academy of Sciencies. Our main goal is to help researchers conduct research on adaptive evolution. We hope this tool will be helpful to you. This tool called "Mutaional Mapping" will mainly use Python to finish. Please be free to contact me (liuyunsong@ioz.ac.cn), if you some advices for this tool or some problems when you using it.  


## DEPENDENCES
Mutational Mapping requires the following dependencies:
- Python (>= 3.7)
- numpy(>=1.19.2)
- scipy(>=1.5.2)
- pprint
- pandas
- math
- random
- datetime
- sys

## WARNING
When this program check your input files, it merely can check whether the specified path exist these files and the files' format, it can not do more accurate work. So you need to check your files by yourself before you using this tool.
1. Firstly, you need to make sure that the Newick format phylogenetic tree file you provide has and only contains one tree (https://en.wikipedia.org/wiki/Newick_format), and this phylogenetic tree must be a perfect binary tree.
For example : 

              "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"    ×
              
              "((A:0.1,B:0.2):0.0,(C:0.3,D:0.4):0.5);"    √
2. The sequences file should be saved in fasta format, ">id\nATCG". The id represents the sequence identifier which you can change according your requirements and ATCG represents the sequence.

3. All input files that need to be analyzed, please strictly follow the sample file format, can not randomly add any content, please refer to the files in the "/Mutational_Mapping/examples/" folder for details.

## USER GUIDE
### TESTING

### PARAMETERS

| parameter | /examples/file | introduction |
| :------------- | :------------- | :------------- |
| --tree | treefile.newick | a perfect binary newick format phylogenetic tree |
| --category | category_IDs.txt | sequences' category |
| --fasta | data.fasta | fasta format data |
| --codon | Required_Codon.txt | |
| --Qrate | GTR_Q_rate.txt | The rate parameter of Q Matrix of GTR model  ( A-C: d, A-G: f, A-T: b, C-G: e, C-T: a, G-T: c) |
| --Gamma | Gamma_Rate.txt | Gamma rate |
| --codon | SeqTimeInfor.txt| |


## Distribution Structure
