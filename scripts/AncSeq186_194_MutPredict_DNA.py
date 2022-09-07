# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:51:08 2021

@author: Liu Yunsong

E-mail: liuyunsong@ioz.ac.cn 

Raising my dick and asking God who fucks!!!
"""


from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestClassifier
import seaborn as sns


Work_Dir = "E:/Project/Flu/Data/H3N2New/Uniq/FilterGene/HA1/"
Output_Dir = "E:/Project/Flu/Script/TestConclusion/"


AllKindsBases_toATCG = {"A":["A"], "T":["T"], "C":["C"], "G":["G"], 
                        "R":["A", "G"], "Y":["C", "T"], "M":["A", "C"], 
                        "K":["G", "T"], "S":["C", "G"], "W":["A", "T"], 
                        "H":["A", "C", "T"], "B":["C", "G", "T"], 
                        "V":["A", "C", "G"], "D":["A", "G", "T"], 
                        "N":["A", "C", "G", "T"], "?":["A", "C", "G", "T"]}


UniqID_Seq = {}
with open(Work_Dir+"HA1.fasta", "r+")as File:
    for per_line in File.readlines():
        if ">" in per_line:
            Name = per_line.replace("\n","").replace(">","")
            Seq = ""
        else:
            Seq += per_line.replace("\n","")
        UniqID_Seq[Name] = Seq 


Required_Codon = [Pos for Pos in range(1,330,1)]
for index, CodonPos in enumerate(Required_Codon):
    for Tip, Seq in UniqID_Seq.items():
        Codon = list(Seq[(CodonPos-1)*3:CodonPos*3])
        if len(AllKindsBases_toATCG[Codon[0]])>1:
            Codon[0] = np.random.choice(list(set(AllKindsBases_toATCG[Codon[0]]) - set("T")))
            UniqID_Seq[Tip] = Seq[:(CodonPos-1)*3] + "".join(Codon) + Seq[CodonPos*3:]
        if len(AllKindsBases_toATCG[Codon[1]])>1:
            Codon[1] = np.random.choice(list(set(AllKindsBases_toATCG[Codon[1]])))
            UniqID_Seq[Tip] = Seq[:(CodonPos-1)*3] + "".join(Codon) + Seq[CodonPos*3:]
        if len(AllKindsBases_toATCG[Codon[2]])>1:
            Codon[2] = np.random.choice(list(set(AllKindsBases_toATCG[Codon[2]])))
            UniqID_Seq[Tip] = Seq[:(CodonPos-1)*3] + "".join(Codon) + Seq[CodonPos*3:]


with open(Work_Dir+"HA1_PassageType_IDList.txt","r+")as File:
    PassageType_IDList = eval(File.readlines()[0])


treefile = Work_Dir + "HA1TrimTree.newick"
with open(treefile, "r+")as File:
    NewTree = File.readlines()[0]


EGG = PassageType_IDList["EGG"]


t = Tree(NewTree)


def get_terminal_descendants(anc_node):
    if len(anc_node.children) != 2:
        return [anc_node]
    else:
        return get_terminal_descendants(anc_node.children[0]) + \
            get_terminal_descendants(anc_node.children[1])


def check_whether_uniq_category(node, node_from, category):
    if node_from == "Left":
        tip_nodes = get_terminal_descendants(node.children[1])
    elif node_from == "Right":
        tip_nodes = get_terminal_descendants(node.children[0])
    for node in tip_nodes:
        if node.name not in category:
            return 1, tip_nodes
    else:
        return -1, tip_nodes


def get_uniq_category_top_anc(node, category):
    children = [node]
    other_children = []
    while True:
        anc = node.up
        if node == anc.children[0]:
            node_from = "Left"
        elif node == anc.children[1]:
            node_from = 'Right'
        check, other_children = check_whether_uniq_category(anc, node_from, category)
        if check == 1:
            return anc, children
        elif check == -1:
            children += other_children
        node = anc


EGG_Nodes = []
EGG_CommonAnc = []
Anc_Nodes = {}


for node in t:
    if node.name in EGG and node not in EGG_Nodes:
        CommonAnc, TipNodes = get_uniq_category_top_anc(node, EGG)
        EGG_Nodes += TipNodes
        EGG_CommonAnc.append(CommonAnc)
        Anc_Nodes[CommonAnc] = TipNodes


retree = NewTree.replace(";", "").replace("(",",").replace(")",",").split(",")
for i, j in enumerate(retree):
    if ":" in j:
        retree[i] = j.split(":")[0]
retree = ",".join(retree)
Anc_Pos = {}
for Anc in EGG_CommonAnc:
    subt = Anc.write(format=1).replace(";","").replace("(",",").replace(")",",").split(",")
    for i, j in enumerate(subt):
        if ":" in j:
            subt[i] = j.split(":")[0]
        else:
            subt[i] = ""
    subt = ",".join(subt)
    pos = retree[:retree.find(subt)+len(subt)].count(",")
    Anc_Pos[Anc] = pos


Anc_CodonPos_DNA = {Anc:{} for Anc, Pos in Anc_Pos.items()}
for i in range(1,330,1):
    with open("E:/Project/Flu/Output/MutationalMapping/HA1/1/"+str(i)+"/MutHistory.txt", "r+")as File:
        his = File.readlines()[0]
        his = his.replace("Codon Position:"+str(i)+"=","").replace(" ","").replace("(",",").replace(")",",").split(",")
    for Anc, Pos in Anc_Pos.items():
        Mut = his[Pos]
        if "-" in Mut:
            Mut = Mut.split("-")[-1]
        else:
            Mut = Mut
        Anc_CodonPos_DNA[Anc].update({i:Mut})


Seq186Yes = []
Seq194Yes = []
Seq186No = []
Seq194No = []
Seq186Un = []
Seq194Un = []
CodonSite_DNAList = {Site:[] for Site in range(1,330,1)}
for Anc, TipNodes in Anc_Nodes.items():
    AncDNA = list(Anc_CodonPos_DNA[Anc].values())
    for CodonSite in range(329):
        CodonSite_DNAList[CodonSite+1].append(AncDNA[CodonSite])
    _186 = []
    _194 = []
    _186Seq = []
    _194Seq = []
    for Node in TipNodes:
        Seq = UniqID_Seq[Node.name]
        ChildDNA = []
        for CodonSite in range(329):
            ChildDNA.append(Seq[CodonSite*3:CodonSite*3+3])
        if AncDNA[185] != ChildDNA[185]:
            _186.append("N")
            _186Seq.append(AncDNA)
        else:
            _186.append("Y")
            _186Seq.append(AncDNA)
        if AncDNA[193] != ChildDNA[193]:
            _194.append("N")
            _194Seq.append(AncDNA)
        else:
            _194.append("Y")
            _194Seq.append(AncDNA)
    if len(set(_186)) == 1:
        if "Y" in _186:
            Seq186Yes += _186Seq
        else:
            Seq186No += _186Seq
    else:
        Seq186Un += _186Seq
    if len(set(_194)) == 1:
        if "Y" in _194:
            Seq194Yes += _194Seq
        else:
            Seq194No += _194Seq
    else:
        Seq194Un += _194Seq


Pos186 = list(range(1,330,1))
Pos194 = list(range(1,330,1))


target186 = []
target194 = []


ValueMatrix186 = []
for Seq in Seq186Yes:
    temp = []
    for pos in Pos186:
        temp.append(Seq[pos-1])
    ValueMatrix186.append(temp)
    target186.append("Y")
for Seq in Seq186No:
    temp = []
    for pos in Pos186:
        temp.append(Seq[pos-1])
    ValueMatrix186.append(temp)
    target186.append("N")
for Seq in Seq186Un:
    temp = []
    for pos in Pos186:
        temp.append(Seq[pos-1])
    ValueMatrix186.append(temp)
    target186.append("Un")


ValueMatrix194 = []
for Seq in Seq194Yes:
    temp = []
    for pos in Pos194:
        temp.append(Seq[pos-1])
    ValueMatrix194.append(temp)
    target194.append("Y")
for Seq in Seq194No:
    temp = []
    for pos in Pos194:
        temp.append(Seq[pos-1])
    ValueMatrix194.append(temp)
    target194.append("N")
for Seq in Seq194Un:
    temp = []
    for pos in Pos194:
        temp.append(Seq[pos-1])
    ValueMatrix194.append(temp)
    target194.append("Un")


DF186 = pd.DataFrame(ValueMatrix186, columns=Pos186)
DF194 = pd.DataFrame(ValueMatrix194, columns=Pos194) 


enc186 = OneHotEncoder(categories="auto").fit(DF186)
result186 = enc186.transform(DF186)


enc194 = OneHotEncoder(categories="auto").fit(DF194)
result194 = enc194.transform(DF194)


Times = 10


###########
## 决策树 ##
###########


tree186 = []
tree_verify186 = []
Pos_Importance186Tree = {Pos:[] for Pos in list(range(1,330,1))}
for i in range(Times):
    print(i)
    feature = Pos186
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(result186, 
                                                    target186, 
                                                    test_size=0.3)
    clf1 = tree.DecisionTreeClassifier()
    clf1 = clf1.fit(Xtrain, Ytrain)
    tree186.append(clf1.score(Xtest, Ytest))
    tree_verify186.append(cross_val_score(clf1,result186,target186,cv=100).mean())
    for Pos, Importance in zip(list(range(1,330,1)), clf1.feature_importances_):
        if Importance > 0:
            Pos_Importance186Tree[Pos].append(Importance)


tree194 = []
tree_verify194 = []
Pos_Importance194Tree = {Pos:[] for Pos in list(range(1,330,1))}
for j in range(Times):
    print(j)
    feature = Pos194
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(result194, 
                                                    target194, 
                                                    test_size=0.3)
    clf2 = tree.DecisionTreeClassifier()
    clf2 = clf2.fit(Xtrain, Ytrain)
    tree194.append(clf2.score(Xtest, Ytest))
    tree_verify194.append(cross_val_score(clf2,result194,target194,cv=100).mean())
    Pos_Importance194 = {}
    for Pos, Importance in zip(list(range(1,330,1)), clf2.feature_importances_):
        if Importance > 0:
            Pos_Importance194Tree[Pos].append(Importance)


#############
## 随机森林 ##
#############


forest186 = []
forest_verify186 = []
Pos_Importance186Forest = {Pos:[] for Pos in list(range(1,330,1))}
for i in range(Times):
    print(i)
    feature = Pos186
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(result186, 
                                                    target186, 
                                                    test_size=0.3)
    clf1 = RandomForestClassifier()
    clf1 = clf1.fit(Xtrain, Ytrain)
    forest186.append(clf1.score(Xtest, Ytest))
    forest_verify186.append(cross_val_score(clf1,result186,target186,cv=100).mean())
    for Pos, Importance in zip(list(range(1,330,1)), clf1.feature_importances_):
        if Importance > 0:
            Pos_Importance186Forest[Pos].append(Importance)


forest194 = []
forest_verify194 = []
Pos_Importance194Forest = {Pos:[] for Pos in list(range(1,330,1))}
for j in range(Times):
    print(j)
    feature = Pos194
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(result194, 
                                                    target194, 
                                                    test_size=0.3)
    clf2 = RandomForestClassifier()
    clf2 = clf2.fit(Xtrain, Ytrain)
    forest194.append(clf2.score(Xtest, Ytest))
    forest_verify194.append(cross_val_score(clf2,result194,target194,cv=100).mean())
    for Pos, Importance in zip(list(range(1,330,1)), clf2.feature_importances_):
        if Importance > 0:
            Pos_Importance194Forest[Pos].append(Importance)


Classifier_DF = pd.DataFrame({"Accuracy":tree186+tree_verify186+tree194+tree_verify194+
                                          forest186+forest_verify186+forest194+forest_verify194,
                              "Methods":["DecisionTreeClassifier\n186"]*Times+
                                        ["DecisionTreeClassifier\ncross val score\n186"]*Times+
                                        ["DecisionTreeClassifier\n194"]*Times+
                                        ["DecisionTreeClassifier\ncross val score\n194"]*Times+
                                        ["RandomForestClassifier\n186"]*Times+
                                        ["RandomForestClassifier\ncross val score\n186"]*Times+
                                        ["RandomForestClassifier\n194"]*Times+
                                        ["RandomForestClassifier\ncross val score\n194"]*Times
                                        })


fig, axs = plt.subplots(1, 1, 
                    sharex=False, sharey=False, 
                    figsize=(25, 15), 
                    dpi=200)


plt.tick_params(labelsize=20)
sns.set_theme(style="ticks", palette="pastel")
sns.boxenplot(x="Accuracy", 
              y="Methods", 
              data=Classifier_DF)
# sns.despine(offset=10, trim=True)
plt.xlabel("")
plt.ylabel("")
plt.savefig("E:/Project/Flu/Script/TestConclusion/Classifier.svg", dpi=500)


























# Times = 200


# tree186 = []
# forest186 = []
# tree_verify186 = []
# forest_verify186 = []
# Pos_Importance186Tree = {Pos:[] for Pos in list(range(1,330,1))}
# Pos_Importance186Forest = {Pos:[] for Pos in list(range(1,330,1))}
# for i in range(Times):
#     print(i, 186)
#     feature = Pos186
#     Xtrain, Xtest, Ytrain, Ytest = train_test_split(result186, 
#                                                     target186, 
#                                                     test_size=0.3)
#     clf1 = tree.DecisionTreeClassifier()
#     clf1 = clf1.fit(Xtrain, Ytrain)
#     tree186.append(clf1.score(Xtest, Ytest))
#     tree_verify186.append(cross_val_score(clf1,result186,target186,cv=30).mean())
#     for Pos, Importance in zip(list(range(1,330,1)), clf1.feature_importances_):
#         Pos_Importance186Tree[Pos].append(Importance)

#     clf2 = RandomForestClassifier()
#     clf2 = clf2.fit(Xtrain, Ytrain)
#     forest186.append(clf2.score(Xtest, Ytest))
#     forest_verify186.append(cross_val_score(clf2,result186,target186,cv=30).mean())
#     for Pos, Importance in zip(list(range(1,330,1)), clf2.feature_importances_):
#         Pos_Importance186Forest[Pos].append(Importance)


# tree194 = []
# tree_verify194 = []
# Pos_Importance194Tree = {Pos:[] for Pos in list(range(1,330,1))}
# forest194 = []
# forest_verify194 = []
# Pos_Importance194Forest = {Pos:[] for Pos in list(range(1,330,1))}
# for j in range(Times):
#     print(j, 194)
#     feature = Pos194
#     Xtrain, Xtest, Ytrain, Ytest = train_test_split(result194, 
#                                                     target194, 
#                                                     test_size=0.3)
#     clf3 = tree.DecisionTreeClassifier()
#     clf3 = clf3.fit(Xtrain, Ytrain)
#     tree194.append(clf3.score(Xtest, Ytest))
#     tree_verify194.append(cross_val_score(clf3,result194,target194,cv=30).mean())
#     Pos_Importance194 = {}
#     for Pos, Importance in zip(list(range(1,330,1)), clf3.feature_importances_):
#         Pos_Importance194Tree[Pos].append(Importance)

#     clf4 = RandomForestClassifier()
#     clf4 = clf4.fit(Xtrain, Ytrain)
#     forest194.append(clf4.score(Xtest, Ytest))
#     forest_verify194.append(cross_val_score(clf4,result194,target194,cv=30).mean())
#     for Pos, Importance in zip(list(range(1,330,1)), clf4.feature_importances_):
#         Pos_Importance194Forest[Pos].append(Importance)


# Classifier_DF = pd.DataFrame({"Accuracy":tree186+\
#                                           tree_verify186+
#                                          tree194+\
#                                           tree_verify194+
#                                          forest186+\
#                                           forest_verify186+
#                                          forest194+
#                                           forest_verify194,
#                               "Methods":["DecisionTreeClassifier\n186"]*Times+\
#                                         ["DecisionTreeClassifier\ncross val score\n186"]*Times+
#                                         ["DecisionTreeClassifier\n194"]*Times+\
#                                         ["DecisionTreeClassifier\ncross val score\n194"]*Times+
#                                         ["RandomForestClassifier\n186"]*Times+\
#                                         ["RandomForestClassifier\ncross val score\n186"]*Times+
#                                         ["RandomForestClassifier\n194"]*Times+
#                                         ["RandomForestClassifier\ncross val score\n194"]*Times
#                                         })
# Classifier_DF = pd.DataFrame({"Accuracy":tree186+
#                                          tree194+
#                                          forest186+
#                                          forest194,
#                               "Methods":["DecisionTreeClassifier\n186"]*Times+
#                                         ["DecisionTreeClassifier\n194"]*Times+
#                                         ["RandomForestClassifier\n186"]*Times+
#                                         ["RandomForestClassifier\n194"]*Times
#                                         })


# fig, axs = plt.subplots(1, 1, 
#                     sharex=False, sharey=False, 
#                     figsize=(25, 15), 
#                     dpi=200)


# plt.tick_params(labelsize=20)
# sns.set_theme(style="ticks", palette="pastel")
# sns.boxenplot(x="Accuracy", 
#               y="Methods", 
#               data=Classifier_DF)
# plt.xlabel("")
# plt.ylabel("")
# plt.savefig("E:/Project/Flu/Script/TestConclusion/Classifier.svg", dpi=500)