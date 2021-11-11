"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    return open(filename,"r").read().replace("\n","")


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    rna=[]
    dna=dna.replace("T","U")
    for i in range(startIndex,len(dna),3):
        rna.append(dna[i:i+3])
        if dna[i:i+3]  in ["UAA","UAG","UGA"]:
            return rna         
    return rna


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    codonDict={}
    acids_codons=json.load(open(filename,"r"))
    for acid,codon in acids_codons.items():
        for i in codon:
            codonDict[i.replace("T","U")]=acid          
    return codonDict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein=["Start"]
    for i in range(1,len(codons)):
        protein.append(codonD[codons[i]])
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna=readFile(dnaFilename)
    proteins=[]
    count,i=0,0
    codonDict=makeCodonDictionary(codonFilename)
    while i<len(dna):
        if dna[i:i+3]=="ATG":
            rna=dnaToRna(dna,i)
            proteins.append(generateProtein(rna,codonDict))
            i=i+(3*len(rna))
        else:
            count=count+1
            i=i+1
    print("\nThere are",len(dna), "bases in the dna")
    print("There are",len(proteins),"synthesized proteins\nand",count,"unused bases")
    return proteins

def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    commonPro=[]
    for i in proteinList1:
        if i in proteinList2 and i not in commonPro:
            commonPro.append(i) 
    return commonPro


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):   
    return [aminoacid for protein in proteinList for aminoacid in protein ]

'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoDic={}
    for i in aaList:
        if i not in aminoDic:
            aminoDic[i]=0
        aminoDic[i]+=1
    return aminoDic


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    totalAmino1,totalAmino2=combineProteins(proteinList1),combineProteins(proteinList2)
    aminoAcidDic1,aminoAcidDic2=aminoAcidDictionary(totalAmino1),aminoAcidDictionary(totalAmino2)
    aminoAcidDiff=[]
    lst=list(aminoAcidDic1.keys())+list(aminoAcidDic2.keys())
    uniqueAmino = list(dict.fromkeys(lst))
    for amino in uniqueAmino:
            if amino not in ["Start","Stop"]:
                if amino in aminoAcidDic1:
                    proteinList1Freq=aminoAcidDic1[amino]/len(totalAmino1)
                else:
                    proteinList1Freq=0
                if amino in aminoAcidDic2:
                    proteinList2Freq=aminoAcidDic2[amino]/len(totalAmino2)
                else:
                    proteinList2Freq=0
                dummylist=[amino,proteinList1Freq,proteinList2Freq]
                if abs(proteinList1Freq-proteinList2Freq)>cutoff and dummylist not in aminoAcidDiff:
                    aminoAcidDiff.append(dummylist)
    return aminoAcidDiff


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:\n")
    commonProteins=[]
    for protein in commonalities:
        if len(protein)>2:
            commonProteins.append(protein[1:len(protein)-1])
    commonProteins.sort(key= lambda proteins:proteins[0])
    for i in commonProteins:
        if len(i)==1:
            print(i[0],"\n")
        else:
            stringText=""
            for amino in i:
                if amino!=i[len(i)-1]:
                    stringText+=amino+"-"
                else:
                    stringText+=amino
            print(stringText,"\n")        
    print("The following amino acids occurred at very different rates in the two DNA sequences:\n")
    for aminoDiff in differences:
        rates=aminoDiff[0]+" : "+str(round(aminoDiff[1]*100,2)) \
        +"% in Seq1, "+str(round(aminoDiff[2]*100,2))+"% in Seq2\n"
        print(rates)
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities =  commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    plist1=[amino for protein in proteinList1 for amino in protein]
    [plist1.append(amino) for protein in proteinList2 for amino in protein]
    plist1=list(dict.fromkeys(plist1))
    return sorted(plist1)

'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    freq=[]
    L=combineProteins(proteinList)
    D=aminoAcidDictionary(L)
    for amino in labels:
        freq.append(D[amino]/len(L)) if amino in D else freq.append(0)
    return freq

'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testGenerateProtein()
    test.testSetupChartData()
    # test.testMakeAminoAcidLabels()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
