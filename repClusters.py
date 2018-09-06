#!/usr/bin/env python

# Author: Jesper Svedberg (jesper.svedberg@ebc.uu.se)
# Version: 0.1


#import sys
import csv
from operator import itemgetter
import argparse
from Bio import SeqIO
from Bio import SeqUtils
import numpy
import scipy.stats as st

SIMPLE_TYPES = "IIIISS"
BTAB_TYPES = "SSISSSIIIIFFIIISISIIIS"
RMOUT_TYPES = "IFFFSIISSSSSSSIS"
MAKER_TYPES = "SSSIISSSS"
BISMARK_CX_TYPES="SISFISS" #Float should be Int. Changed for division
BISMARK_COV_TYPES="SIIFFI"
BED_GENOME_COVERAGE_TYPES="SII"
BED_CONCAT_4_TYPES="SIISIISIISII"
BED_CONCAT_7_TYPES="SIISIIIF" #Actually 8
BED_VALUE_TABLE_TYPES="I"*50
COVERAGE="SF"
BISMARK_SHORT="SIFI"

MET_COL=2
NON_COL=3
COV_CUTOFF = 5

REP_CUTOFF = 0.75

# Sorts a table in list or tuple form after several columns in reversed order
def sortTable(table, columns):
    for col in columns:
        table = sorted(table, key=itemgetter(col))
    return table

# Returns column from matrix
def getColumn(matrix,i):
    f = itemgetter(i)
    return map(f,matrix)


# Converts entries in a table from strings to whatever format that is specified in the type list,
# for instance integers or floats. The type list is now formatted for the MUMmer show-coords
# btab format.
def tableTypeConvert(table, typeList):
    outTable = []
    for line in table:
        newLine = []
        i = 0
        for entry in line:
            if typeList[i] == "S":
                newLine.append(entry)
            elif typeList[i] == "I":
                newLine.append(int(entry))
            elif typeList[i] == "F":
                newLine.append(float(entry))
            i += 1
        
        outTable.append(newLine)
        
    return outTable

def rowTypeConvert(line, typeList):
    newLine = []
    i = 0
    for entry in line:
        if typeList[i] == "S":
            newLine.append(entry)
        elif typeList[i] == "I":
            newLine.append(int(entry))
        elif typeList[i] == "F":
            newLine.append(float(entry))
        i += 1
    return newLine


# Turns a list into a tab delimited string
def tabulateList(inList):
    out = ""
    for i in inList:
        out += str(i) + "\t"
    return out[:-1]

def rmParser(table, lengths):
    outDict = {}
    for line in table:
        if line[4] in outDict:
            outDict[line[4]] = modifyContig(outDict[line[4]], line[5], line[6])
        else:
            outDict[line[4]] = modifyContig(createContig(lengths[line[4]]), line[5], line[6])
    return outDict
        
def createContig(length):
    return [0]*length

def modifyContig(contig, start, stop):
    for i in range(start-1,stop):
        contig[i] = 1
    return contig


def chipParser(table, lengths):
    outDict = {}
    for line in table:
        if line[0] in outDict:
            outDict[line[0]].append(line[1])
        else:
            print line[0]
            outDict[line[0]] = [line[1]]
    return outDict

def chipParserRow(dct, lengths, line):
    if line[0] in dct:
        dct[line[0]].append(line[1])
    else:
        print line[0]
        dct[line[0]] = [line[1]]
    return dct

def metParserRow(dct, lengths, line):
    if line[0] in dct:
        dct[line[0]] = calcMethylation(dct[line[0]], line[1], line[MET_COL], line[NON_COL])
    else:
        print line[0]
        dct[line[0]] = calcMethylation(createContig(lengths[line[0]]), line[1], line[MET_COL], line[NON_COL])
    return dct

def metParser(table, lengths):
    outDict = {}
    for line in table:
        if line[0] in outDict:
            outDict[line[0]] = calcMethylation(outDict[line[0]], line[1], line[MET_COL], line[NON_COL])
        else:
            print line[0]
            outDict[line[0]] = calcMethylation(createContig(lengths[line[0]]), line[1], line[MET_COL], line[NON_COL])
    return outDict

def calcMethylation(contig, pos, met, nonmet):
    if met+nonmet <= COV_CUTOFF:
        div = 0
    elif nonmet == 0:
        if met == 0:
            div = 0
        else:
            div = 1
    else:
        div = met/(met+nonmet)
    contig[pos-1] = div
    return contig

def calcMean(seq):
    if len(seq) == 0:
        return 0
    return sum(seq)/float(len(seq))

def gcCalc(seq):
    outSum = 0
    for i in "GCgcSs":
        outSum += str(seq).count(i)
    return outSum/float(len(seq))



parser = argparse.ArgumentParser(description='Script that calculate things for repeat clusters.')
parser.add_argument('filename', help='RepeatMasker file name')
parser.add_argument('-s','--size',help='Chromosome size. Default: Last line in file.', type=int)
parser.add_argument('-w','--window',help='Sliding window size. Default: 2000', type=int, default=2000)
parser.add_argument('-p','--step',help='Sliding window step length. Default: 2000', type=int, default=2000)
parser.add_argument('-g','--genome',help='Genome fasta file.')
parser.add_argument('-m','--methylation',help='Bismark CX file')
parser.add_argument('-c','--coverage',help='Coverage file. Modified bedtools coverage format.')
parser.add_argument('-o','--contig',help='Specify one contig to limit the analysis to.')
parser.add_argument('-b','--start',help='Start coordinate.')
parser.add_argument('-e','--end',help='Stop coordinate.')
parser.add_argument('-t','--outfile',help='Output file name. Default: [RepeatMasker_file].repClusters.csv')
args = parser.parse_args()


window = args.window
step = args.step

if args.contig:
    onlyContig = args.contig
    print "Only analyzing " + onlyContig
else:
    onlyContig = ""
    
if args.start:
    startcoord = int(args.start)
else:
    startcoord = 0
    
if args.end:
    stopcoord = int(args.end)
else:
    stopcoord = 0

if args.outfile:
    outfile = args.outfile
else:
    outfile = args.filename + ".repClusters.csv"
    

'''
if args.format.upper() == 'CX':
    tabFormat = BISMARK_CX_TYPES
    MET_COL=3
    NON_COL=4
elif args.format.upper() == 'COV':
    tabFormat = BISMARK_COV_TYPES
    MET_COL=4
    NON_COL=5
elif args.format.upper() == 'BED':
    tabFormat = BED_GENOME_COVERAGE_TYPES
    COV_COL=2
    COORD_COL=1
elif args.format.upper() == 'BC4':
    tabFormat = BED_CONCAT_4_TYPES
    COORD_COL=1
    COL1=2
    COL2=5
    COL3=8
    COL4=11
elif args.format.upper() == 'BC7':
    tabFormat = BED_CONCAT_7_TYPES
    COORD_COL=1
    COL1=2
    COL2=5
    SUB=6
    DIV=7
elif args.format.upper() == 'BVT':
    tabFormat = BED_VALUE_TABLE_TYPES
    COORD_COL=1
    COL1=2
    COL2=5
    SUB=6
    DIV=7
else:
    sys.exit("Wrong format specified.")
'''


##############
print "Reading genome..."
contigDict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

lengthDict = {}
for key, val in contigDict.iteritems():
    lengthDict[key] = len(val)

print lengthDict

##############
print "Reading repeats..."
table = [line.rstrip("\n").split() for line in open(args.filename)]

if table[0][0] == "SW":
    del table[0:3]

repDict = rmParser(tableTypeConvert(table, RMOUT_TYPES), lengthDict)

##############
print "Reading chipseq data..."
chipDict = {}

#chipDict = chipParser(tableTypeConvert([line.rstrip("\n").split("\t") for line in open(args.coverage)], COVERAGE), lengthDict)

for line in open(args.coverage):
    if onlyContig:
        if not line.startswith(onlyContig):
            continue
    
    ll = rowTypeConvert(line.rstrip("\n").split("\t"), COVERAGE)

    
    chipDict = chipParserRow(chipDict, lengthDict, ll)

#print lengthDict

##############
print "Reading methylation data..."
#metDict = metParser(tableTypeConvert([line for line in csv.reader(open(args.methylation),delimiter='\t')], BISMARK_SHORT), lengthDict)
metDict = {}
for line in open(args.methylation):
    if onlyContig:
        if not line.startswith(onlyContig):
            continue
        
    ll = rowTypeConvert(line.rstrip("\n").split("\t"), BISMARK_SHORT)
    
    
    metDict = metParserRow(metDict, lengthDict, ll)


##############

wholeGenome = []
repGenome = []
nonrepGenome = []

for key, contig in repDict.iteritems():
    if onlyContig:
        if key != onlyContig:
            continue
    
    print "Contig: " + key
    
    if stopcoord:
        lgLen = stopcoord
    else:
        lgLen = lengthDict[key]
    
    #lgLen = 500000
    
    winCoords = range(startcoord, lgLen-window, step)
    #windows = [sum(contig[i:i+window])/window for i in range(0, lgLen-window, step)]
    
    for i in winCoords:
        if contig[i:i+window]:
            rep = numpy.mean(contig[i:i+window])
            #GC = gcCalc(contigDict[key].seq[i:i+window])
            #met = numpy.mean(metDict[key][i:i+window])
            #chip = numpy.mean(chipDict[key][i:i+window])
        else:
            rep=0
        #GC = SeqUtils.GC(str(contigDict[key][i:i+window]))
    
        if contigDict[key].seq[i:i+window]:
            GC = gcCalc(contigDict[key].seq[i:i+window])
        else:
            GC = 0
        if metDict[key][i:i+window]:
            met = numpy.mean(metDict[key][i:i+window])
        else:
            met = 0
        if chipDict[key][i:i+window]:
            chip = numpy.mean(chipDict[key][i:i+window])
        else:
            chip = 0
        
        
        if rep > REP_CUTOFF:
            repGenome.append([key,i,i+window,rep,GC,met,chip])
        else:
            nonrepGenome.append([key,i,i+window,rep,GC,met,chip])
        
        wholeGenome.append([key,i,i+window,rep,GC,met,chip])
    
    #print chipDict[key][i:i+window]
    #print rep, GC, met, chip
    #print wholeGenome[-1]
    #print "\n".join([str(x) for x in wholeGenome])

datatypes = ["repeats", "GC", "methylation", "H3K9me3"]
areas=["rep","non","all"]
datadict = {x: {} for x in datatypes}
coords = [3,4,5,6]
cDict = dict(zip(datatypes,coords))
aDict={"rep": repGenome, "non": nonrepGenome, "all": wholeGenome}

for d in datatypes:
    for c in areas:
        col = getColumn(aDict[c],cDict[d])
        meancol = numpy.mean(col)
        confcol = st.t.interval(0.95, len(col)-1, loc=meancol, scale=st.sem(col))
        datadict[d][c] = [str(meancol), str(confcol[0]), str(confcol[1])]

with open(outfile, 'w') as ofile:
    ofile.write("genome,data_type,repeats_mean,repeats_conf_low,repeats_conf_high,non_repeats_mean,non_repeats_conf_low,non_repeats_conf_high,whole_genome_mean,whole_genome_conf_low,whole_genome_conf_high\n")
    ofile.write(args.filename + ",windows," + str(len(repGenome)) +",,," + str(len(nonrepGenome)) +",,," + str(len(wholeGenome)) +",,\n")
    for dt in datatypes:
        olist = [args.filename,dt]
        for a in areas:
            for dd in datadict[dt][a]:
                olist.append(dd)
        ofile.write(",".join(olist)+"\n")

'''
repMean = numpy.mean(getColumn(repGenome,3))
GCMean = numpy.mean(getColumn(repGenome,4))
metMean = numpy.mean(getColumn(repGenome,5))
chipMean = numpy.mean(getColumn(repGenome,6))

nonrepMean = numpy.mean(getColumn(nonrepGenome,3))
nonGCMean = numpy.mean(getColumn(nonrepGenome,4))
nonmetMean = numpy.mean(getColumn(nonrepGenome,5))
nonchipMean = numpy.mean(getColumn(nonrepGenome,6))

wholerepMean = numpy.mean(getColumn(wholeGenome,3))
wholeGCMean = numpy.mean(getColumn(wholeGenome,4))
wholemetMean = numpy.mean(getColumn(wholeGenome,5))
wholechipMean = numpy.mean(getColumn(wholeGenome,6))


print "Stats."
print "Window length: " + str(window) + ", Step length: " + str(step)
print "# windows: " + str(len(repGenome)) +"\t" + str(len(nonrepGenome)) +"\t" + str(len(wholeGenome))
print "Repeats: " + str(repMean) + "\t" + str(nonrepMean)+ "\t" + str(wholerepMean)
print "GC: " + str(GCMean) + "\t" + str(nonGCMean)+ "\t" + str(wholeGCMean)
print "Methylation: " + str(metMean) + "\t" + str(nonmetMean)+ "\t" + str(wholemetMean)
print "H3K9me: " + str(chipMean) + "\t" + str(nonchipMean)+ "\t" + str(wholechipMean)
'''
