# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:05:18 2016

@author: Anindya Bhattacharya, Department of Computer Science and Engineering, UCSD, 9500 Gilman Drive, La Jolla, CA 92093-0404
Email:anb013@eng.ucsd.edu

Parameters: python proteinE.py Input_PSM_Result_File_Tab_Separated  SampleID_for_Spectra_Files Output_File
Example 1.: python proteinE.py /home/anindya/VU_known_search/VU_known_tumor_FDRspec01.txt  ./Input/colorectalidlist.txt new_out.txt
Example 2.: python proteinE.py /home/anindya/VU_known_search/VU_known_normal_FDRspec01.txt  ./Input/colorectalidlist.txt new_out.txt

Output: File 1. Output_File : Raw PSM count matrix
        File 2. Normalize_Count_Output_File : Normalize PSM count - Protein Spectra per 100 aminoacid per 1000 spectra (PSK)
        File 1. Rank_on_Normalize_Count_Output_File : Rank on PSK. Highest expression get rank 1.

"""


import sys

data=open(sys.argv[1],"r")
samplelist=open(sys.argv[2],"r")
outputkey=sys.argv[3]


simpleinput=outputkey

pskf="Normalize_Count_"+outputkey

rankf="Rank_on_Normalize_Count_"+outputkey



sample_dic={}

for line in samplelist:
    linew=line.strip().split("\t")
    if sample_dic.has_key(linew[0]):
        idl=sample_dic[linew[0]]
        idl=idl+','+linew[1]
        sample_dic[linew[0]]=idl
    else:
        sample_dic[linew[0]]=linew[1]

id_s={}

for sample in sample_dic.values():
    id_s[sample]=1


unif=open("./Hardcoded/gene_id_length.txt","r")

gene_dic={}

for line in unif:
    line=line.strip().split("\t")
    gene=line[0]
    uniprotid=line[1]
    ensemblid=line[2]
    length=int(line[3])
    gene_dic[gene]={'L':length,'E':ensemblid,'U':uniprotid}


ensf=open("./Hardcoded/mart_export.txt","r")

dic_ens={}

for line in ensf:
    line=line.strip().split("\t")
    if len(line)<2:
        continue
    gene=line[0]
    ensid=line[1]
    dic_ens[ensid]=gene


#initialize count dictionary to 0
dic_count={}
id_l={}

for i,line in enumerate(data):
    if i==0:
        continue  

    linew=line.strip().split("\t")

    sampleid=linew[0]
    if sample_dic.has_key(sampleid):
        sample=sample_dic[sampleid]
    else:
        sample=sampleid

    ensid=linew[9]
    if "|" in ensid:
        if "ENSP" not in ensid:
           ensid=ensid.split("|")[1]
        else:
            ew=ensid.split("|")
            for w in ew:
                if "ENSP" in w:
                    ensid=w
    scannumber=linew[2]
    # peptide=linew[8].translate(None, ':-.*!@#$+0123456789')
    if dic_ens.has_key(ensid):
        gene=dic_ens[ensid]
        uid=gene
    else:
        gene=ensid
        uid=ensid
        if not gene_dic.has_key(gene):
            gene_dic[gene]={'L':1,'E':'-','U':'-'}


    if not id_l.has_key(gene):
        id_l[gene]=1

    if dic_count.has_key(uid):
        count=dic_count[uid][sample]
        dic_count[uid][sample]=count+1
    else:
        
        for s in id_s:
            if dic_count.has_key(uid):
                dic_count[uid][s]=0
            else:            
                tmp_dic={}
                tmp_dic[s]=0 
                dic_count[uid]=tmp_dic
        count=dic_count[uid][sample]
        dic_count[uid][sample]=count+1


###############################Print header####################################

outfile=open(simpleinput,"w")

outfilepsk=open(pskf,"w")

outfilerank=open(rankf,"w")

outfile.write("Gene_Name\tUniProtKB_ID\tEnsembl_Protein_ID\t")     

outfilepsk.write("Gene_Name\tUniProtKB_ID\tEnsembl_Protein_ID\t")     

outfilerank.write("Gene_Name\tUniProtKB_ID\tEnsembl_Protein_ID\t")     


for sid in id_s:
    outfile.write("%s\t" %(sid))
    outfilepsk.write("%s\t" %(sid))
    outfilerank.write("%s\t" %(sid))
    
outfile.write("\n")
outfilepsk.write("\n")
outfilerank.write("\n")


#Datastructure for PSK

sum_dic={}
exp_dic={}


for i,sid in enumerate(id_s):
    sum_dic[sid]=0.0
    

##########End psk datastructure

#strart printing Raw count
        
for uid in id_l:
    s=[]
    sum_s=0
    for sample in id_s:
        sum_s+=int(dic_count[uid][sample])
    if sum_s==0:
        continue
    
    if gene_dic.has_key(uid):
        length=gene_dic[uid]['L']
        uniprotid=gene_dic[uid]['U']
        ensemblid=gene_dic[uid]['E']
    else:
        length=1
        uniprotid=''
        ensemblid=''
    outfile.write("%s\t%s\t%s\t" %(uid,uniprotid,ensemblid))
    
    lengthD=float(length)*0.01
    dic_t={}
    sumrow=0.0

    for sample in id_s:
        val=dic_count[uid][sample]
        outfile.write("%d\t" %(val))
        
        #PSK part
        word=float(val)/lengthD
        sum_dic[sample]=sum_dic[sample]+word
        dic_t[sample]=word
    exp_dic[uid]=dic_t 
    #end psk
    #print new line in raw count output        
    outfile.write("\n")

##############################Normalize PSK Count and Rank data##################################################


def rank_simple(vector):
    return sorted(range(len(vector)), key=vector.__getitem__)

def rankdata(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    reversearray=[n-d for d in newarray]
    return reversearray



dic_new_data={}

idlist=[]

for uid in exp_dic:
    idlist.append(uid)
    
    if gene_dic.has_key(uid):
        uniprotid=gene_dic[uid]['U']
        ensemblid=gene_dic[uid]['E']
    else:
        uniprotid=''
        ensemblid=''

    outfilepsk.write("%s\t%s\t%s\t" %(uid,uniprotid,ensemblid))
    for sample in sum_dic:
        sumval=sum_dic[sample]
        psk=(exp_dic[uid][sample]/sumval)*1000
        outfilepsk.write("%f\t" %(psk))
        if dic_new_data.has_key(sample):
            sl=dic_new_data[sample]
            sl.append(psk)
            dic_new_data[sample]=sl
        else:
            sl=[psk]
            dic_new_data[sample]=sl 
    outfilepsk.write("\n")

for sample in dic_new_data:
    rankdatainput=dic_new_data[sample]
    newdata=rankdata(rankdatainput)
    ez=max(newdata)
    dic_new_data[sample]=[0 if x==ez else x for x in newdata]

for i,uid in enumerate(idlist):
    if gene_dic.has_key(uid):
        uniprotid=gene_dic[uid]['U']
        ensemblid=gene_dic[uid]['E']
    else:
        uniprotid=''
        ensemblid=''
    outfilerank.write("%s\t%s\t%s\t" %(uid,uniprotid,ensemblid))
    
    for sample in dic_new_data:
        dl=dic_new_data[sample]
        outfilerank.write("%f\t" %(dl[i]))
    outfilerank.write("\n")


