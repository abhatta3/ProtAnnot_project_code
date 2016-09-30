import math as m
import sys
#inf=open("ovarianids","r")

#psknormal=open("PSK_Tumor_ovarian_PNNL_knownDB_FDRspec01_Protein_Expression.txt","r")

#ranknormal=open("RANK_Tumor_ovarian_PNNL_knownDB_FDRspec01_Protein_Expression.txt","r")


#countnormal=open("Tumor_ovarian_PNNL_knownDB_FDRspec01_Protein_Expression.txt","r")


#psktumor=open("PSK_Tumor_ovarian_JHU_knownDB_FDRspec01_Protein_Expression.txt","r")

#ranktumor=open("RANK_Tumor_ovarian_JHU_knownDB_FDRspec01_Protein_Expression.txt","r")

#counttumor=open("Tumor_ovarian_JHU_knownDB_FDRspec01_Protein_Expression.txt","r")
#out=open("combined_ovarian.txt","w")


#####################################################
'''
inf=open("selectedid_colorectal.txt","r")

psknormal=open("colorectal/Normalize_Count_colorectal_normal_count.txt","r")
ranknormal=open("colorectal/Rank_on_Normalize_Count_colorectal_normal_count.txt","r")
countnormal=open("colorectal/colorectal_normal_count.txt","r")


psktumor=open("colorectal/Normalize_colorectal_tumor_count.txt","r")
ranktumor=open("colorectal/Rank_on_Normalize_Count_colorectal_tumor.txt","r")
counttumor=open("colorectal/colorectal_tumor_count.txt","r")

out=open("combined_colorectal_beta_actin_control.txt","w")

'''

def get_dic(filen,start):
    with open(filen,"r") as filep:
        t_dic={}
        for i,line in enumerate(filep):
            line=line.strip().split('\t')
            if i==0:
                continue
            t_dic[line[0]]=line[start:]
    return t_dic


def get_reciprocal_of_mean_reciprocal_rank_score(filen,start):
    with open(filen,"r") as filep:
        t_dic={}
        for i,line in enumerate(filep):
            line=line.strip().split('\t')
            if i==0:
                continue
            rcv=0.0
            size=len(line[start:])
            for word in line[start:]:
                word=float(word)
                if word==0:
                    size=size-1
                else:
                    rcv+=1.0/word
            if size<=0:
                rcv=0
            else:
                rcv=rcv/size
                rcv=1/rcv
            t_dic[line[0]]=rcv
    return t_dic

count_data=[]
rank_data=[]
for i,item in enumerate(sys.argv[1:]):
    if item=='-C' or item=='-c':
        count_data.append(get_dic(sys.argv[i+2],3))
    
    if item=='-R' or item=='-r':
        rank_data.append(get_reciprocal_of_mean_reciprocal_rank_score(sys.argv[i+2],3))
        
    if  item=='-O' or item=='-o':
        outputfile=sys.argv[i+2]
        
        

rank_dic={}

for item_dic in rank_data:
    for k in item_dic:
        val=str(item_dic[k])+','+str(item_dic['EGFR'])
        if rank_dic.has_key(k):
            rank_dic[k]+=';'+val
        else:
            rank_dic[k]=val


count_dic={}

for item_dic in count_data:
    for k in item_dic:
        val=item_dic[k]
        valegfr=item_dic['EGFR']
        lv=''
        slv=0
        tc=0
        for i,v in enumerate(val):
            ev=float(valegfr[i])
            v=float(v)
            if ev!=0 and v!=0:
                vd=m.log(v/ev,2)
                if vd!=0:
                    tc+=1
                    slv+=vd
                lv+=str(vd)+','
        if tc>0:
            avglog=slv/tc
        else:
            avglog=0
            
            
        if count_dic.has_key(k):
            count_dic[k]['val']+=';'+str(avglog)
            count_dic[k]['list']+=';'+lv
            
        else:
            count_dic[k]={'val':str(avglog),'list':lv}
            
            
with open(outputfile,"w") as outf:
                
    for gene in count_dic:
        outf.write("%s\t%s\t%s\t" %(gene,count_dic[gene]['val'],count_dic[gene]['list']))
        if rank_dic.has_key(gene):
            rankvalue=rank_dic[gene]
        else:
            rankvalue='-'
        outf.write("%s\n" %(rankvalue))
