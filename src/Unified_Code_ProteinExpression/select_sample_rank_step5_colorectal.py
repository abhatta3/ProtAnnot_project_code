
normalrank=open("colorectal/Rank_on_Normalize_Count_colorectal_normal_count.txt","r")

tumorrank=open("colorectal/Rank_on_Normalize_Count_colorectal_tumor.txt","r")

plist=open("selectedid_colorectal.txt","r")

outfile=open("Colorectal_expression_rank_score_beta_actin_control.txt","w")



protein=[]
for line in plist:
    line=line.strip().split("\t")
    protein.append(line[0])





def get_sample_name_list(rankfilepointer):
    nlist=[]
    for i,line in enumerate(rankfilepointer):
        line=line.strip().split("\t")
        if i==0:
            nlist=line[3:]
            break
    return nlist
    
    
    
headlistnormal=get_sample_name_list(normalrank)    
headlisttumor=get_sample_name_list(tumorrank)    


def get_reciprocal_of_mean_reciprocal_rank_score(rankfilepointer,protein):
    protein_dic={}
    for i,line in enumerate(rankfilepointer):
        line=line.strip().split("\t")
        if i==0:
            continue

        if line[0] in protein:
            rcv=0.0
            size=len(line[3:])
            for word in line[3:]:
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
            protein_dic[line[0]]={'list':line[3:],'score':rcv}
    return protein_dic



protein_dic_normal=get_reciprocal_of_mean_reciprocal_rank_score(normalrank,protein)
protein_dic_tumor=get_reciprocal_of_mean_reciprocal_rank_score(tumorrank,protein)
    



outfile.write("UNI_ID\t")
for sample in headlistnormal:
    outfile.write("%s\t" %(sample)) 
outfile.write("Average_reciprocal_rank\t")

for sample in headlisttumor:
    outfile.write("%s\t" %(sample)) 
outfile.write("Average_reciprocal_rank_tumor_sample\tUP_Regulation_in_Tumor\n")
normalempty=["NA"]*30
tumorempty=["NA"]*95

for prot in protein:
    outfile.write("%s\t" %(prot))
    if protein_dic_normal.has_key(prot):
        normallist=protein_dic_normal[prot]['list']
        nscore=protein_dic_normal[prot]['score']
    else:
        normallist=normalempty
        nscore=0.0

    for rank in normallist:
        outfile.write("%s\t" %(rank)) 
    outfile.write("%f\t" %(nscore))


    if protein_dic_tumor.has_key(prot): 
        tumorlist=protein_dic_tumor[prot]['list']
        tscore=protein_dic_tumor[prot]['score']
    else: 
        tumorlist=tumorempty
        tscore=0.0

    for rank in tumorlist:
        outfile.write("%s\t" %(rank))
    if tscore>0.0: 
        outfile.write("%f\t%f\n" %(tscore,nscore/tscore))
    else:
        outfile.write("%f\t%f\n" %(tscore,tscore))








