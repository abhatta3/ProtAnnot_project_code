inputlf=open("additional_select_colorectal","r")
#"selected_colorectal.txt"
#inputlf=open("samegene_colorectal.txt","r")

inputevents=open("CPTAC_Colon_Proteogenomic_Events.txt","r")
data=open("VUvcf_correct.spl","r" )

output=open("adlVU_qpcr_tumor.txt","w")
outputsample=open("adlVU_sample_encoded_list.txt","w")

#output=open("VU_qpcr_tumor_same_gene.txt","w")
#outputsample=open("VU_sample_encoded_list_same_gene.txt","w")

event_d={}

for line in inputlf:
    line=line.strip().split("\t")
    event_d[line[0]]=1


mut_dic={}
sample_dic={}
for i,line in enumerate(data):
    if i==0:
        line=line.replace('#', '')
        line=line.strip().split(",")
        for l in line:
            l=l.split(':')
            sample_dic[l[1]]=l[0]
    if i<3:
        continue
        
  
    line=line.strip().split("\t")
    
    if len(line)==5:
        chrom=line[0]
        loc=line[1]
        mut=line[2]
        sample=line[3]
        if mut_dic.has_key(chrom):
            mut_dic[chrom][loc]={'mut':mut,'sample':sample}
        else:
            loc_dic={}
            loc_dic[loc]={'mut':mut,'sample':sample}
            mut_dic[chrom]=loc_dic
        

for k in sample_dic:
    outputsample.write("%s\t%s\n" %(sample_dic[k],k))

pep_dic={}

for line in inputevents:
    line=line.strip().split("\t")
    if event_d.has_key(line[0]):
        event=line[0]
        peptide=line[2]
        chrom=line[3]
        tl=line[4].split(";")
        loc=[]
        for itl in tl:
            itl=itl.split("-")
            if len(itl)<2:
                continue
            else:
                loc.append(itl[1])  
  
        if mut_dic.has_key(chrom):
            locdic=mut_dic[chrom]
            for l in loc:
                if locdic.has_key(l):
                    mut=locdic[l]['mut']
                    sample=locdic[l]['sample']
                    if sample_dic.has_key(sample):
                        sname=sample_dic[sample]
                    else:
                        sname='NA'
                    output.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(event,peptide,chrom,l,mut,sample))
                    
        
        
