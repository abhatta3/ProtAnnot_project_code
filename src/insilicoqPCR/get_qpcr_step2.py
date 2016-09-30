import pysam
    
    
data=open("adlVU_qpcr_tumor.txt","r")

mapfile=open("adlVU_sample_encoded_list.txt","r")

outfile=open("testQPCR_report_colorectal.txt","w")

file_dic={}
for n,line in enumerate(mapfile):
    line=line.strip().split('\t')
    idx=int(line[1])
    file_dic[idx]=line[0]


print(n)

outfile.write("Event\tPeptide\tMutation_chromosome\tMutation_Location\tMutation_allele\tQPCR_Datatype\t")
for i in range(1,n+2):
    outfile.write("Sample_%d\t" %(i))
outfile.write("\n")

for line in data:
    line=line.strip().split('\t')
    t_dic={}
    n_dic={}
    for f in file_dic:
        t_dic[f]=0
        n_dic[f]=0

    chrom=line[2]
    loc1=int(line[3])
    loc2=loc1+1
    samplelist=line[5].split(',')
    
    for s in samplelist:
        s=s.split(':')
        idx=int(s[0])
        val=int(s[1])
        t_dic[idx]=val
    
    
    for f in file_dic:
        fname=file_dic[f]
        fnamebam="./VU/bam/"+fname.replace("_reorder_sortsam_realigned.vcf",".bam")
        
        samfile = pysam.AlignmentFile(fnamebam, "rb")
        iter = samfile.fetch(chrom, loc1, loc2)
        x=[i for i,x in enumerate(iter)]
        n_dic[f]=i+1
    

    for l in line[0:-1]:
        outfile.write("%s\t" %(l))
    outfile.write("Reference_Count\t")    
    for i in range(1,n+2):
        if t_dic.has_key(i):
            tval=t_dic[i]
        else:
            tval=0
        if n_dic.has_key(i):
            nval=n_dic[i]
        else:
            nval=0
        dif=nval-tval
        outfile.write("%d\t" %(dif))

    outfile.write("\n")



    for l in line[0:-1]:
        outfile.write("%s\t" %(l))
    outfile.write("Mutation_Count\t")
    for i in range(1,n+2):
        if t_dic.has_key(i):
            tval=t_dic[i]
        else:
            tval=0
        outfile.write("%d\t" %(tval))

    outfile.write("\n")
    
    
