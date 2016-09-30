import os
import re
import protein_annot as pa
import tabix


def siliding_window_score(newseq,mut_pos_start,score_dic,key='hw'):
    mut_hw_string=''
    mut_hw=0
    for ci,chl in enumerate(newseq):
        mut_hw_string+=str(score_dic[chl][key])+','
        if mut_pos_start==ci:
            mut_hw=score_dic[chl][key]

    value=(mut_hw,mut_hw_string)
    return value
    
    

    
def petide_file_scan(inputfile,outputfile,outputfilename,dir_path,annotf,expfile,qpcrfile):
        
    os.mkdir(dir_path+outputfilename+".dir", 0755 );
    path=dir_path+outputfilename+".dir"
    annotfile=open(dir_path+"Datafiles/uniprot.seq.tab","r")
    protter_head="http://wlab.ethz.ch/protter/#"
    protter_tail_1="&tm=auto&mc=lightsalmon&lc=blue&tml=numcount&numbers&legend&n:exp.peps,fc:darkblue=EX.PEPTIDES&n:signal peptide,fc:red,bc:red=UP.SIGNAL&n:disulfide bonds,s:box,fc:greenyellow,bc:greenyellow=UP.DISULFID&n:variants,s:diamond,fc:orange,bc:orange=UP.VARIANT,"
    protter_tail_2="&n:PTMs,s:box,fc:forestgreen,bc:forestgreen=UP.CARBOHYD,UP.MOD_RES&n:exp.peps,cc:white,bc:blue=EX.PEPTIDES&format=svg"
    
    #url="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz"
    #tb = tabix.open(url)
    
    dic_qpcr={}
    
    if qpcrfile!='':
        with open(qpcrfile,"r") as qpf:
            for line in qpf:
                line=line.strip().split('\t')
                dic_qpcr[line[0]]=line[1]


    
    dic_exp_annot={}    
    if expfile!='':
        with open(expfile,"r") as expanf:
            for line in expanf:
                line=line.strip().split('\t')
                dic_exp_annot[line[0]]={'count':line[1],'cont_seq':line[2],'rank':line[3]}


    
    dic_pep_annot={}    
    if annotf!='':
        with open(annotf,"r") as anf:
            for line in anf:
                if "#" not in line:
                    line=line.strip().split('\t')
                    if len(line)<4:
                        continue
                    dic_pep_annot[line[0]]={'specno':len(line[1].split(',')),'rnano':len(line[2].split(',')),'rnad':int(line[3])}

    annotfilescore=open(dir_path+"Datafiles/Score_table.txt","r")
    
    score_dic={}

    for i,line in enumerate(annotfilescore):
        if i<2:
            continue
        line=line.strip().split('\t')
        score_dic[line[0]]={'kd':float(line[1]),'hw':float(line[2]),'charge':int(line[3]),'sa':int(line[4])}
        


    #tmlist={'mutation':'','insertion':'','deletion':'','alternative splice':'','novel exon':''}
    #GOCC=pa.getanotCC(dir_path)


    dic_seq={}
    for line in annotfile:
        line=line.strip().split("\t")
        gene=line[0]
        dseq=line[2]
        if not dic_seq.has_key(gene):
            d_u={}
            d_u[line[1]]=dseq
            dic_seq[gene]=d_u
        else:
            dic_seq[gene][line[1]]=dseq


    UP_dic={} # Dictionary to hold data from uniprot
    
    for i,line in enumerate(inputfile):
        line=line.strip().split('\t')
        term=''
        if i==0:
            for l in line[0:3]:
                outputfile.write("%s\t" %(l))
            outputfile.write("Reference_Peptide\tRef_aa\tMut_aa\tRef_Location\tMut_loc_at_peptide\tTotal_number_of_overlapping_gene\tGene\tProtein\tCellular component\tProtein_position_in_protter_plot\tProtein_Plot\tProteinAtlas\tGeneCards")
            for l in line[3:9]:
                outputfile.write("\t%s" %(l))
            outputfile.write("\tWT_Charge\tMT_Charge\tCharge_Diff\tWT_Surface\tMT_Surface\tSurface_Diff\tKyte-Doolittle_hydrophilic(>0)\tHopp-Woods_hydrophilic(>0)\tdbSNP_Overlap\tSpectra_Samples\tRNA_Samples\tRNA_depth\tMutation_Event_Position\tlog2(Gene/EGFR_Protein_Expression\tlog2_ratios\tGene,EGFR_Rank_Protein_Expression\tInsilico_qPCR;qPCR_Values\n")

            continue

        pep=line[2]
        
        jointpattern=pa.getsearchpeptidepattern(pep)
        
        pattern_left=pa.getsearchpeptidepattern_l(pep)

        pattern_right=pa.getsearchpeptidepattern_r(pep)
       
        
        if len(line)<15:
            glist="NA"
        else:
            glist=line[14]
        if glist=="":
            glist="NA"
        gene=glist.strip().split(';')
        numberofgene=len(gene)-1
        
        
        NewPeptide=pep.translate(None, '-.*!@#$?:')

        for g in gene[0:-1]:
            Ref_Ppeptide=""
            seq=""
            upid=""
            if dic_seq.has_key(g):
                u_d=dic_seq[g]
                for upid in u_d:
                    seq=u_d[upid]
                    lookforpep = re.findall('(?=('+jointpattern+'))',seq)
                    if len(pattern_left)>6:
                        lookforpep_left = re.findall('(?=('+pattern_left+'))',seq)
                        if len(lookforpep_left)>1:
                            lookforpep_left=''
                    else:
                        lookforpep_left=''
                    if len(pattern_right)>6:
                        lookforpep_right = re.findall('(?=('+pattern_right+'))',seq)
                        if len(lookforpep_right)>1:
                            lookforpep_right=''
                    else:
                        lookforpep_right=''
                        
                    if len(pattern_left)>len(pattern_right):
                        lookfor_second=lookforpep_left
                        lookfor_third=lookforpep_right
                    else:
                        lookfor_second=lookforpep_right
                        lookfor_third=lookforpep_left
                    
                    if lookforpep:
                        Ref_Ppeptide = min(lookforpep, key =len)
                        result=pa.compare(Ref_Ppeptide,NewPeptide)
                        locmut=result['loc']
                        raa=result['raa']
                        maa=result['maa']
                    elif lookfor_second and lookfor_second!='':
                        
                        Ref_Ppeptide = min(lookfor_second, key =len)
                        result=pa.compare(Ref_Ppeptide,NewPeptide)
                        locmut=result['loc']
                        raa=result['raa']
                        maa=result['maa']
                    elif lookfor_third and lookfor_third!='':
                        Ref_Ppeptide = min(lookfor_third, key =len)
                        result=pa.compare(Ref_Ppeptide,NewPeptide)
                        locmut=result['loc']
                        raa=result['raa']
                        maa=result['maa']
                    else:
                        continue
                    if line[1]=="mutation":
                        if len(raa)!=1:
                            continue
                        if len(maa)!=1:
                            continue
                    
                    if raa=='':
                        continue
                    if len(Ref_Ppeptide)<=3:
                        continue

                    offset=seq.find(Ref_Ppeptide)+1
                    loc=locmut+offset
                    locend=loc+len(raa)-1
                
                    protter_link=protter_head+'up='+upid+'&peptides='+Ref_Ppeptide+protter_tail_1+str(loc)+'-'+str(locend)+protter_tail_2
                    
                    
                    protter_for_image='http://wlab.ethz.ch/protter/create?'+'up='+upid+'&peptides='+Ref_Ppeptide+protter_tail_1+str(loc)+'-'+str(locend)+protter_tail_2

                
                
                    genecard="http://www.genecards.org/cgi-bin/carddisp.pl?gene="+g
                    proteinatlas="http://www.proteinatlas.org/"+upid+"-"+g+"/cancer"
                    
                    if UP_dic.has_key(upid):
                        localization=UP_dic[upid]['K']
                        term=UP_dic[upid]['C']
                        ED=UP_dic[upid]['E']
                        VD=UP_dic[upid]['V']
                    else:
                        UP_dic=pa.get_uniprot_dic(UP_dic,upid)
                        localization=UP_dic[upid]['K']
                        term=UP_dic[upid]['C']
                        ED=UP_dic[upid]['E']
                        VD=UP_dic[upid]['V']
                        
                    inex=0
                    if "Membrane-Bound" in term:
                        for edic in ED:
                            start_R=ED[edic]['start']
                            end_R=ED[edic]['end']
                            if loc>=start_R and loc<=end_R:
                                inex=2
                                break
                            if locend>=start_R and locend<=end_R:
                                inex=2
                                break
                    if inex==2:
                        peptype='Targetable'
                    else:
                        peptype='Unknown'

                    if inex==0 and "Membrane-Bound" in term:
                        peptype=pa.download_prottersvg(protter_for_image,upid,path)

                    dbsnpstr=""    
                    for vdic in VD:
                        start_v=VD[vdic]['start']
                        des_v=VD[vdic]['des']
                        
                        if start_v>=loc and start_v<=locend:
                            dbsnpstr=dbsnpstr+"#"+des_v

                    for l in line[0:3]:
                        outputfile.write("%s\t" %(l))
                    outputfile.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(Ref_Ppeptide,raa,maa,loc,locmut,numberofgene,g,upid,localization,term,protter_link,proteinatlas,genecard))
                    for l in line[3:9]:
                        outputfile.write("\t%s" %(l))

                    if dic_pep_annot.has_key(pep):
                        specno=dic_pep_annot[pep]['specno']
                        rnano=dic_pep_annot[pep]['rnano']
                        rnad=dic_pep_annot[pep]['rnad']           
                    else:
                        specno=0
                        rnano=0
                        rnad=0   
                            
                    ref_charge=0        
                    for chl in raa:
                        ref_charge+=score_dic[chl]['charge']
                        
                    mut_charge=0        
                    for chl in maa:
                        mut_charge+=score_dic[chl]['charge']
                        
                    charge_dif=int(mut_charge)-int(ref_charge)
                        
                    ref_sa=0        
                    for chl in raa:
                        ref_sa+=score_dic[chl]['sa']
                        
                    mut_sa=0        
                    for chl in maa:
                        mut_sa+=score_dic[chl]['sa']
            
                    sa_dif=mut_sa-ref_sa
            
                    
                    mut_hw,mut_hw_string=siliding_window_score(NewPeptide,locmut,score_dic,key='hw')
                    mut_kd,mut_kd_string=siliding_window_score(NewPeptide,locmut,score_dic,key='kd')
                            
                    
                    outputfile.write("\t%d\t%d\t%d\t%d\t%d\t%d\t%f;%s\t%f;%s\t%s\t%d\t%d\t%d\t%s\t" %(ref_charge,mut_charge,charge_dif,ref_sa, mut_sa, sa_dif, mut_kd,mut_kd_string, mut_hw, mut_hw_string, dbsnpstr,specno, rnano,rnad,peptype))


                    if dic_exp_annot.has_key(g):
                        ecount=dic_exp_annot[g]['count']
                        econt_seq=dic_exp_annot[g]['cont_seq']
                        erank=dic_exp_annot[g]['rank']           
                    else:
                        ecount='-'
                        econt_seq='-'
                        erank='-'     
                    if dic_qpcr.has_key(line[0]):
                        qpcrvalue=dic_qpcr[line[0]]
                    else:
                        qpcrvalue='0;'
                        
                    outputfile.write("%s\t%s\t%s\t%s\n" %(ecount,econt_seq,erank,qpcrvalue))

    os.rmdir(path)
    outputfile.close()
    inputfile.close()


   
    

def petide_file_filter(inputfile,outputfile,outputfilename,dir_path,annotf):
    url="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz"
    tb = tabix.open(url)
    if annotf!='':
        dic_pep_annot={}
        with open(annotf,"r") as anf:
            for line in anf:
                if "#" not in line:
                    line=line.strip().split('\t')
                    if len(line)<4:
                        continue
                    dic_pep_annot[line[0]]={'specno':len(line[1].split(',')),'rnano':len(line[2].split(',')),'rnad':int(line[3])}

    annotfile=open(dir_path+"Datafiles/Score_table.txt","r")
    
    score_dic={}

    for i,line in enumerate(annotfile):
        if i<2:
            continue
        line=line.strip().split('\t')
        score_dic[line[0]]={'kd':float(line[1]),'hw':float(line[2]),'charge':int(line[3]),'sa':int(line[4])}
    
    for i,line in enumerate(inputfile):
        line=line.strip().split('\t')
        if i==0:
            for j,l in enumerate(line):
                if j>21:
                    break
                outputfile.write("%s\t" %(l))
            outputfile.write("WT_Charge\tMT_Charge\tCharge_Diff\tWT_Surface\tMT_Surface\tSurface_Diff\tKyte-Doolittle_Hydrophobicity(<0 = hydrophilic)\tHopp-Woods_Hydrophobicity(>0 = hydrophilic)\tdbSNP_Overlap\tSpectra_Samples\tRNA_Samples\tRNA_depth\n")

            continue

        ref=line[4]
        mut=line[5]
        protter=line[12]
        
        pep=line[2]
        
        if dic_pep_annot.has_key(pep):
            specno=dic_pep_annot[pep]['specno']
            rnano=dic_pep_annot[pep]['rnano']
            rnad=dic_pep_annot[pep]['rnad']           
        else:
            specno=0
            rnano=0
            rnad=0           
            
        if protter!="Membrane-Bound" and protter!="Unknown":
            continue

        
        chrom=line[3].translate(None, 'chrC')
        
        loc_string=line[4]
        dbsnpstr=""

        if "M1" in loc_string:
            loc_string=loc_string.split(';')
            s=len(loc_string)
            for k,item in enumerate(loc_string):
                if item=='M1' and k>0 and k+1<s:
                    start=loc_string[k-1].split('-')[1]
                    end=loc_string[k+1].split('-')[0]
                    try:
                        records = tb.query(chrom, int(start), int(end))
                        for record in records:
                            dbsnpstr=dbsnpstr+"#"+record[7]

                    except Exception as e:
                        print(e)
                    break
                
        ref_charge=0        
        for chl in ref:
            ref_charge+=score_dic[chl]['charge']
            
        mut_charge=0        
        for chl in mut:
            mut_charge+=score_dic[chl]['charge']
            
        charge_dif=int(mut_charge)-int(ref_charge)
            
        ref_sa=0        
        for chl in ref:
            ref_sa+=score_dic[chl]['sa']
            
        mut_sa=0        
        for chl in mut:
            mut_sa+=score_dic[chl]['sa']

        sa_dif=mut_sa-ref_sa

        mut_hw=0        
        for chl in mut:
            mut_hw+=score_dic[chl]['hw']
            

        mut_kd=0        
        for chl in mut:
            mut_kd+=score_dic[chl]['kd']
            
        
  
        for j,l in enumerate(line):
            if j>21:
                break
            outputfile.write("%s\t" %(l))
        
        outputfile.write("%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%d\t%d\t%d\n" %(ref_charge,mut_charge,charge_dif,ref_sa, mut_sa, sa_dif, mut_kd, mut_hw, dbsnpstr,specno, rnano,rnad))
                       


    outputfile.close()
    annotfile.close()
    inputfile.close()

    
    