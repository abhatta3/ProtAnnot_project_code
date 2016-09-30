import sys,os 

try:
    import getopt
except ImportError:
    print('\nThere is no such module getopt installed')
    sys.exit(1)

try:
    import urllib
except ImportError:
    print('\nThere is no such module urllib installed')
    sys.exit(1)




    
    
def getfilename_w_annot(argv):
	filename = {}
	try:
		opts, args = getopt.getopt(argv,"hi:t:a:o:e:p:",["ifile=","taskid=","afile=","ofile=","expfile=","qpcrfile="])
	except getopt.GetoptError:
		print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'
		
	for opt, arg in opts:
		if opt == '-h':
        		print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'
			sys.exit(2)
		elif opt == '-i':
			filename["inputfile"]=arg
		elif opt == '-o':
			filename["outputfile"]=arg
		elif opt == '-t':
			filename["taskid"]=arg
		elif opt == '-a':
			filename["afile"]=arg   
		elif opt == '-e':
			filename["expfile"]=arg   
		elif opt == '-p':
			filename["qpcrfile"]=arg   
   
		else:
        		print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'
			sys.exit(2)
	return filename

    

def getfilename(argv):
	filename = {}
	try:
		opts, args = getopt.getopt(argv,"hi:t:o:",["ifile=","taskid=","ofile="])
	except getopt.GetoptError:
		print 'python pep_reporter.py -i <inputfile> -t <task_id_from_ENOSI_result> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python pep_reporter.py -i <inputfile> -t <task_id_from_ENOSI_result> -o <outputfile>'
			sys.exit(2)
		elif opt == '-i':
			filename["inputfile"]=arg
		elif opt == '-o':
			filename["outputfile"]=arg
		elif opt == '-t':
			filename["taskid"]=arg
		else:
			print 'python pep_reporter.py -i <inputfile> -t <task_id_from_ENOSI_result> -o <outputfile>'
			sys.exit(2)
	return filename


def check_protein_in_membrane(fname):
    
    term='-'
    upperpep=-1
    lowerpep=-1
    external=''
    nterm='-'
    if os.path.isfile(fname):
        with open(fname, 'r') as r_svg:
            input_SVG_lines = r_svg.readlines()
            blue_l=[]
            all_c=[]
            for l,line in enumerate(input_SVG_lines):
                if "legend" in line:
                    break
                if "<circle id=" in line and "#AAAAAA" not in line:
                    all_c.append(line)
                if "#0000ff" in line and "<circle id=" in line and "#AAAAAA" not in line:
                    blue_l.append(line)
                if nterm=='-' and "<text x=" in line:
                    nterm=line.split('>')[1].split('<')[0]
                    
                if "extra" in line or "Extracellular" in line:
                    if "<text x=" in line:
                        external=line
            if external=='':
                return(nterm)
                
            ewords=external.split("\'")
            
            for i,ew in enumerate(ewords):
                if "y=" in ew:
                    eb=float(ewords[i+1])
                    break
            #check for any peptide under extra to be membrane bound

            for item in all_c:
                words=item.split("\'")
                for i,w in enumerate(words):
                    if "cy=" in w:
                        pepb=float(words[i+1])
                        if eb<pepb:
                            lowerpep=1
                            break
                if lowerpep==1:
                    break 
            #check for blue node at top of extra  
            for item in blue_l:
                words=item.split("\'")
                for i,w in enumerate(words):
                    if "cy=" in w:
                        pepb=float(words[i+1])
                        if eb>pepb:
                            upperpep=1
                            break
                if upperpep==1:
                    break
                
    #printing protein class
    if lowerpep!=1 and upperpep==1:
        term='Secreted-Peptide'
    elif lowerpep==1 and upperpep==1:
        term='Targetable-Peptide'
    elif lowerpep==1 and upperpep!=1:
        term='Intra-Cellular-Peptide'
        
    return(term)
    

def compare(Ref_Ppeptide,NewPeptide):
    result={'loc':0,'raa':'','maa':''}
    rstarti=0
    rendi=len(Ref_Ppeptide)-1
    nstarti=0
    nendi=len(NewPeptide)-1
    p1=0
    p2=0
    for xi,x in enumerate(NewPeptide):
        if xi>len(Ref_Ppeptide)-1:
            break
        if rendi-xi<0:
            break
        if x != Ref_Ppeptide[xi] and p1==0:
            result['loc']=xi
            rstarti=xi
            nstarti=xi
            p1=1
        if NewPeptide[nendi-xi] != Ref_Ppeptide[rendi-xi] and p2==0:
            nendi=nendi-xi
            rendi=rendi-xi             
            p2=1
        if p1==1 and p2==1:
            if rstarti<rendi:
                result['raa']=Ref_Ppeptide[rstarti:rendi+1]
            else:
                result['raa']=Ref_Ppeptide[rendi:rstarti+1]
 
            if nstarti<nendi:
                result['maa']=NewPeptide[nstarti:nendi+1]
            else:
                result['maa']=NewPeptide[nendi:nstarti+1]
            break    
        
    return result



def compare_old(Ref_Ppeptide,NewPeptide):
    result={'loc':0,'raa':'','maa':''}
    for xi,x in enumerate(NewPeptide):
        if xi>=len(Ref_Ppeptide):
            break
        if x != Ref_Ppeptide[xi]:
            result['loc']=xi
            result['raa']=Ref_Ppeptide[xi]
            result['maa']=x
            break    
    return result


def getsearchpeptidepattern(Peptide):
    Peptide=Peptide.translate(None, '-*!@#$?')
    allpart=Peptide.split(".")
    Peptide=allpart[1]

    while ":" in Peptide:
        p=Peptide.index(":")
        if p > 0 and p+1<len(Peptide):
            newPeptide=Peptide[0:p-1]+".+?"+Peptide[p+2:]
            Peptide=newPeptide
        elif (p+1)==len(Peptide):
            newPeptide=Peptide[0:p-1]+".+?"
            Peptide=newPeptide
            break
        elif p==0:
            newPeptide=".+?"+Peptide[p+2:]
            Peptide=newPeptide
    rpeptide=allpart[0]+Peptide+allpart[2]
    
    return rpeptide
 


def getsearchpeptidepattern_l(Peptide):
    Peptide=Peptide.translate(None, '-*!@#$?')
    if Peptide.find('.')!=-1:
        allpart=Peptide.split(".")
        Peptide=allpart[1]
    
    
    if Peptide.find(':')==-1:
        Peptide=Peptide[0:-1]+".+?"
    else:
        p=Peptide.index(":")
        Peptide=Peptide.split(':')[0]
        Peptide=Peptide[0:p-1]+".+?"
    rpeptide=allpart[0]+Peptide
        
    return rpeptide


def getsearchpeptidepattern_r(Peptide):
    Peptide=Peptide.translate(None, '-*!@#$?')
    if Peptide.find('.')!=-1:
        allpart=Peptide.split(".")
        Peptide=allpart[1]
        
    if Peptide.find(':')==-1:
        Peptide=".+?"+Peptide[1:]
    else:
        p=Peptide.index(":")
        Peptideall=Peptide.split(':')
        Peptide=Peptideall[1]
        Peptide=".+?"+Peptide[p+1:]
        if len(Peptideall)>2:
            for item in Peptideall[2:]:
                if len(item)>1:
                    Peptide=Peptide[0:-1]+".+?"+item[1:]
                else:
                    Peptide=Peptide[0:-1]+".+?"
                    
    rpeptide=Peptide+allpart[2]
        
    return rpeptide




def download_prottersvg(url,fname,path):
    """
    download a protter plot
    """
    os.chdir(path)
    term=''
    try:
        image=urllib.URLopener()
        image.retrieve(url,fname)
        term=check_protein_in_membrane(fname)
        os.remove(fname)
    except Exception as e:
        print e
        
    os.chdir('../')
    if term=='' or term=='-':
        term='Unknown-Pep-localization'
    
    return(term)


def get_uniprot_dic(UP_dic,UPID):
    kword=""
    transmem=0
    extra=0
    vc=0
    extracell={}
    variant={}
    clt='Unknown'
    try:
        data = urllib.urlopen("http://www.uniprot.org/uniprot/"+UPID+".txt").read()
        for line in data.split('\n'):
            tline=line
            line=line.split()
            if len(line)>=4:
                if "FT" == line[0] or "KW" == line[0]:
                    if "KW" == line[0]:
                        kword+="-".join(line[1:])
                    if "FT" == line[0]:
                        if "Extracellular" in tline:
                            extracell[extra]={'start':int(line[2]),'end':int(line[3])}
                            extra+=1
                        if "TRANSMEM" in line[1] and transmem==0:
                            transmem=1
                        if "VARIANT" in line[1]:
                            variant[vc]={'start':int(line[2]),'end':int(line[3]),'des':"".join(line[4:])}
                            vc+=1
                            
        
        if transmem>0:
            clt='Membrane-Bound'
        elif extra>0 and transmem==0:
            clt='Secreted'
        else:
            clt='Intracellular'
    except Exception as e:
        print(e)

    UP_dic[UPID]={'K':kword,'C':clt,'E':extracell,'V':variant}
        
    return UP_dic

    
        
        
        

