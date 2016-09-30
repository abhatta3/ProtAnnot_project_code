# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 14:39:29 2016

@author: anindya
"""


import sys 

import os

try:
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
except ImportError:
    print('\nThere is no such module matplotlib installed')
    sys.exit(1)

from build_css_html_table_v5_new import CSSHTMLBuildFilter

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'


#if __name__ == "__main__":


try:
    inputfile=sys.argv[1]
    taskid=sys.argv[2]
    

except Exception as e:  
	print 'python Selection_of_membrane_bound.py -i <inputfile> -t <task_id_from_ENOSI_result>'
	sys.exit(2)





print("Input parameters are in correct format...Computation is in progress....")

filteredfile=inputfile+"_filtered"

eventfile=inputfile+"_filtered_event"

genefile=inputfile+"_filtered_gene"


event_dic={}
gene_dic={}



    


def create_figure_file(y,pdf,title='',label=[]):
    
    N = len(y)
    x = range(N)
    width = 1/1.5
    plt.bar(x, y, width, color="blue")
    if len(label)==len(x):
        plt.xticks(x, label)
    
    plt.title(title)
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


        
with open(inputfile,'r') as data, open(filteredfile,'w') as filtered_file, open(eventfile,"w") as eventf, open(genefile,"w") as genef, PdfPages(filteredfile+'_plot_Results.pdf') as pdf :
    
    for i,item in enumerate(data):
        line=item.strip().split('\t')
        
        if i==0:
            for word in line:
                filtered_file.write("%s\t" %(word))
            filtered_file.write("\n")
            continue
        
        if 'Targetable' not in item:
            continue
        
        if 'IG gene' in item:
            continue
        
        if not event_dic.has_key(line[0]):
            event_dic[line[0]]=1
        if not gene_dic.has_key(line[9]):
            gene_dic[line[9]]=1
            
        try:
            peptide=line[2].translate(None, '-.*!@#$?:')
            loc=int(line[7])
            label=[]
            for k,pw in enumerate(peptide):
                if k==loc:
                    pw=pw+'*'
                label.append(pw)
                
            if ';' in line[29] and ',' in line[29]:
                dt=line[29].split(';')[1]
                dt=dt.split(',')
                y1=[]
                for y in dt:
                    if '.' in y:
                        y=float(y)
                        y1.append(y)
                if len(y1)>0:
                    create_figure_file(y1,pdf,title='Event-'+line[0]+'-Hydrophilic',label=label)
                    
            if ',' in line[36]:
                alls=line[36].split(';')
                lsize=len(alls)
                if lsize>1:
                    sname=['Tumor']*(lsize-1)+['Normal']*1
                else:
                    sname=['Tumor']
                    
                
                for j,ds in enumerate(alls):
                    dt=ds.split(',')
                    y2=[]
                    for y in dt:
                        if '.' in y:
                            y=float(y)
                            if y!=0:
                                y2.append(y)
                    if len(y2)>0:
                        create_figure_file(y2,pdf,title='Event-'+line[0]+'-Sample_set-'+sname[j]+'-log2(gene/EGFR)')

            if ';' in line[38] and ',' in line[38]:
                dt=line[38].split(';')[1]
                dt=dt.split(',')
                y3=[]
                for y in dt:
                    if '.' in y:
                        y=float(y)
                        if y!=0:
                            y3.append(y)
                if len(y3)>0:
                    create_figure_file(y3,pdf,title='Event-'+line[0]+'-Insilico_qPCR')

        except Exception as e:
            print(e)


            
        for word in line:
            filtered_file.write("%s\t" %(word))
        filtered_file.write("\n")
    for item in event_dic:
        eventf.write("%s\n" %(item))

    for item in gene_dic:
        genef.write("%s\n" %(item))



CSSHTMLBuildFilter(filteredfile,taskid)

