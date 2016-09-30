# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:05:18 2016

@author: Anindya Bhattacharya, Department of Computer Science and Engineering, UCSD, 9500 Gilman Drive, La Jolla, CA 92093-0404
Email:anb013@eng.ucsd.edu

Parameters: python pep_reporter.py Input_File_Tab_Separated_from_Events  Output_File execution_tag_from_ENOSI_for_event_spectra_link_creation

"""

import sys 

import os

import pep_scan as psc
import protein_annot as pa

from build_css_html_table_v5_new import CSSHTMLBuildFilter



dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'


#if __name__ == "__main__":


try:
    names=pa.getfilename_w_annot(sys.argv[1:])
    if names.has_key("inputfile"):
        inputfilename=names["inputfile"]
    else:
        print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'

        sys.exit(2)
        
    if names.has_key("outputfile"):
        outputfilename=names["outputfile"]
    else:
        print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'
        sys.exit(2)
    if names.has_key("taskid"):
        taskid=names["taskid"]
    else:
        taskid=''
    if names.has_key("afile"):
        annotf=names["afile"]
    else:
        annotf=''
        
    if names.has_key("expfile"):
        expfile=names["expfile"]
    else:
        expfile=''

    if names.has_key("qpcrfile"):
        qpcrfile=names["qpcrfile"]
    else:
        qpcrfile=''
        
        
except Exception as e:  
	print 'python pep_reporter.py -i <Input file> -t <Task ID from ENOSI> -o <Output file> -a <Optional Event Description File> -e <Optional Expression Data> -p <Optional qPCR file>'
	sys.exit(5)



inputfile=open(inputfilename,"r")

outputfile=open(outputfilename,"w")

htmlinputfilename=os.path.realpath(outputfilename)

psc.petide_file_scan(inputfile,outputfile,outputfilename,dir_path,annotf,expfile,qpcrfile)

CSSHTMLBuildFilter(htmlinputfilename,taskid)


