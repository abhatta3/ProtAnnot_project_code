Format:

python ./src/pep_reporter.py -i [inputfile from ENOSI] -o [outputfilename] -t [task id for events from enosi on proteosafe output page]

For colorectal:

nohup python ./src/pep_reporter.py -i ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/Colon_Proteogenomic_Events.txt -o Colorectal -t 0d298cacb6d4486b9a94ab1762e707d0 -a ./src/Datafiles/Additional_data/Peptide_list_sampleID_Colon.txt  -e ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/colon_protein_expression.txt -p ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/QPCR_colorectal.txt &



For ovarian:

nohup python ./src/pep_reporter.py -i ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/Ovarian_Proteogenomic_Events.txt -o Ovarian -t f0622b5022da4d43ad1e069f40429c35  -a ./src/Datafiles/Additional_data/Peptide_list_sampleID_Ovarian.txt  -e ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/ovarian_protein_expression.txt -p ./Example_input_files/Enosi_main_result_file_example_Colon_Ovarian/QPCR_ovarian.txt  &




For colorectal selection:

python ./src/Selection_of_membrane_bound.py  Colorectal  0d298cacb6d4486b9a94ab1762e707d0 


For ovarian selection:

python ./src/Selection_of_membrane_bound.py Ovarian f0622b5022da4d43ad1e069f40429c35





