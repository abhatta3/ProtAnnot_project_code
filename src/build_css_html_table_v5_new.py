#CPTAC_Colon_Proteogenomic_Events_membrane = '12_14_2015_CPTAC_Colon_Proteogenomic_Events_membrane.txt'
#CSSHTMLBuild(CPTAC_Colon_Proteogenomic_Events_membrane)


def CSSHTMLBuildFilter(inputfile,spect_task_id=''):
    css_table_header = '''<!DOCTYPE html>\n\
    <html>\n\
    <head>\n\

    <link href="http://cdn.bootcss.com/bootstrap/3.3.0/css/bootstrap.min.css" type="text/css" rel="stylesheet">\n\
    <script type="text/javascript" src="https://dl.dropboxusercontent.com/u/63775276/sorttable.js"></script>\n\
    <style type="text/css">\n\
    table.sortable {\n\
        font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;\n\
        width: 100%;\n\
        border-collapse: collapse;\n\
    }\n\
    \n\
    table.sortable td, #mutated_peptides th {\n\
        font-size: 1em;\n\
        border: 1px solid #ddd;\n\
        padding: 3px 7px 2px 7px;\n\
    }\n\
    \n\
    table.sortable tr:nth-child(odd){background-color: #f2f2f2}\n\
    \n\
    table.sortable th {\n\
        font-size: 1em;\n\
        text-align: center;\n\
        padding-top: 5px;\n\
        padding-bottom: 4px;\n\
        background-color: #4CAF50;\n\
        color: #ffffff;\n\
        overflow: hidden;\n\
        white-space: nowrap;\n\
    }\n\
    \n\
    table.sortable td.center {\n\
        text-align: center;\n\
    }\n\
    \n\
    table.sortable td.canceratlas {\n\
        color: #8cd68f;\n\
        background-color: #8cd68f;\n\
        text-align: center;\n\
    }\n\annotfile
    \n\
    </style>\n\


    </head>\n\
    <body>\n\
    \n
    <form class="form-horizontal" role="form" style="width: 80%;margin: auto;margin-top:50px;">\n\

    <div class="form-group">\n\
        <div class="col-sm-12">\n\
            <h2 class="text-center">Select Items</h2>\n\
        </div>\n\
    </div>\n\
    <div class="form-group" style="border: 1px solid #ddd">\n\
        <table id="table2"  class="table col-sm-12">\n\
        </table>\n\
    </div>\n\
    </form>\n\
    <script src="http://code.jquery.com/jquery-1.11.1.min.js"></script>\n\
    <script type="text/javascript">\n\
    
    (function ($) {
    $.fn.multiSelect = function (options) {
        $.fn.multiSelect.init($(this), options);
    };

    $.extend($.fn.multiSelect, {
        defaults: {
            actcls: 'active',
            selector: 'tbody tr', 
            except: ['tbody'], 
            statics: ['.static'], 
            callback: false 
        },
        first: null, 
        last: null, 
        init: function (scope, options) {
            this.scope = scope;
            this.options = $.extend({}, this.defaults, options);
            this.initEvent();
        },
        checkStatics: function (dom) {
            for (var i in this.options.statics) {
                if (dom.is(this.options.statics[i])) {
                    return true;
                }
            }
        },
        initEvent: function () {
            var self = this,
                scope = self.scope,
                options = self.options,
                callback = options.callback,
                actcls = options.actcls;

            scope.on('click.mSelect', options.selector, function (e) {
                if (!e.shiftKey && self.checkStatics($(this))) {
                    return;
                }

                if ($(this).hasClass(actcls)) {
                    $(this).removeClass(actcls);
                } else {
                    $(this).addClass(actcls);
                }

                if (e.shiftKey && self.last) {
                    if (!self.first) {
                        self.first = self.last;
                    }
                    var start = self.first.index();
                    var end = $(this).index();
                    if (start > end) {
                        var temp = start;
                        start = end;
                        end = temp;
                    }
                    $(options.selector, scope).removeClass(actcls).slice(start, end + 1).each(function () {
                        if (!self.checkStatics($(this))) {
                            $(this).addClass(actcls);
                        }
                    });
                    window.getSelection ? window.getSelection().removeAllRanges() : document.selection.empty();
                } else if (!e.ctrlKey && !e.metaKey) {
                    $(this).siblings().removeClass(actcls);
                }
                self.last = $(this);
                $.isFunction(callback) && callback($(options.selector + '.' + actcls, scope));
            });
            $(document).on('click.mSelect', function (e) {
                for (var i in options.except) {
                    var except = options.except[i];
                    if ($(e.target).is(except) || $(e.target).parents(except).size()) {
                        return;
                    }
                }
                scope.find(options.selector).each(function () {
                    if (!self.checkStatics($(this))) {
                        $(this).removeClass(actcls);
                    }
                });
                $.isFunction(callback) && callback($(options.selector + '.' + actcls, scope));
            });

            $(document).on('keydown.mSelect', function (e) {
                if ((e.keyCode == 65) && (e.metaKey || e.ctrlKey)) {
                    $(options.selector, scope).each(function () {
                        if (!self.checkStatics($(this))) {
                            $(this).addClass(actcls);
                        }
                    });
                    $.isFunction(callback) && callback($(options.selector + '.' + actcls, scope));
                    e.preventDefault();
                    return false;
                }
            });
            $(document).on('keyup.mSelect', function (e) {
                if (e.keyCode == 16) {
                    self.first = null;
                }
            });
        }
    });
 })(jQuery);
    </script>\n\
    <script type="text/javascript">\n\
    $(function () {\n\
        $('#table1').multiSelect({\n\
            actcls: 'info', \n\
            selector: 'tbody tr', \n\
            except: ['tbody'], \n\
            statics: ['.danger', '[data-no="1"]'], \n\
            callback: function (items) {\n\
                $('#table2').empty().append(items.clone().removeClass('info').addClass('success'));\n\
            }\n\
        });\n\
    })\n\
    </script>\n\
    <script type="text/javascript">\n\

    var _gaq = _gaq || [];\n\
    _gaq.push(['_setAccount', 'UA-36251023-1']);\n\
    _gaq.push(['_setDomainName', 'jqueryscript.net']);\n\
    _gaq.push(['_trackPageview']);\n\

    (function() {\n\
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;\n\
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';\n\
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);\n\
    })();\n\
    </script>\n\

       <div class="form-group">\n\
            <div class="col-sm-12">\n\
                <span class="help-block">Sort table with Left Click at Header. @ Click event to see spectra. You can select with <b>Ctrl</b>, <b>Shift</b>, <b>Common</b> Or <b>Ctrl + A</b>.</span>\n\
            </div>\n\
        </div>\n\


        <div class="form-group">\n\
            <div class="col-sm-12">\n\
                <h1 class="text-center">Proteogenomics Annotation Data Table Towards Target Selection</h1>\n\
            </div>\n\
        </div>\n\
        <div class="form-group" style="border: 1px solid #ddd">\n'''

    table_column_name ='''<table id="table1"  class="sortable">\n\
    	
    <thead>\n\
      <tr>\n\
        <th>Event ID<br>[@]</th>\n\
        <th>Event Type</th>\n\
        <th>Peptide</th>\n\
        <th>Ref_aa</th>\n\
        <th>Mut_aa</th>\n\
        <th>Ref_Location</th>\n\

        <th>Mut_loc_at_peptide</th>\n\
        <th>Total Number of<br>Overlapping Gene</th>\n\

        <th>Overlapping<br>Gene</th>\n\
        <th>Overlapping Protein<br>(UniProtKB_AC)</th>\n\


        <th>UniProtKB <br> Key Words</th>\n\
        <th>Protein_position</th>\n\

        <th>Protein Atlas<br>(Cancer Atlas)</th>\n\

        <th>Chromosome</th>\n\

        <th>Genomic_Location</th>\n\

        <th>#Novel peptides<br>in overlapping gene</th>\n\
        <th>#Known peptides<br>in overlapping gene</th>\n\

        <th>Total<br>Spec Count</th>\n\
        <th>Number of possible<br>genomic locations</th>\n\

        <th>WT_Charge</th>\n\
        <th>MT_Charge</th>\n\
        <th>MT-WT_Charge</th>\n\
        <th>WT_Surface</th>\n\
        <th>MT_Surface</th>\n\
        <th>MT-WT_Surface</th>\n\
        <th>Kyte-Doolittle<br>Hydrophilic(>0)</th>\n\
        <th>Hopp-Woods<br>Hydrophilic(>0)</th>\n\
        <th>dbSNP_Overlap</th>\n\
        <th>#Spectra_Samples</th>\n\
        <th>#RNA_Samples</th>\n\
        <th>Max_Read_Depth</th>\n\
        <th>Mutation_Event_Position</th>\n\
        <th>log2(Gene/EGFR)<br>Protein Expression</th>\n\
        <th>Gene Rank, EGFR Rank <br> Protein Expression</th>\n\
        <th>In-Silico PCR</th>\n\ 

      </tr>\n\
     </thead>\n\
     <tbody>\n'''
     
     

    #<th>Reference_Peptide</th>\n\
    htmlfile=inputfile+'.html' 
        
    with open(inputfile,'r') as data, open(htmlfile,'w') as html_file:
        
        data_lines = data.readlines()[1:]
        html_file.write(css_table_header)

        html_file.write(table_column_name)

            
        for item in data_lines:
            
            item_list = item.split('\t')
            
            html_file.write('  <tr>\n')

            for i, item in enumerate(item_list):
                #skip reference peptide. It is in the text file
                if i==3:
                    continue
                if i>38:
                    break;

                #Event_ID
                if i!=15 and i!=13 and i!=36:
                    dt = item.strip('\n')
                    if dt=='':
                        dt='-'
                        
                    if i==0:
                        if spect_task_id!='':
                            dt='<a href="http://proteomics2.ucsd.edu/ProteoSAFe/result.jsp?task='+spect_task_id+'&view=group_by_event#%7B%22%23Num_lowerinput%22%3A%22'+dt+'%22%2C%22%23Num_upperinput%22%3A%22'+dt+'%22%7D" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==9:
                        dt='<a href="'+item_list[15]+'" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==10:
                        dt='<a href="http://www.uniprot.org/uniprot/'+dt+'" target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==12: 
                        dt='<a href="'+item_list[13]+'"  target="_blank"><FONT COLOR="#8cd68f">'+dt+'</FONT></a>'
                    elif i==14:
                        dt='<a href="'+dt+'" target="_blank"><FONT COLOR="#8cd68f">Atlas</FONT></a>'
                    if i==28:
                        dt=dt.split(';')[0]
                        hval=float(dt)
                        if hval>0:
                            hlabel="Hydrophobic;"
                            dt=hlabel+dt
                            
                    if i==35 or i==37:
                        combinedt=''
                        for part in dt.split(';'):
                            combinedt+=part+'<br>'
                        dt=combinedt
                            

                    if i==29:
                        dt=dt.split(';')[0]
                        hval=float(dt)
                        if hval>0:
                            hlabel="Hydrophilic;"
                            dt=hlabel+dt

                    if i==38:
                        dt=dt.split(';')[0]
                            
                    if i==30:
                        dt=dt.split('#')
                        dbl=''
                        for d in dt:
                            if dbl=='':
                                
                                if 'dbSNP:rs' in d:
                                    snpid=d.split(':')[1].translate(None, '.)rs')
                                    d='<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='+snpid+'">'+d+'</a>'
                                
                                dbl=d
                            else:
                                
                                if 'dbSNP:rs' in d:
                                    snpid=d.split(':')[1].translate(None, '.)rs')
                                    d='<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='+snpid+'">'+d+'</a>'
                                dbl=dbl+'<br><br>'+d
                        dt=dbl
                    
                    html_file.write('    <td class="center">%s</td>\n'%(dt))
                        
            html_file.write('  </tr>\n')
        html_file.write('</tbody>\n</table>\n\n</div>\n</body>\n</html>\n')









