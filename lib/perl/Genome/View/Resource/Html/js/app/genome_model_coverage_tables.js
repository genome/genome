// INCLUDED FROM:
// xsl/html/coverage/genome_model.xsl

$(document).ready(function() {
    $('#alignment-lister').dataTable( {
                    "oLanguage": { "sUrl": "/res/tgi.txt" },
                     "sDom": 'T<"clear">lfrtip',
                         "oTableTools": {
                             "sSwfPath": "/res/js/pkg/TableTools/media/swf/copy_cvs_xls_pdf.swf"
                         },
                         "sPaginationType": "full_numbers",
                         "aoColumns": [  null,
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                      ]
                         });
    $('#coverage-depth-lister').dataTable( {
                    "oLanguage": { "sUrl": "/res/tgi.txt" },
                     "sDom": 'T<"clear">lfrtip',
                         "oTableTools": {
                             "sSwfPath": "/res/js/pkg/TableTools/media/swf/copy_cvs_xls_pdf.swf"
                         },
                         "sPaginationType": "full_numbers",
                         "aoColumns": [  null,
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                      ]
                     });
    $('#coverage-summary-lister').dataTable( {
                    "oLanguage": { "sUrl": "/res/tgi.txt" },
                     "sDom": 'T<"clear">lfrtip',
                         "oTableTools": {
                             "sSwfPath": "/res/js/pkg/TableTools/media/swf/copy_cvs_xls_pdf.swf"
                         },
                         "sPaginationType": "full_numbers",
                         "aoColumns": [  null,
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                         {"sType": "percent"},
                                      ]
                     });
    $('#enrichment-factor-lister').dataTable( {
                    "oLanguage": { "sUrl": "/res/tgi.txt" },
                     "sDom": 'T<"clear">lfrtip',
                         "oTableTools": {
                         "sSwfPath": "/res/js/pkg/TableTools/media/swf/copy_cvs_xls_pdf.swf"
                     },
                         "sPaginationType": "full_numbers",
                         "aoColumns": [  null,
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                         {"sType": "formatted-num"},
                                      ]
                     });
});
