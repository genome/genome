package Genome::Model::Build::View::Annotation::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::View::Annotation::Html {
    is => 'Genome::View::Status::Html',
    has_constant => [
        perspective => {
            value => 'annotation',
        },
    ],
};

sub _generate_content {
    my $self = shift;
    my $build = $self->subject;

    my $return_value = <<'HTML'

<!doctype html>
<html>
  <head>
    <meta charset="utf-8">
    <title>Annotation Search</title>
      <script src="/res/js/pkg/jquery.js"></script>
      <script src="/res/js/pkg/jquery-ui.js"></script>

      <script src="/res/js/pkg/dataTables/media/js/jquery.dataTables.min.js"></script>
    <link href="/res/css/jquery-ui-overrides.css" type="text/css" rel="stylesheet" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/dataTables.css" type="text/css" media="screen, projection" />
    <link rel="shortcut icon" href="/res/img/gc_favicon.png" type="image/png" />
    <link rel="stylesheet" href="/res/css/blueprint/screen.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/blueprint/print.css" type="text/css" media="print" />
    <link rel="stylesheet" href="/res/css/master.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/buttons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/icons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/forms.css" type="text/css" media="screen, projection" />
    <link type="text/css" href="/res/js/pkg/jquery-ui-1.8.1.custom/css/gsc-theme/jquery-ui-1.8.1.custom.css" rel="stylesheet" />
    <link rel="stylesheet" href="/res/css/genome_model_build_annotation.css" type="text/css" media="screen, projection" />
<!--    <style type="text/css" media="screen">
		
		/*
		 * Override styles needed due to the mix of three different CSS sources! For proper examples
		 * please see the themes example in the 'Examples' section of this site
		 */
/*      .dataTables_info { padding-top: 0; }
        .dataTables_paginate { padding-top: 0; }
        .css_right { float: right; }
        #example_wrapper .fg-toolbar { font-size: 0.8em }
        #theme_links span { float: left; padding: 2px 10px; }
*/	</style> -->
  </head>
  <body>
  <div class="page">
<!--template: /html/common_includes/components.xsl name:control_bar_view-->
      <div xmlns:rest="urn:rest" xmlns:str="http://xsltsl.org/string" class="control_bar view shadow">
        <div class="control_bar_menu" id="bar_menu">
<!--template: /html/common_includes/components.xsl name:control_bar_cache_info-->
          <div class="cache_info">
            <div class="cache_time">
              <p>
          View generated<br /><strong><span id="updatedOn"></span></strong></p>
            </div>
            <div class="cache_refresh">
              <a class="btn_cache_refresh rounded" id="refreshCache" title="Refresh Cache"></a>
            </div>
          </div>
<!--template: /html/common_includes/components.xsl name:control_bar_menu-->
          <ul class="app_menu">
            <li>
              <a href="/view/genome/status.html" class="app btn shadow"><div class="icon"><img src="/res/img/icons/app_deprecated_search_16.png" width="16" height="16" /></div>
          Deprecated Search
        </a>
            </li>
            <li>
              <a href="/view/genome/search/status.html" class="app btn shadow"><div class="icon"><img src="/res/img/icons/app_analysis_search_16.png" width="16" height="16" /></div>
          Analysis Search
        </a>
            </li>
          </ul>
        </div>
        <div class="control_bar_base" id="bar_base"> </div>
      </div>
<!--template: /html/common_includes/components.xsl:view_header-->
      <div xmlns:rest="urn:rest" xmlns:str="http://xsltsl.org/string" class="header rounded-bottom gradient-grey shadow">
        <div class="container">
          <div class="title span-24 last">
            <h1 class="no_icon" style="margin-left: 0;">Query Results</h1>
          </div>
        </div>
      </div>
      <div class="content rounded shadow">
        <div class="container">
    <div class="search">
        <form action="javascript:return false;" id="searchForm">
            <input type="text" id="geneName"/>
            <input type="submit" value="Search" id="geneSearch"/><br />
			<input type="radio" name="searchType" value="gene" checked/> Gene<br />
			<input type="radio" name="searchType" value="grep" /> Full-text grep
        </form>
    </div>
    <div class"progress">
        <div id="progressIndex">
            Idle
        </div>
        <div id="progressSearch">
            Idle
        </div>
    </div>
    <div id="annotationTableHolder" style="width: 1250px;">
        <table cellpadding="0" cellspacing="0" border="0" class="display" id="annotationTable">
            <thead>
                <tr>
                    <th rowspan="2">&nbsp;</th>
                    <th rowspan="2">Location</th>
                    <th rowspan="2">Chr</th>
                    <th rowspan="2">Start</th>
                    <th rowspan="2">Stop</th>
                    <th rowspan="2">Reference</th>
                    <th rowspan="2">Variant</th>
                    <th rowspan="2">Variant Type</th>
                    <th rowspan="2">Gene</th>
                    <th rowspan="2">Strand</th>
                    <th rowspan="2">c pos</th>
                    <th rowspan="2">amino acid change</th>
                    <th rowspan="2">ucsc cons</th>
                    <th rowspan="2">domain</th>
                    <th rowspan="2">all domains</th>
                    <th rowspan="2">deletion substructures</th>
                    <th colspan="7" rowspan="1" class="ui-state-default">Transcript</th>
                </tr>
                <tr>
                    <th>Name</th>
                    <th>Species</th>
                    <th>Source</th>
                    <th>Version</th>
                    <th>Status</th>
                    <th>Error</th>
                    <th>trv Type</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td colspan="21" class="dataTables_empty">Loading data from server</td>
                </tr>
            </tbody>
            <!-- <tfoot>
                <tr>
                    <th>MINUS</th>
                    <th>1</th>
                    <th>2</th>
                    <th>3</th>
                    <th>4</th>
                    <th>5</th>
                    <th>6</th>
                    <th>7</th>
                    <th>8</th>
                    <th>9</th>
                    <th>10</th>
                    <th>11</th>
                    <th>12</th>
                    <th>13</th>
                    <th>14</th>
                    <th>15</th>
                    <th>16</th>
                    <th>17</th>
                    <th>18</th>
                    <th>19</th>
                    <th>20</th>
                    <th>21</th>
                </tr>
            </tfoot> -->
        </table>
    </div>
    <script>
    
        function fnFormatDetails(oTable, nTr) {
            var aData = oTable.fnGetData(nTr);
            var sOut = '';
            
            sOut += '<div class="details_chr_title">chr</div>';
            sOut += '<div class="details_chr_val">'+aData[2]+'</div>';
            
            sOut += '<div class="details_top">';
            
            sOut += '<div class="details_key_info_title">';
            sOut += 'Start:<br/>Stop:';
            sOut += '</div>';
            
            sOut += '<div class="details_key_info_val">';
            sOut += aData[3] + '<br/>' + aData[4];
            sOut += '</div>';
            
            sOut += '<div class="details_key_info_title">';
            sOut += 'Reference:<br/>Variant:';
            sOut += "</div>";
            
            sOut += '<div class="details_key_info_val">';
            sOut += aData[5] + '<br/>' + aData[6];
            sOut += '</div>';
            
            sOut += '<div class="details_subchr">';
            if (aData[7] == 'INS') {
                sOut += 'Insertion';
            } else if (aData[7] == 'DEL') {
                sOut += 'Deletion';
            } else if (aData[7] == 'SNV') {
                sOut += 'SNV';
            } else {
                sOut += aData[7];
            }
            sOut += ' in <span class="details_gene">' + aData[8] + '</span>';
            sOut += '</div>';
            sOut += '</div>';
            
            
            sOut += '<div class="details_transcript">';
            sOut += '<div class="details_transcript_title">';
            sOut += 'Transcript';
            sOut += '</div>';

            sOut += '<div class="details_transcript_subtitle">';
            sOut += 'Name:<br/>Species:<br/>Source:<br/>Version:<br/>Status:<br/>Error:';
            sOut += '</div>';

            sOut += '<div class="details_transcript_val">';
            sOut += aData[16] + '<br/>' + aData[17] + '<br/>' + aData[18] + '<br/>' + aData[19] + '<br/>' + aData[20]+ '<br/>' + aData[21];
            sOut += '</div>';

            sOut += '</div>';

            sOut += '<div class="details_other_title">';
            sOut += 'Other';
            sOut += '</div>';

            sOut += '<div class="details_other_subtitle">';
            sOut += 'trv Type:<br/>strand:<br/>c Pos:<br/>Amino Acid Change:<br/>UCSC cons:<br/>Deletion Substructures:';
            sOut += '</div>';

            sOut += '<div class="details_other_val">';
            sOut += aData[22] + '<br/>' + aData[9] + '<br/>' + aData[10] + '<br/>' + aData[11] + '<br/>' + aData[12] + '<br/>' + aData[15];
            sOut += '</div>';

            sOut += '<div class="details_domain_title">';
            sOut += 'Domains:'
            sOut += '</div>';

            sOut += '<div class="details_domain_val">';
            
            var domains = new Array();
            domains = aData[13].split(",");
            if (domains[0] == '-' || domains[0] == 'NULL') {
                sOut += domains[0];
            } else {
                sOut += '<ul class="details_list">';
                for (var domain in domains) {
                    sOut += '<li>' + domains[domain] + '</li>';
                }
                sOut += '</ul>';
            }
            sOut += '</div>';

            sOut += '<div class="details_domain_title">';
            sOut += 'All Domains:';
            sOut += '</div>';

            sOut += '<div class="details_domain_val">';

            var alldomains = new Array();
            alldomains = aData[14].split(",");
            if (alldomains[0] == '-' || alldomains[0] == 'NULL') {
                sOut += alldomains[0];
            } else {
                sOut += '<ul class="details_list">';
                for (var alldomain in alldomains) {
                    sOut += '<li>' + alldomains[alldomain] + '</li>';
                }
                sOut += '</ul>';
            }
            sOut += '</div>';

            return sOut;

        }
        
        function processAnnotationLine(row) {            
            var expand = '<img src="/res/img/build_annotation/details_open.png">';
            var loc = row[0] + ":" + row[1] + "-" + row[2]; // chr:start-stop
            var newrow = new Array();
            newrow.push(expand);
            newrow.push(loc);
            newrow.push(row[0]); // chr
            newrow.push(row[1]); // start
            newrow.push(row[2]); // stop
            newrow.push(row[3]); // ref
            newrow.push(row[4]); // variant
            newrow.push(row[5]); // variant type
            newrow.push(row[6]); // gene
            newrow.push(row[11]); // strand
            newrow.push(row[14]); // c_pos
            newrow.push(row[15]); // amino acid change
            newrow.push(row[16]); // ucsc cons
            newrow.push(row[17]); // domain
            newrow.push(row[18]); // all domains
            newrow.push(row[19]); // deletion substructures
            newrow.push(row[7]); // transcript name
            newrow.push(row[8]); // transcript species
            newrow.push(row[9]); // transcript source
            newrow.push(row[10]); // transcript version
            newrow.push(row[12]); // transcript status
            newrow.push(row[20]); // transcript error
            newrow.push(row[13]); // trv type
            
            // row.unshift(pos);
            // var bp = row[3] + " / " + row[4];
            // var type = row[5];
            // var gene = row[6];
            // var loc = row[7];
            // var species = row[8];
            // var genedb = row[9];
            // var dbloc = row[10];
            
            // return [expand, pos, bp, type, gene, loc, species, genedb, dbloc];
            return newrow;
        }
        
        $(document).ready(function(){
            $.ajaxSetup ({ cache: false });
            
            var build_id =
HTML
 . $build->id . <<'HTML'
;            var base_url = "/view/genome/model/build/annotation.json"
            
            var build_type = 1; 
            
            var ajax_load = "<em>loading</em>";
            var ajax_grep = "<em>grepping</em>";
            var ajax_grab = "<em>grabbing</em>";
            var ajax_done = "<strong>complete</strong>";
            var ajax_error = "<strong>error</strong>";
            
            $("#progressIndex").html(ajax_load);
            
            $.getJSON(
                base_url,
                {
                    "build_id": build_id,
                    "-request_index": build_type
                },  
                function(json) {
                    index_status = json
                    if ((build_type == 2 && index_status['built'] == 1) || (index_status['exists'])) {
                        $("#progressIndex").html(ajax_done);                        
                    } else {
                        $("#progressIndex").html(ajax_error)
                    }
                }  
            );
                        
            // Initialize DataTables
            theTable = $("#annotationTable").dataTable({
                // "bFilter": false,
                "bProcessing": true,
                "bServerSide": true,
                "sAjaxSource": base_url,
                "sPaginationType": "full_numbers",
                "bJQueryUI": true,
                "aaSorting": [[1, "asc" ]],
                "aoColumnDefs": [
                    { "bVisible": false, "aTargets": [ 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 21, 22 ] },
                    { "bSearchable": false, "aTargets" : [ 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]},
                    { "bSortable": false, "aTargets": [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]},
                    { "sWidth": "200px", "aTargets": [1] }
                ],
                "fnServerData": function ( sSource, aoData, fnCallback ) {
                    // flatten the list, and get the filter string if there is one
                    var aoDataFlat = new Array();
                    // var filterString = "";
                    for (var pair in aoData) {
                        // if (aoData[pair]["name"] == "sSearch") {
                        //     filterString = aoData[pair]["value"];
                        // }
                        aoDataFlat.push(aoData[pair]["name"]);
                        aoDataFlat.push(aoData[pair]["value"]);
                    }
					aoDataFlat.push("sSearchType");
					aoDataFlat.push($('input:radio[name=searchType]:checked').val())
					
	// alert($('input:radio[name=searchType]:checked').val())
					
                    var aoDataString = aoDataFlat.join(",");
                    var aoDataNew = new Array();
                    aoDataNew.push( {
                        "name": "-datatables_params",
                        "value": aoDataString 
                    }, {
                        "name": "build_id",
                        "value": build_id
                    } );
                   // alert(""+sSource+" "+aoDataNew); 
                    $.getJSON( sSource, aoDataNew, function (json) {
                        if ($("#progressSearch").html() == ajax_grab || $("#progressSearch").html() == ajax_grep) {
                            $("#progressSearch").html(ajax_done);
                        }
                        
                        var processedData = new Array();
                        for (var curRow in json['aaData']) {
                            processedData.push(processAnnotationLine(json['aaData'][curRow]));
                        }
                        json['aaData'] = processedData;
                                                
                        /* Inform DataTables */
                        fnCallback(json);
                    } );
                }
            });
            
            // Bind search form to DataTables
            $("#searchForm").submit(function() {
                query = $("#geneName").val().toUpperCase();
                theTable.fnFilter(query)
            });
            
            // Add event listener for opening and closing details
            $('#annotationTable tbody td img').live('click', function () {
        		var nTr = this.parentNode.parentNode;
        		if ( this.src.match('details_close') ) {
        			/* This row is already open - close it */
        			this.src = this.src.replace("details_close","details_open");
        			theTable.fnClose( nTr );
        		} else {
        			/* Open this row */
        			this.src = this.src.replace("details_open","details_close");
        			theTable.fnOpen( nTr, fnFormatDetails(theTable, nTr), 'details' );
        		}
        	} );
        });
    </script>
    </div>
    </div>
  </body>
</html>

HTML
;
    return $return_value;
}
