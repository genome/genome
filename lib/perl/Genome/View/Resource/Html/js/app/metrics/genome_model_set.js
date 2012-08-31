

   
$(document).ready(function() {

    //TableToolsInit.sSwfPath = "/res/js/pkg/dataTables-1.7/extras/TableTools/media/swf/ZeroClipboard.swf";

    var ajaxData;
    var source = "/viewajax/genome/model/set/metrics.json?";

    // lean on getTarget function defined in view html to find out the set search url
    source += getTarget();

    $(document).data('updatedOn', new Date(1320090254000));

    $.ajax({
        "dataType": 'json',
        "type": "GET",
        "url": source,
        "error" : function (jqXHR, textStatus, errorThrown) {
            alert("failed to get data for metrics table");
        },
        "success": function (data, textStatus, jqXHR) {
            var columns = data.aoColumns;
            var header_row = $('#table-header');
            for (i=0;i<columns.length;i++) {
                header_row.append($('<th>' + columns[i].mDataProp + '</th>'));
            }
            data.bJQueryUI = true;
            data.oLanguage = { "sUrl": "/res/tgi.txt" };
            data.sScrollX = '100%';
            data.sDom = 'T<"clear">lfrtip';
            data.oTableTools = {
                             "sSwfPath": "/res/js/pkg/TableTools/media/swf/copy_cvs_xls_pdf.swf",
                             "aButtons": [
//                                    { "sExtends" : "print",
//                                      "fnClick": function() { alert ('click');
//                                      }
//                                    },
                                    { "sExtends" : "csv",
                                      },
                                    { "sExtends" : "pdf",
                                      },
                            ]
                         };
            $('#loading-message').hide();
            $('#model_metrics').dataTable(data);
        }
    });

});




