

$(document).ready(function() {


    $(document).ready(function(){
        $('#projects').dataTable({
            "sScrollX": "100%",
            "sScrollInner": "110%",
            "bJQueryUI": true,
            "sPaginationType": "full_numbers",
            "bStateSave": true,
            "iDisplayLength": 25,
            "sDom": '<"H"l<"toolbar">f>t<"F"ip>',
            "oLanguage": { "sUrl": "/res/tgi.txt" },
        });

        $("div.toolbar").html('');
     });


/*
    var url = "/view/genome/sys/user/detail.json?id=" + window.email;
    $.ajax({
        "dataType": 'json',
        "type": "GET",
        "url": url,
        "error" : function (jqXHR, textStatus, errorThrown) {
            alert("failed to get data for projects table");
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
                                    { "sExtends" : "csv",
                                      },
                                    { "sExtends" : "pdf",
                                      },
                            ]
                         };
            $('#projects').dataTable(data);
        }
    });
*/

});



