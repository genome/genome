$(document).ready(function() {


    var url = "/viewajax/genome/sample/set/detail.json?project_id in=" + window.project_id

alert("url = " + url);
return;

    oTable = $('#samples').dataTable({
        "oLanguage": { "sUrl": "/res/tgi.txt" },
        "sPaginationType": "full_numbers",
        "bSort": false,
        "bProcessing": false,
        "sAjaxSource": url,
        "aoColumnDefs": [
            { "bVisible": false, "aTargets": [ 0 ] }
        ],
        "aaSortingFixed": [[ 0, 'asc' ]],
        "aaSorting": [[ 1, 'asc' ]],
        "sDom": 'lfr<"giveHeight"t>ip'
    });

} );
