$(document).ready(function() {
    oTable = $('#individuals').dataTable({
        "oLanguage": { "sUrl": "/res/tgi.txt" },
        "sPaginationType": "full_numbers",
        "bSort": false,
        "bProcessing": false,
        "sAjaxSource": "/viewajax/genome/individual/detail.json?id=" + window.individual_id,
        "fnDrawCallback": function ( oSettings ) {
            if ( oSettings.aiDisplay.length == 0 )
            {
                return;
            }

            // sort individual by nomenclature (this is probably not the best method,
            // but couldn't figure out how to sort the json before dataTables gets
            // it into its oSettiongs.aoData array.
            oSettings.aoData.sort(function(a,b) {
                if (a._aData[0].toLowerCase() == b._aData[0].toLowerCase()) {
                    return 0;
                }
                return a._aData[0].toLowerCase() > b._aData[0].toLowerCase() ? 1 : -1;
            });

            var nTrs = $('#individuals tbody tr');
            nTrs.sort(function(a,b) {
                return a.innerText == b.innerText;
            });
    
            var iColspan = nTrs[0].getElementsByTagName('td').length;
            var sLastGroup = "";
            for ( var i=0 ; i<nTrs.length ; i++ )
            {
                var iDisplayIndex = oSettings._iDisplayStart + i;
                var sGroup = oSettings.aoData[ oSettings.aiDisplay[iDisplayIndex] ]._aData[0];
                var nomenclatureId= oSettings.aoData[ oSettings.aiDisplay[iDisplayIndex] ]._aData[1];

                if ( sGroup != sLastGroup )
                {
                    var nGroup = document.createElement( 'tr' );
                    var nCell = document.createElement( 'td' );
                    nCell.colSpan = iColspan;
                    nCell.className = "group";
                    if (nomenclatureId == 'null') {
                        nCell.innerHTML = '[no defined nomenclature]'
                    } else {
                        nCell.innerHTML = '<a class="link" href="/view/genome/nomenclature/set/create.html#id=' 
                                            + nomenclatureId
                                            + '">' 
                                            + sGroup 
                                            + "</a>";
                    }
                    nGroup.appendChild( nCell );
                    nTrs[i].parentNode.insertBefore( nGroup, nTrs[i] );
                    sLastGroup = sGroup;
                }
            }

        },
        "aoColumnDefs": [
            { "bVisible": false, "aTargets": [ 0 ] },
            { "bVisible": false, "aTargets": [ 1 ] }
        ],
        "aaSortingFixed": [[ 0, 'asc' ]],
        "aaSorting": [[ 1, 'asc' ]],
        "sDom": 'lfr<"giveHeight"t>ip'
    });

} );
