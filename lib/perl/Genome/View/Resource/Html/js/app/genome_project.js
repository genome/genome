

$(document).ready( function() {

    var parts_table = $('#parts_table');
    var toggleCheckboxes = function(v) { $(".partCheckbox").prop("checked",v); };

    var setButtonStates = function() {
        var n = $(".partCheckbox:checked").length;
        if (n > 0) {
            $(".actionButton").removeAttr("disabled");
        } else {
            $(".actionButton").attr("disabled","false");
        }
    };

    var removeSelectedProjectParts = function() {

        var ids = new Array();
        $(".partCheckbox:checked").each(function(i, p) {
            ids.push(p.id);
        });
        $(".partCheckbox:checked").attr("disabled",true);

        var removeUrl = "/view/genome/project-part";
        $.ajax({
            "dataType": "json",
            "type": "DELETE",
            "url": removeUrl,
            "data": { "ids": ids },
            "error" : function (jqXHR, textStatus, errorThrown) {
                alert("request to delete parts failed"); 
            },
            "success" : function(data, textStatus, jqXHR) {
                $(".partCheckbox:checked").each(function(index,checkbox) {
                    var td = checkbox.parentElement;
                    var tr = td.parentElement;
                    parts_table.fnDeleteRow(tr,null,true); 
                });
            }
        });
    };

    var afterTableUpdated = function() {
    
        $("#selectAll").click(function() {
            var n = $("#selectAll:checked").length;
            if (n > 0) { toggleCheckboxes(true); }
            else { toggleCheckboxes(false); }

            setButtonStates(); // must come after check boxes are set
        });

        // confirm removing items from project
        $( "#dialog-confirm" ).dialog({
                resizable: false,
                height:140,
                modal: true,
                autoOpen: false,
                buttons: {
                    Cancel: function() {
                        $( this ).dialog( "close" );
                    },
                    "Remove": function() {
                        removeSelectedProjectParts();
                        $( this ).dialog( "close" );
                    }
                }
        });

        $("#removeButton").attr("disabled","true");
        $("#removeButton").click(function(e) {
            var items = $(".partCheckbox:checked");
            if (items.length !== 0) {
                $("#dialog-confirm").dialog("open");
            }
        });

        // activate/inactivate actions buttons depending on whats selected
        $(".partCheckbox").click(function() {
            setButtonStates();
        });
    }; // end of afterTableUpdated()



    // datatable must come before dialog-confirm because trigger is on select box
    parts_table.dataTable({
        "bAutoWidth": false,
        "bStateSave": true,
        "bJQueryUI": true,
        "aoColumns": [ null, null, null ],
        "oLanguage" : { "sUrl": "/res/tgi.txt" },
        "fnInitComplete": function() {
            $("#parts_table_length").append('<br/><div class="action_box"><input id="selectAll" type="checkbox"/> Select All &nbsp; <button id="removeButton" class="actionButton">remove</button></div>');
            afterTableUpdated();
            $("select[name='parts_table_length']").change(afterTableUpdated);
        }
    });

});





