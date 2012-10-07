$(document).ready( function() {
    var projectPartUrl = "/view/genome/project-part";

    var labelWatchButton = function() {
        $.ajax({
            "dataType": "json",
            "type": "GET",
            "url": projectPartUrl + "/default.json?role=watcher&entity_id=" + $("#hiddenAuthUserID").html(),
            "error" : function (jqXHR, textStatus, errorThrown) {
                console.log(errorThrown);
                if (jqXHR.status == 404) {
                    $("#watchButton").removeAttr("disabled");
                    $("#watchButton").button("option", "label", "Watch");
                } else {
                    alert("Failed to check watcher status.");
                    console.log(errorThrown);
                }
            },
            "success" : function(data, textStatus, jqXHR) {
                $("#watchButton").removeAttr("disabled");
                $("#watchButton").button("option", "label", "Unwatch");
            }
        });
    };

    var toggleWatch = function() {
        var entityID = $("#hiddenAuthUserID").html();
        var entityClassName = 'Genome::Sys::User';
        var role = 'watcher';
        var projectID = $("#projectID").html();
        var label = $("#watchButton").button("option", "label");
        var data = {
            'entity_class_name': entityClassName,
            'entity_id': entityID,
            'role': role,
            'project_id': projectID,
        };
        var type;
        if (label == "Unwatch") {
            type = 'DELETE';
        } else if (label == "Watch") {
            type = 'PUT';
        } else {
            console.log("Watcher status (" + label + ") unknown, cannot toggle.");
        }
        $.ajax({
            url: projectPartUrl,
            type: type,
            data: data,
            context: document.body,
            success: function(data, textStatus, jqXHR) {
                location.reload();
            },
            error: function(jqXHR, textStatus, errorThrown) {
                var verb = type == 'DELETE' ? 'remove' : 'add';
                alert("Failed to " + verb + " you as a watcher.");
                console.log(errorThrown);
            }
        });
    }

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

        $.ajax({
            "dataType": "json",
            "type": "DELETE",
            "url": projectPartUrl,
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

        $("#watchButton").button();
        $("#watchButton").on('click', function(e) {
            $("#watchButton").attr("disabled", "disabled");
            toggleWatch();
            labelWatchButton();
        });
        $("#hiddenAuthUserID").on('change', labelWatchButton);

    }; // end of afterTableUpdated()



    // datatable must come before dialog-confirm because trigger is on select box
    parts_table.dataTable({
        "bAutoWidth": false,
        "bStateSave": true,
        "bJQueryUI": true,
        "aoColumns": [ null, null, null ],
        "oLanguage" : { "sUrl": "/res/tgi.txt" },
        "fnInitComplete": function() {
            $("#parts_table_filter").append('<br/><button id="watchButton">...</button>');
            $("#parts_table_length").append('<br/><div class="action_box"><input id="selectAll" type="checkbox"/> Select All &nbsp; <button id="removeButton" class="actionButton">remove</button></div>');
            afterTableUpdated();
            $("select[name='parts_table_length']").change(afterTableUpdated);
        }
    });

});





