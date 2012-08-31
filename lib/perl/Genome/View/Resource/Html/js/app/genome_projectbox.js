

$(document).ready(function(){
    $( "#errorAccordion" ).accordion({
        collapsible: true,
             active: false,
    });
    $("#errorAccordion").hide();
});


function assignProjectPart(droppableThing, draggableThing) {

$('#loadingStatus').text("updating...");
$('#errorMsg').html('');
$('#errorAccordion').hide();

var objParts = draggableThing.children('div .icon_grip').attr('object_id').split('---');
var projectPartClass = objParts[0];
var projectPartId = objParts[1];
var projectId = droppableThing.getAttribute('project_id');
var url = "/view/genome/project-part";

    $.ajax({
        url: "/view/genome/project-part",
        type: "PUT",
        data: {
            'entity_class_name': projectPartClass,
            'entity_id': projectPartId,
            'project_id': projectId,
            'label': 'web'
        },
        context: document.body,
        success: function(data, textStatus, jqXHR) {

            updateProjectBox();
            var d = $.parseJSON(data);

            if(d.id == undefined) {
                $('#errorMsg').html(d.error);
                var e = $("#errorAccordion");
                e.show();
                e.accordion("option","active",false);
                $('html,body').animate({ scrollTop: $(document.body).offset().top }, 200);
                updateProjectBox();
            }
        },
        error: function(jqXHR, textStatus, errorThrown) {
            $('#errorMsg').html("Failed to assign project part- ajax request failed");
            var e = $("#errorAccordion");
            e.show();
            e.accordion("option","active",false);
            $('html,body').animate({ scrollTop: $(document.body).offset().top }, 200);
            updateProjectBox();
        }
    });

}

function memberSort(a, b) {

    a = a.name.toUpperCase();
    b = b.name.toUpperCase();

    if (a == b) {
        return 0;
    } else if(a > b) {
        return 1;
    } else if(a < b) {
        return -1;
    }
}

function updateProjectBox(userName) {

    var duration = 100;
    $('#loadingStatus').text("updating...");
    var label = userName + "'s";
    if (! userName || userName == '-' ) {
        userName = $("#authUser").html();
        label = "My";
    }

    if (! userName || userName == '-' ) {
        alert("Failed to get userame");
        return false;
    }

    var projectsURL = "/view/genome/project/set/status.json?user_ids=" + userName + '@genome.wustl.edu';
    $.ajax({
        url: projectsURL,
        type: "GET",
        context: document.body,
        success: function(data, textStatus, jqXHR) {

            $('#projectBox li').remove();
            if (data.members != undefined && data.members.length > 0) {

                var href = '<a href="/view/genome/sys/user/status.html?id=' + userName + '@genome.wustl.edu">';
                var projectBoxLabel = label + ' projects (' + href + 'see all ' + data.count + '</a>)';
                $('#myProjectsCount').html(projectBoxLabel);

                var sortedMembers = data.members.sort(memberSort);
                sortedMembers.map(function(e, i) {
                                    var count = e.parts_count || '?';
                                    $('#projectBox').append('<li><div class="projectContainer clickable" project_id="' 
                                                            + e.id + '">' 
                                                            + '<a href="/view/genome/project/status.html?id='
                                                            + e.id + '">' 
                                                            + e.fixed_size_name + '</a> (' + count + ')' 
                                                            + '</div></li>');
                                    $('.projectContainer').droppable({
                                        activeClass: "ui-state-default",
                                        hoverClass: "ui-state-hover",
                                        tolerance: "pointer",
                                        drop: function( event, ui ) {
                                            assignProjectPart(this, ui.draggable);
                                        }
                                    });

                                });
            } else {
                $('#myProjectsCount').text(label);
                $('#projectBox').append('<li><div class="projectContainer">None (yet!)<div></li>');
            }
            $('#loadingStatus').text('');
            $("#createProjectLink").show(duration);
        },
        error: function(e) {
           $('#loadingStatus').text('failed.');
           $("#createProjectLink").show(duration);
        }
    });
}



