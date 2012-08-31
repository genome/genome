
$('document').ready(function() {
    $('#authUser').change(function() {
        handleProjectBox($(this).html());
    });
});

function handleProjectBox(userName) {
    var duration = 100;

/*
    $("a:not(#createProjectLink)").click(function() {
        $(".search_result").draggable("option", "disabled", true);
    });
*/

    updateProjectBox(userName); // see genome_projectbox.js
    setInterval(function(){ updateProjectBox(userName); }, 600000); // 10 minutes

    $(".search_result").draggable({
        handle: ".icon_grip",
        revert: "invalid",
        appendTo: "body",
        helper: "clone",
    });


// The stuff below is about the "create project" feature

    $("#createProjectDiv").hide();

    $("#createProjectLink").click(function() {
        $("#createProjectDiv").show(duration);
        $("#createProjectLink").hide(duration);
        $("#createProjectName").focus();
        $("#createProjectName").val('');
    });

    $("#createProject").submit(function() {

        $("#errorAccordion").hide();
        $('#loadingStatus').text("updating...");
        $("#submitCreateProject").attr("disabled","disabled");
        $("#createProjectDiv").hide(duration);

        var formData = $(this).serializeArray();
        // console.log(formData);

        $.ajax({
            url: "/view/genome/project",
            type: "PUT",
            data: formData,
            context: document.body,
            success: function(data, textStatus, jqXHR) {

                // reset state of form
                $("#submitCreateProject").removeAttr("disabled");
                $("#createProject > input").map( 
                        function(i,e) { 
                            if(e.name != '') { 
                                e.value = ''; 
                            } 
                        });
           
                var d = $.parseJSON(data); 
                if(d.id == undefined) {
                    $('#loadingStatus').text('failed.');
                    $('#errorMsg').html(d.error);
                    var e = $("#errorAccordion");
                    e.show(); e.accordion("option","active",false);
                    $('html,body').animate({ scrollTop: $(document.body).offset().top }, 200);
                }

                updateProjectBox(userName);
            },
            error: function(jqXHR, textStatus, errorThrown) {
                alert("failed to create project (ajax request failed)");
                $("#submitCreateProject").removeAttr("disabled");
                $("#createProjectLink").show(duration);
            }
        });

        return false;    
    });

}



