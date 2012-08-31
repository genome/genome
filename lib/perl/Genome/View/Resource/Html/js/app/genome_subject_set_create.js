$(document).ready(function() {

    $('#submit-button').click(function() {

        // validate
        var fakeFilename = $('#file').val(); 
        if ( ! fakeFilename ) {
            alert('Please select the file containing the subject data');
            return false;
        }

        var filename = fakeFilename.substring(1 + fakeFilename.lastIndexOf("\\"));

        if (filename.indexOf(".csv") > -1) {
            $('#filetype').val("csv");
        } else {
            alert("CSV (comma separated values) is the only format currently supported");
            return false;
        }

        if ( $("input[name='subclass_name']:checked").length <= 0 ) {
            alert('Please select sample or patient');
            return false;
        }

        if (! $('#nomenclature').val() ) {
            alert("Please specify a nomenclature");
            return false;
        }

        if (! $('#project').val() ) {
            alert("Please specify a project name");
            return false;
        }

        $('create-samples').submit();
    });

    if (window.location.hash) {
      m = window.location.hash.split('=');
      if (m[0] === '#nomenclature_name') {
        $("#nomenclature").val(m[1])
      }
    }



    var cache = [];
    // Arguments are image paths relative to the current page.
    $.preLoadImages = function() {
    var args_len = arguments.length;
        for (var i = args_len; i--;) {
            var cacheImage = document.createElement('img');
            cacheImage.src = arguments[i];
            cache.push(cacheImage);
        }
    }

    // nomenclature autocomplete stuff
    nAutocompleteUrl = "/view/genome/nomenclature/set/autocomplete.json"
    $.ajax({
        url: nAutocompleteUrl,
        type: "GET",
        context: document.body,
        success: function(data, textStatus, jqXHR) {
            var nAutoCompleteData = data.names;
            $("#nomenclature").autocomplete({
                source: nAutoCompleteData
            });
        },
        error: function() {
            alert("ajax request for nomenclature names failed");
        }
    });

});

