$(document).ready(function() {

    var flagstat = $("div#flagstat_table").dialog({
        title: "Flagstat Report",
        width: 500,
        height: 450,
        autoOpen: false
    });

    // tweak titlebar styles
    var dWidget = flagstat.dialog("widget");
    dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

    $('a#flagstat_button').click(
        function() {
            $('div#flagstat_table').dialog('open');
            return false;
        }
    );

    $('a.popup-ajax-input-diff').click(function() {
        var ptitle = this.title;
        var popup_selector = '#input_diff_' + this.id;
        var popup = $("#input_diff_" + this.id);

        popup.dialog({
            title: ptitle,
            width: 450,
            height: 250
        });

        var dWidget = popup.dialog("widget");

        // tweak titlebar styles
        dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

        return false;
     });

      $('a.notes-popup').click(function() {
          var ptitle = this.title;
          var popup_selector = '#notes_subject_' + this.id;
          var popup = $("#notes_subject_" + this.id);

          popup.dialog({
              title: ptitle,
              width: 450,
              height: 250
          });

          var dWidget = popup.dialog("widget");

          // tweak titlebar styles
          dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

          return false;
      });

});