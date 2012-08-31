$(function (){
      $('a.popup-ajax').click(function() {
          var url = this.href;
          var ptitle = this.title;
          var popup = $('<div id="event_popup"><p style="text-align: center; margin-top: 50px;"><img src="/res/img/spinner_ffffff.gif" width="16" height="16" align="absmiddle"/> Loading ' + ptitle + ' ...</p></div>').appendTo('body');

          popup.dialog({
              title: ptitle,
              width: 550,
              height: 450
          });

          var dWidget = popup.dialog("widget");

          // tweak titlebar styles
          dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

          popup.load(url);

          return false;
      });
});