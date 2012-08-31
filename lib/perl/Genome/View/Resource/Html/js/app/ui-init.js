
$('document').ready(function() {

        $.ajaxSetup({
            cache: true,
            beforeSend: function(xhr) {
                if ($(document).data('updatedOn')) {
                    var d = $(document).data('updatedOn').getTime();
                    var elapsed = Math.ceil((Date.now() - d) / 1000);
                    xhr.setRequestHeader("X-Max-Age", elapsed.toString());
                    return true;
                }
            }
        });

        // set user name span at top of page
        $.ajax({
            url: "/view/whoami",
            dataType: "json",
            type: "GET",
            context: document.body,
            success: function(data, textStatus, jqXHR) {
                $("#authUser").html(data.username);
                $("#authUserLink").attr("href","/view/genome/sys/user/status.html?id=" + data.username + '@genome.wustl.edu');
                $("#authUser").trigger('change');
            },
            error: function(jqXHR, textStatus, errorThrown) {
                alert("ajax query to get username failed");
            }
        });

        // apply jQueryUI styles to button elements
        $("a.button, input:submit, button").button();

        // draw last updated time

        if (document.cookie.indexOf("cacheon=1") >= 0) {
            $('#updatedOn').text($(document).data('updatedOn').toString()).easydate();

            $('#refreshCache').click(function() {
                var url = location.pathname.substr(5) + location.search;

                $.ajax({
                    url: '/cachetrigger' + url,
                    success: function(data) {
                        location.reload();
                    }
                });

                $(this).parent().parent().find('.cache_time p').replaceWith("<p style='margin-top: 12px;'><strong>Loading...</strong></p>");

                return false;
            });
        } else {
            $('.cache_info').hide();
        }

        if ($('#objects').length) {
        // init masonry for view object container
            $('#objects').masonry(
                {
                    columnWidth: $('.span_8_box_masonry').width()+10,
                    singleMode: true,
                    itemSelector: '.span_8_box_masonry'
                }
            );

            $('#objects').masonry(
                {
                    columnWidth: $('.span_12_box_masonry').width()+10,
                    singleMode: true,
                    itemSelector: '.span_12_box_masonry'
                }
            );
        }

        // set up tasks popup window & button
        $('a#view_all_tasks').click(function() {
            var ptitle = this.title;
            var popup = $("#tasks_table");

            popup.dialog({
                title: ptitle,
                width: 750,
                height: 300
            });

            var dWidget = popup.dialog("widget");

            // tweak titlebar styles
            dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

            return false;
        });

        //perspective switcher
        $('a#perspective_switcher').click(function() {
            var url = this.href;
            var ptitle = this.title;
            var popup = $('<div id="perspective_switcher_popup"><p style="text-align: center; margin-top: 50px;"><img src="/res/img/spinner_ffffff.gif" width="16" height="16" align="absmiddle"/> Loading ' + ptitle + ' ...</p></div>').appendTo('body');

            popup.dialog({
                title: ptitle,
                width: 550,
                height: 450
            });

            var dWidget = popup.dialog("widget");
            dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");

            $.get(url, function(response, status, xhr) {
                var r = $(response).find("ul");
                r.find("a").each(function() {
                    $(this).attr('href', ($(this).attr('href') + location.search));
                });
                popup.html(r);
            });

            return false;
        });
/*
        // set up control bar state & behavior
        $('#bar_menu, #bar_menu ul').hide();

        var barClosed = 1;

        $('#bar_base').mouseenter(function() {
            if (barClosed) {
                $('#bar_menu')
                    .show('fast', function() {
                        barClosed = 0;
                        $('#bar_menu ul').fadeIn('fast');
                    })
                    .mouseleave(function() {
                        var mouseBackOver = 0;
                        $(this).mouseenter(function(){ mouseBackOver = 1; });

                        // wait for a second to see if user hovers over menu again
                        setTimeout(function() {
                            if (!mouseBackOver) {
                            $('#bar_menu ul')
                                .fadeOut('fast', function(){
                                    $(this).parent()
                                        .hide('fast', function() {
                                            barClosed = 1;
                                        })
                                        .unbind();
                                });
                            }
                        }, 1000);
                    });
            }
        });
*/

});

