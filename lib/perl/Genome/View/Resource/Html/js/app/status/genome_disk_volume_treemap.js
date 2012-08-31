function render_treemap(data, width, height, build_info, total_allocation) {
    var legend_height = 30;
    var allocation_meter_height = 30;

    // the following two arrays contain the same info as the object that follows them,
    // but I coldn't figure out how to do what needed to be done w/ just the object.
    // Hence should probably be refactored b/c it just looks inelegant...

    var statuses = [
        "No Build",
        "Succeeded",
        "Done",
        "Running",
        "Scheduled",
        "New",
        "Failed",
        "Crashed",
        "Abandoned"
    ];

    var status_colors = [
        "#DDDDDD",
        "#d9ffc7",
        "#d9ffc7",
        "#f8f7a3",
        "#ffac2a",
        "#ffac2a",
        "#ff9292",
        "#ff9292",
        "#666666"
    ];

    var status_colors_h = {
        "--": "#DDD",
        "Succeeded": "#d9ffc7",
        "Done": "#d9ffc7",
        "Running": "#f8f7a3",
        "Scheduled": "#ffac2a",
        "New": "#ffac2a",
        "Failed": "#ff9292",
        "Crashed": "#ff9292",
        "Abandoned": "#666666"
    };

    var total_allocated = [
        total_allocation['total_kb'] - total_allocation['unallocated_kb']
    ];

    // group allocation data by owner_class_name
    var ndata = pv.nest(data)
        .key(function(d) { return d.owner_class_name; })
        .key(function(d) { return d.display_name; })
        .rollup(function(d) { return Number(d[0].kilobytes_requested); });

    // group allocation data by allocation id
    var adata = pv.nest(data)
        .key(function(d) { return d.id; })
        .map();

    // console.dir(window.allocation_build_info);
    // console.dir(data);
    // console.dir(ndata);
    // console.dir(window.allocations);
    // console.dir(adata);
    var nodes = pv.dom(ndata).root("allocations").nodes();

    var vis = new pv.Panel()
        .width(width)
        .height(height)
        .fillStyle("lightyellow");

    var treemap = vis.add(pv.Layout.Treemap)
        .top(allocation_meter_height)
        .height(height - legend_height - allocation_meter_height)
        .width(width)
        .nodes(nodes)
        .round(true)
        .mode("squarify")
        .strokeStyle("#999");

    var leafTipsy = pv.Behavior.tipsy({gravity: "w", fade: true});

    treemap.leaf.add(pv.Panel)
        .def("hover", -1)
        .fillStyle(function(d) { return this.hover() == this.index ? "white" : status_colors_h[build_info[d.nodeName].status]; })
        .strokeStyle("#999")
        .lineWidth(1)
        .antialias(false)
        .title(function(d) { return "ID: " + d.nodeName; })
        .cursor("pointer")
        .event("mouseover", function() { return this.hover(this.index); } )
        .event("mouseout", function() { return this.hover(-1); })
        .event("click", function(d) { renderPopup(d); });

    treemap.label.add(pv.Label)
        .visible(function(d) { return getVisible(d); })
        .font(function(d) { return getFont(d); })
        .textStyle(function(d) { return getTextStyle(d); });

    var xalloc = pv.Scale.linear(0, total_allocation['total_kb']).range(0, width);

    var allocation_meter = vis.add(pv.Panel)
        .width(width)
        .height(allocation_meter_height)
        .fillStyle("#e5e5e3")
        .anchor("right").add(pv.Label)
            .text("")
            .textStyle("#333");

    allocation_meter.add(pv.Bar)
        .data(total_allocated)
        .top(0)
        .left(0)
        .height(30)
        .width(function(d)  { return xalloc(d); })
        .fillStyle(function(d) { 
            var perc_allocated = Math.round((d/total_allocation['total_kb'])*100);
            // console.log("perc_allocated: " + perc_allocated);
            var fillColor;
            if (perc_allocated < 80) {
                fillColor = "#878775";
            } else if (80 <= perc_allocated < 90) {
                fillColor = "#fdba05";
            } else if (90 <= perc_allocated <= 100) {
                fillColor = "#bf3228";                
            }
            return fillColor;
        })
        .anchor("right").add(pv.Label)
            .text(function(d) { return "allocated (" + Math.round((d/total_allocation['total_kb'])*100) + "%)"; })
            .textStyle("#FFF");

    var legend = vis.add(pv.Panel)
        .width(width)
        .bottom(0)
        .height(legend_height)
        .fillStyle("#FFFFFF");

    legend.add(pv.Dot)
        .data(statuses)
        .left(function() { return 15 + this.index * 80; })
        .size(30)
        .top(15)
        .shape("square")
        .fillStyle(function() { return status_colors[this.index]; })
        .strokeStyle("#999")
        .lineWidth(1)
        .antialias(false)
        .anchor("right").add(pv.Label)
            .text(function() { return statuses[this.index]; });

    vis.render();

    function renderPopup(d) {
        // console.group("renderPopup");
        var allocation_id = d.nodeName,
        popup_content = '<table class="boxy_info" cellpadding="0" cellspacing="0" border="0" width="300"><tbody>',
        build_obj = window.allocation_build_info[allocation_id],
        alloc_obj = adata[allocation_id][0];

        // console.log("alloc_obj: " + alloc_obj);
        // console.log("allocation_id: " + allocation_id);
        
        // console.dir(alloc_obj[0]);

        popup_content += '<tr><td colspan="2"><h3 class="subhead">allocation</h3></td></tr>';

        for (var key in alloc_obj) {
            popup_content += '<tr><td class="label" style="width: 50%">' + key + ':</td><td class="value" style="width: 50%">' + alloc_obj[key] + "</td></tr>";
        }

        popup_content += '<tr><td colspan="2"><h3 class="subhead">build</h3></td></tr>';

        for (var key in build_obj) {
            popup_content += '<tr><td class="label" style="width: 50%">' + key + ':</td><td class="value" style="width: 50%">' + build_obj[key] + "</td></tr>";
        }

        popup_content += '</tbody></table>';
        var popup = $('<div id="allocation_popup_' + allocation_id + '">' + popup_content + '</div>').appendTo('body');
        popup.dialog({
            title: "allocation " + allocation_id + " info",
            width: 500,
            height: 600
        });
        var dWidget = popup.dialog("widget");
        // tweak titlebar styles
        dWidget.find(".ui-dialog-titlebar").removeClass("ui-corner-all").addClass("ui-corner-top");
        // console.groupEnd();
        return false;
    }

    function getVisible(d) {
        // only show label if panel is large enough to contain it
        // and this is a child node
        var showLabel;
            if (d.parentNode && (d.dx > 50 || d.dy > 50)){
            showLabel = true;
        } else {
            showLabel = false;
        }
        
        return showLabel;
    }

    function getFont(d) {
        // show owner class names in a larger font
        var font;
        if ( d.depth == 1 ) {
            font = "bold 14px sans-serif";
        } else {
            font = "10px sans-serif";
        }
        return font;
    }

    function getTextStyle (d) {
        // give the owner class names a bit of an alpha transparency
        var textStyle;

        if ( d.depth == 1 ) {
            textStyle = "rgba(0,0,0,0.5";
        } else {
            textStyle = "rgba(0,0,0,1)";
        }
        return textStyle;        
    }

}

