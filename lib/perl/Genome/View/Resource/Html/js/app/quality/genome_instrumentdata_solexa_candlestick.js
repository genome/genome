function render_candlestick_graph (data, v_width, v_height) {
    /* parameters */
    var vis_width = v_width,
        vis_height = v_height,
        vis_top = 9.5,
        vis_bottom = 20,
        vis_left = 30,
        vis_right = 30,
        quartile_width = 6,
        whisker_width = 1.5;

    var vis = new pv.Panel();
    vis.width(vis_width)
    .height(vis_height)
    .top(vis_top)
    .bottom(vis_bottom)
    .left(vis_left)
    .right(vis_right);

    /* scales */
    var xscale = pv.Scale.linear(data, function(d) { return Math.min(d.column); }, function(d) { return Math.max(d.column); }).range(5, vis_width),
        yscale = pv.Scale.linear(-15, 45).range(0, vis_height);

    /* y axis */
    vis.add(pv.Rule)
    .data(yscale.ticks())
    .bottom(yscale)
    .left(-5).right(-5)
    .strokeStyle("#f1f1f1")
    .anchor("left").add(pv.Label);

    /* x axis */
    vis.add(pv.Rule)
    .data([1].concat(xscale.ticks()))
    .strokeStyle("#f1f1f1")
    .left(xscale)
    .anchor("bottom").add(pv.Label);

    /* candlestick type graph */
    var tooltip_txt = function(d) { return "Column: " + d.column + "  Quality Min: " + d.quality_min + "  Quality Max: " + d.quality_max + "  Q1: " + d.quartile_q1 + "  Q3: " + d.quartile_q3 + "  Qmed: " + d.quartile_med + "  IQR: " + d.quartile_iqr };

    vis.add(pv.Rule) // whisker line
     .data(data)
     .left(function(d) { return xscale(d.column); })
     .bottom(function(d) { return yscale(d.whisker_left); })
     .height(function(d) { return yscale(d.whisker_right) - yscale(d.whisker_left); })
     .strokeStyle("#d2d0a5")
     .lineWidth(whisker_width)
     .title(tooltip_txt)
    .add(pv.Rule) // whisker bottom bracket
     .bottom(function(d) { return yscale(d.whisker_left); })
     .height(whisker_width)
     .lineWidth(5)
     .title(tooltip_txt)
     .add(pv.Rule) // whisker top bracket
     .bottom(function(d) { return yscale(d.whisker_right); })
     .height(whisker_width)
     .lineWidth(5)
     .title(tooltip_txt)
    .add(pv.Rule) // quartile
     .bottom(function(d) { return yscale(d.quartile_q1); })
     .height(function(d) { return yscale(d.quartile_q3) - yscale(d.quartile_q1); })
     .lineWidth(quartile_width)
     .strokeStyle("#60604b")
     .title(tooltip_txt)
    .add(pv.Rule) // quartile med
     .bottom(function(d) { return yscale(d.quartile_med); })
     .height(whisker_width)
     .lineWidth(quartile_width)
     .strokeStyle("#f3b028")
     .title(tooltip_txt);

    vis.render();

}
