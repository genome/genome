function render_bar_graph (data, v_width, v_height) {
    /* parameters */
    var vis_width = v_width,
        vis_height = v_height,
        vis_top = 9.5,
        vis_bottom = 20,
        vis_left = 30,
        vis_right = 30,
        bar_width = .86;

    var vis = new pv.Panel();
    vis.width(vis_width)
    .height(vis_height)
    .top(vis_top)
    .bottom(vis_bottom)
    .left(vis_left)
    .right(vis_right);

    /* scales */
    var xscale = pv.Scale.linear(data, function(d) { return Math.min(d.column); }, function(d) { return Math.max(d.column); }).range(5, vis_width),
        yscale = pv.Scale.linear(data, function(d) { return 0; }, function(d) { return d.count; }).range(0, vis_height);

    /* y axis */
    vis.add(pv.Rule)
    .data(yscale.ticks())
    .bottom(yscale)
    .left(-5).right(-5)
    .strokeStyle("#f1f1f1")
    .anchor("left")
        .left(8)
        .add(pv.Label)
            .text(function(d) { return d.toExponential(); });

    /* x axis */
    vis.add(pv.Rule)
    .data([1].concat(xscale.ticks()))
    .strokeStyle("#f1f1f1")
    .left(xscale)
    .anchor("bottom").add(pv.Label);

	 /* stacked bar graph */
    var tooltip_txt = function(d) { return "Column: " + d.column + "  Count: " + d.count + "  A: " + d.count_a + "  C: " + d.count_c + "  G: " + d.count_g + "  T: " + d.count_t + "  N: " + d.count_n; };

    vis.add(pv.Panel)
    .data(["count_a", "count_c", "count_g", "count_t", "count_n"])
    .add(pv.Bar)
    .data(data)
    .width(xscale(bar_width))
    .bottom(pv.Layout.stack())
    .height(function(d, t) { return yscale(d[t]); })
    .left(function(d) { return xscale(d.column); })
    .title(tooltip_txt)
	 .fillStyle(function(d,t) {
		  var bcolor;
		  switch(t) {
		  case "count_a":
				bcolor = "#7d598c";
				break;
		  case "count_c":
				bcolor = "#bb4250";
				break;
		  case "count_g":
				bcolor = "#90c86f";
				break;
		  case "count_t":
				bcolor = "#f6c460";
				break;
		  case "count_n":
				bcolor = "#999";
				break;
		  }
		  return bcolor;
	 });

    vis.render();
}
