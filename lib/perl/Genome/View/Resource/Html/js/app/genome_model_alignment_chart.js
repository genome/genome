// INCLUDED FROM:
// xsl/html/coverage/genome_model.xsl
function render_alignment_chart (aSummary) {

    aSummary.sort(aSummarySort);

    var metrics = [
        "unique_target_aligned_bp",
        "duplicate_target_aligned_bp",
        "unique_off_target_aligned_bp",
        "unique_off_target_aligned_bp_500",
        "duplicate_off_target_aligned_bp",
        "total_unaligned_bp"
    ];

    var metrics_short = [
        "unique on target",
        "duplicate on target",
        "unique off target (wingspan 500)",
        "unique off target",
        "duplicate off target",
        "unaligned"
    ];

    var c = pv.colors("#1f77b4", "#aec7e8", "#ff7f0e", "#fcb34c", "#fcd297", "#CCCCCC" );

    // create column arrays
    var summary_data = [];
    for (var subject in aSummary) {
        var summary_col = [];
        for (var i in metrics) {
            summary_col.push(aSummary[subject][metrics[i]]);
        }
        summary_data.push(summary_col);
    }

    // determine max to calculate the width of the chart
    var max = pv.max(aSummary, function(d) { return d.total_bp });

    // protovis' Stack layout likes the data w/ one array per layer instead of one per column,
    // therefore we must transpose the conversion results
    var summary_data_t = pv.transpose(summary_data);

    var alignment_w = 290,
    alignment_h = 16 * aSummary.length,
    alignment_x = pv.Scale.linear(0, max).range(0, alignment_w-10),
    alignment_y = pv.Scale.ordinal(pv.range(aSummary.length)).splitBanded(0, alignment_h, .90);

    // adjust the color palette
    // c = pv.Colors.category20();

    var alignment_vis = new pv.Panel()
        .width(alignment_w)
        .height(alignment_h)
        .bottom(0)
        .left(0)
        .right(10)
        .top(110);

    var bar = alignment_vis.add(pv.Layout.Stack)
        .layers(summary_data_t)
        .orient("left-top")
        .x(function() { return alignment_y(this.index) })
        .y(alignment_x)
        .layer.add(pv.Bar)
        .height(alignment_y.range().band)
        .fillStyle(function(d) { return c(this.parent.index) })
        .title(function(d) { return "model: " + aSummary[this.index].model_name  + "; subject: " + aSummary[this.index].subject_name + "; " + metrics_short[this.parent.index] + ": " + addCommas(d) } );

    bar.anchor("right").add(pv.Label)
        .visible(function() { return this.parent.index == 0 }) // only show a label when the bar is wide enough
        .textStyle("white")
        .text(function(d) { return round(d/1000000000,2).toFixed(2); });

    // bar.anchor("left").add(pv.Label)
    //     .visible(function() !this.parent.index)
    //     .textMargin(5)
    //     .textAlign("right")
    //     .text(function() aSummary[this.index].subject_name );

    // add x axis rules
    alignment_vis.add(pv.Rule)
        .data(alignment_x.ticks(5))
        .left(alignment_x)
        .strokeStyle(function(d) { return (d == 0 || d == max) ? "#AAA" : "rgba(255,255,255,.3)" })
        .add(pv.Rule)
        .top(0)
        .height(5)
        .strokeStyle("rgba(255,255,255,.3)")
        .anchor("top").add(pv.Label)
        .text(function(d) { return d == 0 ? "" : d/1000000000 });

    // add target rule
    alignment_vis.add(pv.Rule)
        .data([6000000000])
        .left(alignment_x)
        .strokeStyle("#F00")
        .add(pv.Rule)
        .top(-10)
        .height(15)
        .strokeStyle("#F00")
        .anchor("top").add(pv.Label)
        .text(function(d) { return d/1000000000 });

    // legend

    alignment_vis.add(pv.Panel)
        .top(-100)
        .left(5)
        .add(pv.Dot)
        .data(metrics_short)
        .top(function() { return this.index * 15 })
        .size(8)
        .shape("square")
        .strokeStyle(null)
        .fillStyle(function(d) { return c(this.index) })
        .anchor("right").add(pv.Label);

    // x axis label
    alignment_vis.add(pv.Label)
        .left(115)
        .font("bold 14px sans-serif")
        .top(-25)
        .text("sequence (Gb)");

    alignment_vis.render();

    function addCommas(nStr) {
        nStr += '';
        var x = nStr.split('.');
        var x1 = x[0];
        var x2 = x.length > 1 ? '.' + x[1] : '';
        var rgx = /(\d+)(\d{3})/;
        while (rgx.test(x1)) {
            x1 = x1.replace(rgx, '$1' + ',' + '$2');
        }

        return x1 + x2;
    }

    function round(rnum, rlength) {
        return Math.round(rnum*Math.pow(10,rlength))/Math.pow(10,rlength);
    }

    function aSummarySort(a,b) {
         // if (a.subject_name < b.subject_name) { return -1; };
        // if (a.subject_name > b.subject_name) { return 1; };

        if (a.model_name < b.model_name) { return -1; };
        if (a.model_name > b.model_name) { return 1; };

        if (a.lane_count < b.lane_count) { return -1; };
        if (a.lane_count > b.lane_count) { return 1; };

        if (a.id < b.id) { return -1; };
        if (a.id > b.id) { return 1; };

        return 0;
    }
}
