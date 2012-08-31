// INCLUDED FROM:
// xsl/html/coverage/genome_model.xsl
function render_enrichment_chart (raw_enrichment_data) {

    raw_enrichment_data.sort(enrichmentSort);

    /* get chart data into nice arrays */
    var models = new Array();
    var factors = pv.keys(raw_enrichment_data[0]['enrichment_factors']);

    //convert full % to stacked %
    var enrichment_data = []; // this will store the stacked version of the enrichment data
    var enrichment_data_full = []; // we'll need this to show the full instead of the stacked % in rollovers

    for (var model in raw_enrichment_data) {
        models.push(raw_enrichment_data[model].subject_name);
        
        var efactors    = raw_enrichment_data[model].enrichment_factors;
        var unique      = Number(efactors['unique_on_target']);
        var total       = Number(efactors['total_on_target']);
        var theoretical = Number(efactors['theoretical_max']);

        var st_unique = unique;
        var st_total = round(total - unique, 1);
        var st_theo = round(theoretical - total, 1);

        var full_factor_pc = [unique, total, theoretical];
        var stacked_factor_pc = [st_unique, st_total, st_theo];

        enrichment_data.push(stacked_factor_pc);
        enrichment_data_full.push(full_factor_pc);

    }

    // protovis' Stack layout likes the data w/ one array per layer instead of one per column,
    // therefore we must transpose the conversion results
    var enrichment_data_t = pv.transpose(enrichment_data);

    // determine max to calculate the width of the chart
    var max = pv.max(raw_enrichment_data, function(d) { return d.enrichment_factors.theoretical_max; });

    var enrichment_w = 180,
        enrichment_h = 16 * models.length,
        enrichment_x = pv.Scale.linear(0, max).range(0, enrichment_w-10),
        enrichment_y = pv.Scale.ordinal(pv.range(models.length)).splitBanded(0, enrichment_h, .90);

    var c = pv.colors("#865c21", "#f4ac44", "#866d4c");

    var enrichment_vis = new pv.Panel()
        .width(enrichment_w)
        .height(enrichment_h)
        .bottom(0)
        .left(0)
        .right(0)
        .top(110);

    var bar = enrichment_vis.add(pv.Layout.Stack)
        .layers(enrichment_data)
        .orient("left-top")
        .x(function() { return enrichment_y(this.index); } )
        .y(enrichment_x)
        .layer.add(pv.Bar)
        .height(enrichment_y.range().band)
        .fillStyle(function(d) { return c(this.parent.index); })
        .title(function(d) { return enrichment_data_full[this.index][this.parent.index].toFixed(1); });        

    bar.anchor("right").add(pv.Label)
        .visible(function(d) { return this.parent.index == 1 ? false : true; }) // only show a label on first and last bars.
        .textStyle("white")
        .text(function(d) { return enrichment_data_full[this.index][this.parent.index].toFixed(1); });

    // legend
    enrichment_vis.add(pv.Panel)
        .top(-100)
        .left(4)
        .add(pv.Dot)
        .data(factors)
        .top(function() { return this.index * 15; } )
        .size(8)
        .shape("square")
        .strokeStyle(null)
        .fillStyle(function(d) { return c(this.index); })
        .anchor("right").add(pv.Label)
        .text(function(d) { return d.replace(/_/g, " "); });

    // x axis label
    enrichment_vis.add(pv.Label)
        .left(18)
        .font("bold 14px sans-serif")
        .top(-25)
        .text("enrichment factor");

    enrichment_vis.render();

    function round(rnum, rlength) {
        return Math.round(rnum*Math.pow(10,rlength))/Math.pow(10,rlength);
    };

    function enrichmentSort(a,b) {
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
