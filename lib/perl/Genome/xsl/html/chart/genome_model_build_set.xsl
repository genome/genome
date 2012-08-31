<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_model_build_set" match="object[./types[./isa[@type='Genome::Model::Build::Set']]]">

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="set_header">
      <xsl:with-param name="display_name" select="'Query Results'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div class="set_query rounded span-24 last">
          <div class="padding10">
            <strong>Query: </strong> <xsl:value-of select="aspect[@name='rule_display']/value" />
          </div>
        </div>

        <xsl:call-template name="genome_model_build_set_chart"/>

        <hr class="space" style="height: 10px; margin: 0;"/>

        <div class="span-24 last">
          <table width="100%" cellpadding="0" cellspacing="0" border="0" id="set" class="dataTable">
            <thead>
              <xsl:apply-templates select="aspect[@name='members']" mode="set_header" />
            </thead>
            <tbody>
              <xsl:apply-templates select="aspect[@name='members']" mode="set_row" />
            </tbody>
          </table>
        </div>
        <xsl:apply-templates select="aspect[@name='members']/object[1]" mode="set_table_init" />
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="genome_model_build_set_chart">
    <script type="text/javascript" src="/res/js/pkg/protovis.js"></script>
    <script type="text/javascript">
    </script>


    <script type="text/javascript+protovis">
      <xsl:text disable-output-escaping="yes">
        <![CDATA[

function dropDown(data) {
    var options_metrics = '<option>Metric<\/option>'
    $.each(data, function(index,value){
      options_metrics += '<option value="' + value + '">' + value + '<\/option>'
    })
    $("select#metric").html(options_metrics)
}

function renderGraph(data,field) {

    /* Sizing and scales. */
    var w = 400,
        h = 200,
        x = pv.Scale.linear(0, data.length - 1).range(0, w),
        y = pv.Scale.linear(0, pv.max(data, function(d) d.y)).range(0, h);

    /* The root panel. */
    var vis = new pv.Panel()
        .width(w)
        .height(h)
        .bottom(100)
        .left(100)
        .right(10)
        .top(5);

    /* X-axis ticks. */
    vis.add(pv.Rule)
        .data(data)
        .left(function(d) x(this.index))
        .strokeStyle("#eee")
        .add(pv.Rule)
        .bottom(0)
        .height(0)
        .strokeStyle("#000")
        .anchor("bottom").add(pv.Label)
        .textAngle(Math.PI/2)
        .textAlign("left")
        .text(function(d) { return d.x } );
        //.bottom(-25)

    vis.add(pv.Label)
        .left(200)
        .bottom(-95)
        .textAlign("center")
        .text(field + " vs. build ID")

    /* Y-axis ticks. */
    vis.add(pv.Rule)
        .data(y.ticks())
        .bottom(y)
        .strokeStyle(function(d) d ? "#eee" : "#000")
        .anchor("left").add(pv.Label)
        .text(y.tickFormat);

    /* The line. */
    vis.add(pv.Line)
        .data(data)
        .left(function(d) x(this.index))
        .bottom(function(d) y(d.y))
        .lineWidth(3);

    vis.canvas('fig')
    vis.render();
}

var jsonData = {};
var jsonFields = [];

function renderData(data,field) {
         var chartData = data["members"].filter(function(m) {
             if (m.hasOwnProperty('metrics')) {
               return true
             }
         }).map(function(m) {
             for ( var metric in m['metrics'] ) {
               if ( m['metrics'][metric].name == field ) {
                 return { x: m['metrics'][metric].build_id, y: m['metrics'][metric].value }
               }
             }
         })
         renderGraph(chartData,field);
}

function redrawGraph(menuform) {
  selecteditem = menuform.metric.selectedIndex ;
  metric = menuform.metric.options[ selecteditem ].value
  renderData(jsonData,metric)
}

$(document).ready(function () {

    $.ajax({
        url: location.href.replace('.html','.json'),
        dataType: 'json',
        success: function(data) {
            // store our initial json data
            jsonData = data;
            // draw our default graph with default walltime metric
            renderData(data,'sample.time');
            // create a hash of metric names for use in dropdown menu
            var fields = {};
            for ( var member in data['members'] ) {
              if ( data['members'][member].hasOwnProperty('metrics') ) {
                $.each( data['members'][member]['metrics'].map(function(m) {
                          return m.name
                        }), function(index, value) {
                        fields[value] = 1
                      })
              }
            }

            jsonFields = []
            $.each(fields, function(index, value) {
                jsonFields.push(index)
              }
            )
            dropDown(jsonFields.sort());
        }
    })

    $('select#metric').change(function() {
      dropDown(jsonFields.sort());
    })
})
        ]]>
      </xsl:text>
    </script>

  <div id="fig"></div>
  <form id="menu" action="" name="menu">
  <select name="metric" id="metric" onchange="redrawGraph(this.form)">
    <option>Metric</option>
  </select>
  </form>

  </xsl:template>

</xsl:stylesheet>
