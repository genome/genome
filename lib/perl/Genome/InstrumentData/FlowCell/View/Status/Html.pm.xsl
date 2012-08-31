<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:template match="/">

    <html>
      <head>
        <title>Flow Cell <xsl:value-of select="//flow-cell/@id"/></title>

        <link rel="shortcut icon" href="/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png" />

        <link rel="stylesheet" href="/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen" />
        <link rel="stylesheet" href="/res/old/report_resources/apipe_dashboard/css/tablesorter.css" type="text/css" media="screen" />
        <script type="text/javascript" src="/res/old/report_resources/jquery/jquery.js"></script>
        <script type="text/javascript" src="/res/old/report_resources/jquery/jquery.tablesorter.min.js"></script>
        <script type="text/javascript">
          $(document).ready(function() {
          $("#lane_list").tablesorter({
          // sort on first column, ascending
          sortList: [[0,0]]
          });
          });
        </script>
        <script type="text/javascript" src="/res/old/report_resources/protovis-3.1/protovis-d3.1.js"></script>
        <!-- set up lane data and chart -->
        <script type="text/javascript">
          window.indexData = [];

          id = window.indexData;

          <xsl:for-each select="//lane">

            <xsl:variable name="count" select="position() - 1" />
            id[<xsl:value-of select="$count"/>] = {};
            id[<xsl:value-of select="$count"/>].percent = [
            <xsl:for-each select="index">
              <xsl:value-of select="percent"/>,
              </xsl:for-each>];

              id[<xsl:value-of select="$count"/>].sequence = [
              <xsl:for-each select="index">
                "<xsl:value-of select="sequence"/>",
                </xsl:for-each>];

          </xsl:for-each>

          <![CDATA[

// determine max value of index percentages
// so that all charts have the same x scale

window.pMax = 0;
console.log("id.length: " + id.length);
for (var i=0,len=id.length;i<len;i++) {
for (j in id[i].percent) {
if (id[i].percent[j] > window.pMax) {
window.pMax = id[i].percent[j];
}
}
}

window.chartMax = Math.ceil(window.pMax);
window.chartMin = 0;
          ]]>
        </script>

        <style type="text/css">
          div.lane_chart {
          width: 200px;
          float: left;
          margin-left: 15px;
          margin-right: 15px;
          }

          div.legend_block {
          width: 910px;
          float: left;
          background: #EFEFEF;
          padding: 5px;
          }

          div.legend_block div.label {
          font-weight: bold;
          float: left;
          width: 50px;
          }

          div.legend_block div.legend {
          float: left;
          width: 800px;
          }

          div.charts_block {
          float: left;
          width: 920px;
          margin-bottom: 20px;
          }

          p.axis_label_x {
          text-align: center;
          color: #AAA;
          font-size: 11px;
          }

        </style>

      </head>

      <body>
        <div class="container">
          <div class="background">
            <div class="page_header">
              <table cellpadding="0" cellspacing="0" border="0">
                <tr>
                  <td>
                    <a href="status.cgi" alt="Go to Search Page" title="Go to Search Page"><img src="/res/old/report_resources/apipe_dashboard/images/gc_header_logo2.png" width="44" height="45" align="absmiddle" /></a>
                  </td>
                  <td>
                    <h1>Analysis Reports v0.1</h1>
                  </td>
                </tr>
              </table>
            </div>

            <div class="page_padding">
              <h2 class="page_title icon_instrument_data">Flow Cell <xsl:value-of select="//flow-cell/@id"/> Status</h2>
              <table cellpadding="0" cellspacing="0" border="0" class="info_table_group">
                <tr>
                  <td>
                    <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                      <colgroup>
                        <col/>
                        <col width="100%"/>
                      </colgroup>
                      <tr><td class="label">Flow Cell ID:</td><td class="value"><xsl:value-of select="//flow-cell/@id"/></td></tr>

                      <tr><td class="label">Run Type:</td><td class="value"><xsl:value-of select="//production/@run-type"/></td></tr>
                      <tr><td class="label">Machine:</td><td class="value"><xsl:value-of select="//production/@machine-name"/></td></tr>
                    </table>
                  </td>
                  <td>
                    <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                      <colgroup>
                        <col/>
                        <col width="100%"/>
                      </colgroup>
                      <tr><td class="label">Run Name:</td><td class="value"><xsl:value-of select="//production/@run-name"/></td></tr>
                      <tr><td class="label">Date Started:</td><td class="value"><xsl:value-of select="//production/@date-started"/></td></tr>
                      <tr><td class="label">Group:</td><td class="value"><xsl:value-of select="//production/@group-name"/></td></tr>
                    </table>
                  </td>
                </tr>
              </table>
              <hr/>

              <h2>lanes</h2>
              <table id="lane_list" class="list tablesorter" width="100%" cellspacing="0" cellpadding="0" border="0">
                <colgroup>
                  <col />
                  <col />
                </colgroup>
                <thead>
                  <tr>
                    <th>lane</th>
                    <th>id</th>
                    <th>reports</th>
                    <th>resources</th>
                  </tr>
                </thead>
                <tbody>
                  <xsl:choose>
                    <xsl:when test="count(//instrument-data) > 0">
                      <xsl:for-each select="//instrument-data">
                        <xsl:sort select="@subset_name" data-type="text" order="ascending"/>
                        <tr>
                          <td>
                            <xsl:value-of select="@subset_name"/>
                          </td>
                          <td><xsl:value-of select="@id"/></td>
                          <td>
                            <xsl:choose>
                              <xsl:when test="report">
                                <xsl:for-each select="report">
                                  <xsl:variable select="@name" name="report_name_full"/>
                                  <a><xsl:attribute name="class">btn_link</xsl:attribute><xsl:attribute name="href">/cgi-bin/dashboard/flow_cell_report.cgi?instrument-data-id=<xsl:value-of select="../@id"/>&amp;report-name=<xsl:value-of select="@name"/></xsl:attribute><xsl:value-of select="substring-before($report_name_full,'.')"/></a><xsl:text> </xsl:text>
                                </xsl:for-each>
                              </xsl:when>
                              <xsl:otherwise>
                                <span class="note">No reports available for this lane.</span>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td>
                            <xsl:choose>
                              <xsl:when test="@gerald-directory">
                                <a><xsl:attribute name="class">btn_link</xsl:attribute><xsl:attribute name="href">https://gscweb<xsl:value-of select="@gerald-directory"/></xsl:attribute>gerald directory</a>
                              </xsl:when>
                              <xsl:otherwise>
                                <span class="note">No resources found.</span>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                        </tr>
                      </xsl:for-each>
                    </xsl:when>
                    <xsl:otherwise>
                      <tr>
                        <td colspan="5">
                          <strong>No available lanes for this flow cell.</strong>
                        </td>
                      </tr>
                    </xsl:otherwise>
                  </xsl:choose>
                </tbody>
              </table>
              <xsl:if test="/flow-cell/illumina-lane-index">
                <h2>lane index report</h2>
                <div class="legend_block">
                  <div class="label">
                    Legend:
                  </div>
                  <div class="legend">
                    <script type="text/javascript+protovis">
                      var w=700,
                      h=15,
                      c = pv.Colors.category10();

                      var legend = new pv.Panel()
                      .width(w)
                      .height(h);


                      legend.add(pv.Dot)
                      .data(window.indexData[0].sequence)
                      .top(8)
                      .left(function() 15 + this.index * 65)
                      .shape("square")
                      .size(16)
                      .strokeStyle(null)
                      .fillStyle(function(d) c(window.indexData[0].sequence[this.index]))
                      .anchor("right").add(pv.Label);

                      legend.render();
                    </script>
                  </div>
                </div>
                <div class="charts_block">
                  <xsl:for-each select="//report">
                    <xsl:apply-templates select="lane" />
                  </xsl:for-each>
                </div>

              </xsl:if>
            </div>
          </div>
        </div>
      </body>
    </html>

  </xsl:template>
  <xsl:template match="lane">
    <div class="lane_chart">
      <h3>Lane <xsl:value-of select="@number"/></h3>

      <!--
          <script type="text/javascript+protovis">
          new pv.Panel()
          .width(150)
          .height(150)
          .add(pv.Bar)
          .data(window.indexData.lane<xsl:value-of select="@number"/>.percent)
          .bottom(0)
          .width(20)
          .height(function(d) d * 80)
          .left(function() this.index * 25)
          .root.render();
          </script>
      -->
      <script type="text/javascript+protovis">
        var y_bars = window.indexData[0].sequence.length,
        w = 200,
        h = y_bars * 15,
        x = pv.Scale.linear(window.chartMin, window.chartMax).range(0, w),
        y = pv.Scale.ordinal(pv.range(y_bars)).splitBanded(0, h, 4/5),
        c = pv.Colors.category10();

        var vis =new pv.Panel()
        .width(w)
        .height(h)
        .bottom(20)
        .left(1)
        .right(10)
        .top(5);

        var bar = vis.add(pv.Bar)
        .data(window.indexData[<xsl:value-of select="@number  - 1"/>].percent)
        .top(function() y(this.index))
        .height(y.range().band)
        .left(0)
        .width(x)
        .fillStyle(function(d) c(window.indexData[<xsl:value-of select="@number  - 1"/>].sequence[this.index]));

        bar.anchor("right").add(pv.Label)
        .textStyle("white")
        .text(function(d) d.toFixed(2));

        // bar.anchor("left").add(pv.Label)
        //     .textMargin(5)
        //     .textAlign("right")
        //     .text(function() window.indexData[<xsl:value-of select="@number  - 1"/>].sequence[this.index]);

        vis.add(pv.Rule)
        .data(x.ticks())
        .left(function(d) Math.round(x(d)) - .5)
        .strokeStyle(function(d) d ? "rgba(255,255,255,.3)" : "#999")
        .add(pv.Rule)
        .bottom(0)
        .height(5)
        .strokeStyle("#999")
        .anchor("bottom").add(pv.Label)
        .text(function(d) d.toFixed(1))
        .textStyle("#999");

        vis.render();
      </script>
      <p class="axis_label_x">Index %</p>
    </div>

  </xsl:template>
</xsl:stylesheet>
