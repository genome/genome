<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_instrumentdata_solexa_quality" match="/report">
    <xsl:comment>file: html/quality/genome_instrumentdata_solexa.xsl name:genome_instrumentdata_solexa_quality</xsl:comment>
    <script type="text/javascript" src="/res/js/pkg/protovis-r3.1.js"></script>
    <script type="text/javascript" src="/res/js/app/quality/genome_instrumentdata_solexa_bar.js"></script>
    <script type="text/javascript" src="/res/js/app/quality/genome_instrumentdata_solexa_candlestick.js"></script>


    <script type="text/javascript">
      <xsl:for-each select="//quality-stats/read-set">

        var <xsl:value-of select="@read-set-name"/>_data = [
        <xsl:for-each select="cycle">
          {
          column: <xsl:value-of select="column"/>,
          count: <xsl:value-of select="count"/>,
          quality_min: <xsl:value-of select="min"/>,
          quality_max: <xsl:value-of select="max"/>,
          quality_sum: <xsl:value-of select="sum"/>,
          quality_mean: <xsl:value-of select="mean"/>,
          quartile_q1: <xsl:value-of select="Q1"/>,
          quartile_med: <xsl:value-of select="med"/>,
          quartile_q3: <xsl:value-of select="Q3"/>,
          quartile_iqr: <xsl:value-of select="IQR"/>,
          whisker_left: <xsl:value-of select="lW"/>,
          whisker_right: <xsl:value-of select="rW"/>,
          count_a: <xsl:value-of select="A-Count"/>,
          count_c: <xsl:value-of select="C-Count"/>,
          count_g: <xsl:value-of select="G-Count"/>,
          count_t: <xsl:value-of select="T-Count"/>,
          count_n: <xsl:value-of select="N-Count"/>,
          read_set_name: "<xsl:value-of select="../@read-set-name"/>"
          },
        </xsl:for-each>
        ];
      </xsl:for-each>
    </script>

    <style type="text/css" media="screen">

      div.content_padding {
      padding: 0 10px 20px 10px;
      }

      div.graph_block {
      width: 460px;
      float: left;
      }

      table.key {
      width: auto;
      margin: 0;
      padding: 0;
      height: 17px;
      }

      table.key td {
      margin: 0;
      padding: 0;
      height: 17px;
      }

      table.key td.title{
      font-size: 10px;
      line-height: 10px;
      font-weight: bold;
      padding-right: 10px;
      }

      table.key td.graphic {
      width: 12px;
      height: 12px;
      padding 0;
      margin: 0;
      }

      table.key td.value {
      font-size: 10px;
      line-height: 10px;
      padding-left: 3px;
      padding-right: 10px;
      }

    </style>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Instrument Data Quality:'" />
      <xsl:with-param name="display_name" select="//instrument-data-info/id" />
      <xsl:with-param name="icon" select="'genome_instrumentdata_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div class="box rounded">
          <div style="width: 49%; margin-right; 2%; display: inline-block">
            <table border="0" cellpadding="0" cellspacing="0" class="name-value" style="margin:0;">
              <tr>
                <td class="name">instrument data ID:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/id"/></td>
              </tr>

              <tr>
                <td class="name">flow cell ID:</td>
                <td class="value"><a class="mini btn"><xsl:attribute name="href">/view/genome/instrument-data/flow-cell/status.html?id=<xsl:value-of select="//report/@flow_cell_id"/></xsl:attribute><xsl:value-of select="//report/@flow_cell_id"/></a> (lane <xsl:value-of select="//report/@lane"/>)</td>
              </tr>

              <tr>
                <td class="name">platform:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/sequencing-platform"/></td>
              </tr>

              <tr>
                <td class="name">run name:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/run-name"/></td>
              </tr>

              <tr>
                <td class="name">library name:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/library-name"/></td>
              </tr>


            </table>
          </div>

          <div style="width: 49%; display: inline-block;">
            <table border="0" cellpadding="0" cellspacing="0" class="name-value" style="margin:0;">
              <tr>
                <td class="name">sample name:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/sample-name"/></td>
              </tr>

              <tr>
                <td class="name">subset name:</td>
                <td class="value"><xsl:value-of select="//instrument-data-info/subset-name"/></td>
              </tr>

              <tr>
                <td class="name">analysis version:</td>
                <td class="value"><xsl:value-of select="//report/@analysis_software_version"/></td>
              </tr>

              <tr>
                <td class="name">read length:</td>
                <td class="value"><xsl:value-of select="//report/@read_length"/></td>
              </tr>

              <tr>
                <td class="name">clusters:</td>
                <td class="value"><xsl:value-of select="format-number(//report/@clusters, '###,###,###')"/></td>
              </tr>

            </table>

          </div>
        </div>


        <div class="span-24 last">
          <div class="box_header span-24 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-24 last">
              Quality Stats and Nucleotide Distribution <span style="font-size: 75%; color: #666">(hover over columns for details)</span>
            </h3></div>
          </div>
          <div class="box_content rounded-bottom span-24 last">
            <div style="background-color: #FFF; padding: 15px;float:left;border-bottom: 1px solid #acaca3;margin-bottom:10px;">
              <xsl:for-each select="//quality-stats/read-set">
                <div class="graph_block">
                  <h3><xsl:value-of select="@read-set-name"/> quality stats</h3>
                  <table cellpadding="0" cellspacing="0" border="0" class="key">
                    <tr>
                      <td class="title">Legend:</td>

                      <td class="graphic" style="background-color: #d2d0a5;"></td>
                      <td class="value">whiskers</td>

                      <td class="graphic" style="background-color: #60604b;"></td>
                      <td class="value">high/low quartiles</td>

                      <td class="graphic" style="background-color: #f3b028;"></td>
                      <td class="value">med. quartile</td>
                    </tr>
                  </table>
                  <table width="100%" cellpadding="0" cellspacing="0" border="0">
                    <tr>
                      <td valign="middle">
                        <img src="/res/img/legacy/axis_label_v_quality_sm.png" width="17" height="44" alt="Quality"/>
                      </td>
                      <td>
                        <script type="text/javascript">
                          render_candlestick_graph(<xsl:value-of select="@read-set-name"/>_data, 385, 400);
                        </script>
                      </td>
                    </tr>
                    <tr>
                      <td></td>
                      <td><div align="center"><img src="/res/img/legacy/axis_label_h_col_count.png" width="93" height="16" alt="Cycle/Column"/></div></td>
                    </tr>
                  </table>
                </div>

                <div class="graph_block" style="">
                  <h3><xsl:value-of select="@read-set-name"/> nucleotide distribution</h3>
                  <table cellpadding="0" cellspacing="0" border="0" class="key">
                    <tr>
                      <td class="title">Legend:</td>

                      <td class="graphic" style="background-color: #7d598c;"></td>
                      <td class="value">A</td>

                      <td class="graphic" style="background-color: #bb4250;"></td>
                      <td class="value">C</td>

                      <td class="graphic" style="background-color: #90c86f;"></td>
                      <td class="value">G</td>

                      <td class="graphic" style="background-color: #f6c460;"></td>
                      <td class="value">T</td>

                      <td class="graphic" style="background-color: #999;"></td>
                      <td class="value">N</td>

                    </tr>
                  </table>

                  <table width="100%" cellpadding="0" cellspacing="0" border="0">
                    <tr>
                      <td>
                        <script type="text/javascript">
                          render_bar_graph(<xsl:value-of select="@read-set-name"/>_data, 385, 400);
                        </script>
                      </td>
                      <td valign="middle">
                        <img src="/res/img/legacy/axis_label_v_read_count.png" width="13" height="71" alt="Read Count"/>
                      </td>
                    </tr>
                    <tr>
                      <td><div align="center"><img src="/res/img/legacy/axis_label_h_col_count.png" width="93" height="16" alt="Cycle/Column"/></div></td>
                      <td></td>
                    </tr>
                  </table>
                </div>

              </xsl:for-each>
            </div>
          </div>
        </div>
      </div>
    </div>

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>
</xsl:stylesheet>
