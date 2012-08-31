<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
  <xsl:template match="/">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
  <head>
    <title><xsl:value-of select="//report-meta/name"/></title> 
    <link href="layout.css" rel="stylesheet" type="text/css"></link>
    <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png"/>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen"/>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/css/tablesorter.css" type="text/css" media="screen" />
    <script src="https://imp.gsc.wustl.edu/res/old/report_resources/db_snp_concordance/js/jquery.js" type="text/javascript"></script>
    <script src="https://imp.gsc.wustl.edu/res/old/report_resources/db_snp_concordance/js/jquery.flot.js" type="text/javascript"></script>

    <!-- progression graph -->
    <script  type="text/javascript">
      $(function() {
        var datasets = {
          "lane error": {
            yaxis: 3,
            label: "Lane Error",          
            bars: { show: true },
            points: { show: false },
            lines: { show: false },
            data: [<xsl:for-each select="//progression/lane">[
                    <xsl:value-of select="lane-name"/>,
                    <xsl:value-of select="lane-error-average"/>,
                    <xsl:value-of select="lane-error-standard-deviation"/>],
                   </xsl:for-each>]
          },          
          "coverages": {
            yaxis: 1,
            label: "Total Coverage",
            data: [<xsl:for-each select="//progression/lane">[
                    <xsl:value-of select="lane-name"/>,
                    <xsl:value-of select="percent-coverage"/>,
                    <xsl:value-of select="bases-covered"/>,
                    <xsl:value-of select="reference-bases"/>],
                   </xsl:for-each>] 
          },
          "gains": {
            yaxis: 2,
            label: "Coverage Gains",
            data: [<xsl:for-each select="//progression/lane">[
                    <xsl:value-of select="lane-name"/>,
                    <xsl:value-of select="percent-gain"/>],
                   </xsl:for-each>]
          }
        };
        
      <xsl:text disable-output-escaping="yes">
        //<![CDATA[
        // hard-code color indices to prevent them from shifting as
        // plots are turned on/off
        var i = 0;
        $.each(datasets, function(key, val) {
            val.color = i;
            ++i;
        });
        
        // insert checkboxes 
        var choiceContainer = $("#progression_choices");
        $.each(datasets, function(key, val) {
            choiceContainer.append(
                '<div class="plot_checkbox"><input type="checkbox" name="' + key + '" checked="checked" id="' + key + '"/><label for="'+ key + '">' + val.label + '</label></div>'
            );
        });
        
	  choiceContainer.find("input").click(plotAccordingToChoices);
	  
        function showTooltip(x, y, contents) {
          $('<div id="tooltip">' + contents + '</div>').css( {
	          position: 'absolute',
	          display: 'none',
	          top: y + 5,
	          left: x + 5,
	          border: '1px solid #fdd',
	          padding: '2px',
	          'background-color': '#EEE',
	          opacity: 0.80
          }).appendTo("body").fadeIn(200);
        }

        var previousPoint = null;
        $("#progression").bind("plothover", function (event, pos, item) {
            $("#x").text(pos.x.toFixed(2));
            $("#y").text(pos.y.toFixed(2));
            
            if (item) {
                if (previousPoint != item.datapoint) {
                    previousPoint = item.datapoint;
                    
                    $("#tooltip").remove();
                    var x = item.datapoint[0].toFixed(2),
                    y = item.datapoint[1].toFixed(2);
                    
                    switch (item.series.label) {                    
                      case "Lane Error":
                        var avg_lane_error = item.datapoint[1]
                        var std_deviation = item.datapoint[2];
                        showTooltip(item.pageX, item.pageY, "<strong>avg. lane error:</strong> " + avg_lane_error + "<br/><strong>std. deviation:</strong> " + std_deviation);
                      break;
                        
                      case "Total Coverage":
                        var bases_covered = item.datapoint[2]
                        var reference_bases = item.datapoint[3];
                        showTooltip(item.pageX, item.pageY, "<strong>bases covered:</strong> " + addCommas(bases_covered) + "<br/><strong>reference bases:</strong> " + addCommas(reference_bases));
                      break;
                      
                      case "Coverage Gains":
                        var coverage_gain = item.datapoint[1]
                        showTooltip(item.pageX, item.pageY, "<strong>coverage gain:</strong> " + coverage_gain + "%");
                      break;
                    }
                }
            }
            else {
                $("#tooltip").remove();
                previousPoint = null;
            }
        });

        function addPercentage(number) {
          return number + "%";
        }

        function addCommas(number) {
              number = '' + number;
              if (number.length > 3) {
                  var mod = number.length % 3;
                  var output = (mod > 0 ? (number.substring(0,mod)) : '');
                  for (i=0 ; i < Math.floor(number.length / 3); i++) {
                      if ((mod == 0) && (i == 0)) {
                          output += number.substring(mod+ 3 * i, mod + 3 * i + 3);
                      } else {
                          output+= ',' + number.substring(mod + 3 * i, mod + 3 * i + 3);
                      }
                  }
                  return (output);
              } else { return number; };
        }

        function plotAccordingToChoices() {
            var data = [];
            
            choiceContainer.find("input:checked").each(function () {
                var key = $(this).attr("name");
                if (key && datasets[key])
                    data.push(datasets[key]);
            });
            
            if (data.length > 0)
                $.plot($("#progression"),data,
                    {
                      legend: { show: true, position: "ne", margin: 40 },
                      lines: { show: true },
                      points: { show: true },
                      xaxis: { min: 1, tickDecimals: 0 },
                      yaxis: { tickSize: 1, tickFormatter: function (v, axis) { return addPercentage(v); } },
                      y2axis: { tickSize: 1, tickFormatter: function (v, axis) { return addPercentage(v); } },
                      y3axis: { min: 0, max: 1 },
                      grid: { hoverable: true, clickable: true }
                });
        }
        
        plotAccordingToChoices();
      });
      /* ]]> */
      </xsl:text>
    </script>
    
    <!-- depth by relative position graph -->
    <script  type="text/javascript">
      $(function() {
        var datasets = {
          <xsl:for-each select="//relative-coverage/relative-coverage-bin">
          "<xsl:value-of select="@size"/>": {
            label: "<xsl:value-of select="@size"/>",
            data: [<xsl:for-each select="reading">[<xsl:value-of select="relative-position"/>,<xsl:value-of select="depth"/>],</xsl:for-each>] 
          },
          </xsl:for-each>
        };
        
      <xsl:text disable-output-escaping="yes">
        //<![CDATA[
        // hard-code color indices to prevent them from shifting as
        // plots are turned on/off
        var i = 0;
        $.each(datasets, function(key, val) {
            val.color = i;
            ++i;
        });
        
        // insert checkboxes 
        var choiceContainer = $("#depth_by_relative_position_choices");
        $.each(datasets, function(key, val) {
            choiceContainer.append(
                '<div class="plot_checkbox"><input type="checkbox" name="' + key + '" checked="checked" id="' + key + '"/><label for="'+ key + '">' + val.label + '</label></div>'
            );
        });
        
	  choiceContainer.find("input").click(plotAccordingToChoices);
	  
        function showTooltip(x, y, contents) {
          $('<div id="tooltip">' + contents + '</div>').css( {
	          position: 'absolute',
	          display: 'none',
	          top: y + 5,
	          left: x + 5,
	          border: '1px solid #fdd',
	          padding: '2px',
	          'background-color': '#EEE',
	          opacity: 0.80
          }).appendTo("body").fadeIn(200);
        }

        var previousPoint = null;
        $("#depth_by_relative_position").bind("plothover", function (event, pos, item) {
            $("#x").text(pos.x.toFixed(2));
            $("#y").text(pos.y.toFixed(2));
            
            if (item) {
                if (previousPoint != item.datapoint) {
                    previousPoint = item.datapoint;
                    
                    $("#tooltip").remove();
                    var x = item.datapoint[0].toFixed(2),
                    y = item.datapoint[1].toFixed(2);
                    var depth = item.datapoint[1]
                    
                    showTooltip(item.pageX, item.pageY, "<strong>depth:</strong> " + addCommas(depth)); 
                }
            }
            else {
                $("#tooltip").remove();
                previousPoint = null;
            }
        });

        function addPercentage(number) {
          return number + "%";
        }

        function addCommas(number) {
              number = '' + number;
              if (number.length > 3) {
                  var mod = number.length % 3;
                  var output = (mod > 0 ? (number.substring(0,mod)) : '');
                  for (i=0 ; i < Math.floor(number.length / 3); i++) {
                      if ((mod == 0) && (i == 0)) {
                          output += number.substring(mod+ 3 * i, mod + 3 * i + 3);
                      } else {
                          output+= ',' + number.substring(mod + 3 * i, mod + 3 * i + 3);
                      }
                  }
                  return (output);
              } else { return number; };
        }


        function plotAccordingToChoices() {
            var data = [];
            
            choiceContainer.find("input:checked").each(function () {
                var key = $(this).attr("name");
                if (key && datasets[key])
                    data.push(datasets[key]);
            });
            
            if (data.length > 0)
                $.plot($("#depth_by_relative_position"),data,
                    {
                      legend: { show: true, position: "ne", margin: 20},
                      lines: { show: true },
                      points: { show: false },
                      xaxis: {  },
                      yaxis: { tickFormatter: function (v, axis) { return addCommas(v); } },
                      grid: { hoverable: true, clickable: true }
                });
        }
        
        plotAccordingToChoices();
      });
      /* ]]> */
      </xsl:text>
    </script>
    
    <style type="text/css" media="screen">
      div.graph {
          width: 100%;
          height: 450px;
      }

      div.content_padding {
          padding: 0 10px 20px 10px;
      }

      div.choices {
          width: 98%;
          padding: 1%;
          margin-top: 15px;
          background: #EFEFED;
          border: 1px solid #CCC;
          float: left;
          font-size: 85%;
      }
      
      div#choices p {
          margin: 0 0 5px 0;
          padding: 0;
      }
      
      div.plot_checkbox {
          display: inline-block;
          float: left;
          padding: 0 10px;
      }

      span.small {
        font-size: 85%;
        color: #999;
      }

    </style>
  </head>
  <body>
    <div class="container"><div class="background">
      <h1 class="page_title"><xsl:value-of select="//report-meta/name"/></h1>
      <div class="page_padding">
       <table cellpadding="0" cellspacing="0" border="0" class="info_table_group">
          <tr>
            <td>
              <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                <tr><td class="label">Model ID:</td><td class="value"><a><xsl:attribute name="href"><xsl:text>https://imp.gsc.wustl.edu/view/Genome/Model/status.html?id=</xsl:text><xsl:value-of select="//model-info/id"/></xsl:attribute><xsl:value-of select="//model-info/id"/></a></td></tr> 
                <tr><td class="label">Model Name:</td><td class="value"><xsl:value-of select="//model-info/name"/></td></tr>
                <tr><td class="label">Model Type:</td><td class="value"><xsl:value-of select="//model-info/type-name"/></td></tr>
                <tr><td class="label">Subject Name:</td><td class="value"><xsl:value-of select="//model-info/subject-name"/></td></tr>
                <tr><td class="label">Subject Type:</td><td class="value"><xsl:value-of select="//model-info/subject-type"/></td></tr>
                <tr><td class="label">Read Aligner:</td><td class="value"><xsl:value-of select="//model-info/read-aligner-name"/> v<xsl:value-of select="//model-info/read-aligner-version"/></td></tr>
              </table>
            </td>
            <td>
              <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                <tr><td class="label">Build ID:</td><td class="value"><a><xsl:attribute name="href"><xsl:text>https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=</xsl:text><xsl:value-of select="//model-info/build-id"/></xsl:attribute><xsl:value-of select="//model-info/build-id"/></a></td></tr>
                <tr><td class="label">Data Directory:</td><td class="value"><a><xsl:attribute name="href"><xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="//model-info/data-directory"/></xsl:attribute>view build directory</a>
</td></tr>                
                <tr><td class="label">Processing Profile:</td><td class="value"><xsl:value-of select="//model-info/processing-profile-name"/></td></tr> 
                <tr><td class="label">Sequencing Platform:</td><td class="value"><xsl:value-of select="//model-info/sequencing-platform"/></td></tr>
                <tr><td class="label">DNA Type:</td><td class="value"><xsl:value-of select="//model-info/dna-type"/></td></tr>
                <tr><td class="label">Reference Sequence:</td><td class="value"><xsl:value-of select="//model-info/reference-sequence-name"/></td></tr>
              </table>
            </td>
          </tr>
        </table>
        <h2 class="report_section" style="margin-bottom: 0;">cDNA Aligment Summary</h2>
        <table id="cdna_alignment_summary" class="list" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
          <colgroup>
            <col width="40%" />
            <col />
            <col />
            <col />
          </colgroup>
          <thead>
            <tr>
              <th>category</th>
              <th class="last">items</th>
              <th class="last">total</th>
              <th class="last">%</th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="//summary/category">
              <tr>
                <td><xsl:value-of select="category-name"/></td>
              	<td class="last"><xsl:variable name="category-items" select="category-items"/><xsl:value-of select="format-number($category-items, '#,##0')"/></td>
              	<td class="last"><xsl:variable name="total" select="total"/><xsl:value-of select="format-number($total, '#,##0')"/></td>
              	<td class="last"><xsl:value-of select="percent-category"/>%</td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>

        <h2 class="report_section" style="margin-bottom: 0;">Alignment Breakdown</h2>
        <table id="alignemnt_breakdown" class="list" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
          <colgroup>
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />            
          </colgroup>
          <thead>
            <tr>
              <th>category</th>
              <th class="last">reads</th>
              <th class="last">total reads</th>
              <th class="last">read %</th>
              <th class="last">poly AT</th>
              <th class="last">poly AT%</th>
              <th class="last">poly A</th>
              <th class="last">poly A%</th>
              <th class="last">poly T</th>
              <th class="last">poly T%</th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="//breakdown/category">
              <tr>
                <td><xsl:value-of select="category-name"/></td>

              	<td class="last"><xsl:variable name="reads" select="reads"/><xsl:value-of select="format-number($reads, '#,##0')"/></td>
                <td class="last"><xsl:variable name="percent-total-reads" select="percent-total-reads"/><xsl:value-of select="format-number($percent-total-reads, '#,##0')"/></td>
                <td class="last"><xsl:value-of select="percent-category"/>%</td>
                <td class="last"><xsl:variable name="poly-at" select="poly-at"/><xsl:value-of select="format-number($poly-at, '#,##0')"/></td>
                <td class="last"><xsl:value-of select="percent-poly-at"/>%</td>
                <td class="last"><xsl:variable name="poly-a" select="poly-a"/><xsl:value-of select="format-number($poly-a, '#,##0')"/></td>
                <td class="last"><xsl:value-of select="percent-poly-a"/>%</td>
                <td class="last"><xsl:variable name="poly-t" select="poly-t"/><xsl:value-of select="format-number($poly-t, '#,##0')"/></td>
                <td class="last"><xsl:value-of select="percent-poly-t"/>%</td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>

        <h2 class="report_section">Coverage Progression by Lane</h2>
          <table width="100%" cellpadding="5" cellspacing="0">
            <tr>
              <td valign="middle"><img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/axis_label_y_total_coverage.png" width="25" height="141" alt="Total Coverage (%)" style="margin-right: 10px;"/></td>
              <td align="center" valign="middle" width="100%">
                <div id="progression" class="graph"></div>
              </td>
              <td valign="middle"><img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/axis_label_y2_coverage_gain.png" width="25" height="140" alt="Coverage Gain (%)" style="margin-right: 10px;"/></td>
            </tr>            
            <tr>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
              <td valign="middle" align="center"><img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/axis_label_x_lane.png" width="47" height="20" alt="Gigabases" style="margin-top: 10px;"/></td>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>
            <tr>
              <td colspan="3">
                <div id="progression_choices" class="choices">
                  <p>
                    <strong>Show:</strong>
                  </p>
                </div>  
              </td>
            </tr>
          </table>

        <h2 class="report_section">Depth by Relative Position</h2>
          <table width="100%" cellpadding="5" cellspacing="0">
            <tr>
              <td valign="middle"><img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/axis_label_y_depth.png" width="25" height="55" alt="Depth" style="margin-right: 10px;"/></td>
              <td align="center" valign="middle" width="100%">
                <div id="depth_by_relative_position" class="graph"></div>
              </td>
              <td valign="middle"><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>            
            <tr>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
              <td valign="middle" align="center"><img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/axis_label_x_relative_position.png" width="134" height="17" alt="Relative Position" style="margin-top: 10px;"/></td>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>
            <tr>
              <td colspan="3">
                <div id="depth_by_relative_position_choices" class="choices">
                  <p>
                    <strong>Show:</strong>
                  </p>
                </div>  
              </td>
            </tr>
          </table>

       <h2 class="report_section" style="margin-bottom: 0;">Breadth and Depth of Coverage Summary</h2>
        <table id="breadth_depth_coverage_summary" class="list" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
          <colgroup>
            <col />
            <col />
            <col />
            <col />
            <col />
            <col />
          </colgroup>
          <thead>
            <tr>
              <th>reference size</th>
              <th>minimum depth</th>
              <th class="last">minimum coverage %</th>
              <th class="last">covered references</th>
              <th class="last">total references</th>
              <th class="last">covered %</th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="//coverage-bins/coverage-bin">
              <tr>
                <td><xsl:value-of select="reference-size"/><br/><span class="small">(<xsl:value-of select="minimum-basepair"/> <xsl:value-of select="maximum-basepair"/>)</span></td>
              	<td><xsl:value-of select="minimum-depth"/></td>
                <td class="last"><xsl:value-of select="minimum-percent-coverage"/>%</td>
                <td class="last"><xsl:variable name="covered-references" select="covered-references"/><xsl:value-of select="format-number($covered-references, '#,##0')"/></td>
              	<td class="last"><xsl:variable name="total-references" select="total-references"/><xsl:value-of select="format-number($total-references, '#,##0')"/></td>
                <td class="last"><xsl:value-of select="percent-covered"/>%</td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>

      </div>
    </div>
  </div>
</body>

</html>

  </xsl:template>
</xsl:stylesheet>
