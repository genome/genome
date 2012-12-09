<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
  <xsl:template match="/">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
    <title><xsl:value-of select="//report-meta/description"/></title> 
    <link href="layout.css" rel="stylesheet" type="text/css"></link>
    <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/static/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png"/>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/static/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen"/>
    <script src="https://imp.gsc.wustl.edu/static/report_resources/db_snp_concordance/js/jquery.js" type="text/javascript"></script>
    <script src="https://imp.gsc.wustl.edu/static/report_resources/db_snp_concordance/js/jquery.flot.js" type="text/javascript"></script>
    <script  type="text/javascript">
      $(function() {
      var datasets = {
      <xsl:for-each select="//qualities">
        "<xsl:value-of select="@label"/>": {
            label: "<xsl:value-of select="@label"/>",
            hoverable: true,
            clickable: true,
            shadowSize: 0,
            data:[<xsl:for-each select="quality">[ <xsl:value-of select="position"/>, <xsl:value-of select="value"/> ], </xsl:for-each>]
        },          
      </xsl:for-each>
      };

      <xsl:text disable-output-escaping="yes">
              //<![CDATA[
              // hard-code color indices to prevent them from shifting as
              // countries are turned on/off
              var i = 0;
          $.each(datasets, function(key, val) {
              val.color = i;
              ++i;
          });
          
          // insert checkboxes 
          var choiceContainer = $("#choices");
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
		          'background-color': '#fee',
		          opacity: 0.80
	          }).appendTo("body").fadeIn(200);
          }

          var previousPoint = null;
          $("#placeholder").bind("plothover", function (event, pos, item) {
              $("#x").text(pos.x.toFixed(2));
              $("#y").text(pos.y.toFixed(2));
              
              if (item) {
                  if (previousPoint != item.datapoint) {
                      previousPoint = item.datapoint;
                      
                      $("#tooltip").remove();
                      var x = item.datapoint[0].toFixed(2),
                      y = item.datapoint[1].toFixed(2);
                      
                      showTooltip(item.pageX, item.pageY, item.series.label + " at " + x + " = " + y);
                  }
              }
              else {
                  $("#tooltip").remove();
                  previousPoint = null;
              }
          });

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
                  $.plot($("#placeholder"), data, {
                      legend: { show: true, position: "ne", margin: 20},
                      lines: { show: true },
                      points: { show: false },
                      xaxis: { tickFormatter: function (v, axis) { return addCommas(v); } },
                      yaxis: { ticks: 11, min: 0, max: 100 },
                      grid: { hoverable: true, clickable: true }
                  });
          }
          
          plotAccordingToChoices();
      });
      /* ]]> */
      </xsl:text>
    </script>
    <style type="text/css" media="screen">
      div#placeholder {
          width: 100%;
          height: 450px;
      }

      div.content_padding {
          padding: 0 10px 20px 10px;
      }

      div#choices {
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

    </style>

  </head>
  <body>
    <div class="container"><div class="background">
      <h1 class="page_title"><xsl:value-of select="//report-meta/description"/></h1>
      <div class="page_padding">
        
        <table cellspacing="0" cellpadding="0" border="0" class="info_table_group">
          <tr>
            <td>
              <table cellspacing="0" cellpadding="0" border="0" class="info_table">
                <tr>
                  <td class="label">build id</td>
                  <td class="value"><xsl:value-of select="//report-meta/generator-params/build-id"/></td>
                </tr>
				<tr>
				  <td class="label">data directory</td>
				  <td class="value"><a><xsl:attribute name="href">https://<xsl:value-of select="//model-info/data-directory"/></xsl:attribute><xsl:text>View Directory</xsl:text></a></td>
				</tr>
              </table>
            </td>
            <td>
              <table cellspacing="0" cellpadding="0" border="0" class="info_table">
                <colgroup>
                  <col width="50%"/>
                  <col width="50%"/>
                </colgroup>
                <tr>
                  <td class="label">report generated</td>
                  <td class="value"><xsl:value-of select="//report-meta/date"/></td>
                </tr>
                <tr>
                  <td class="label">generator</td>
                  <td class="value"><xsl:value-of select="//report-meta/generator"/></td>
                </tr>
              </table>
            </td>
          </tr>
        </table>

        <h2 class="report_section">Model Information</h2>
        <table cellspacing="0" cellpadding="0" border="0" class="info_table_group">
          <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">model name</td>
				  <td class="value"><xsl:value-of select="//model-info/name"/></td>
				</tr>
				<tr>
				  <td class="label">subject name</td>
				  <td class="value"><xsl:value-of select="//model-info/subject-name"/></td>
				</tr>
				<tr>
				  <td class="label">subject type</td>
				  <td class="value"><xsl:value-of select="//model-info/subject-type"/></td>
				</tr>
				<tr>
				  <td class="label">purpose</td>
				  <td class="value"><xsl:value-of select="//model-info/purpose"/></td>
				</tr>
				<tr>
				  <td class="label">assembler</td>
				  <td class="value"><xsl:value-of select="//model-info/assembler"/></td>
				</tr>
				<tr>
				  <td class="label">sequencing platform</td>
				  <td class="value"><xsl:value-of select="//model-info/sequencing-platform"/></td>
				</tr>
				<tr>
				  <td class="label">sequencing center</td>
				  <td class="value"><xsl:value-of select="//model-info/sequencing-center"/></td>
				</tr>
			  </table>
			</td>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">region of interest</td>
				  <td class="value"><xsl:value-of select="//model-info/region-of-interest"/></td>
				</tr>
				<tr>
				  <td class="label">assembly size</td>
				  <td class="value"><xsl:value-of select="//model-info/assembly-size"/></td>
				</tr>
				<tr>
				  <td class="label">reverse seq primer</td>
				  <td class="value"><xsl:value-of select="//model-info/primer-seq-reverse"/></td>
				</tr>
				<tr>
				  <td class="label">reverse amp primer</td>
				  <td class="value"><xsl:value-of select="//model-info/primer-amp-reverse"/></td>
				</tr>
				<tr>
				  <td class="label">forward seq primer</td>
				  <td class="value"><xsl:value-of select="//model-info/primer-seq-forward"/></td>
				</tr>
				<tr>
				  <td class="label">forward amp primer</td>
				  <td class="value"><xsl:value-of select="//model-info/primer-amp-forward"/></td>
				</tr>
			  </table>
			</td>
          </tr>
        </table>

        <h2 class="report_section">Assembly Statistics</h2>
        <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
          <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">attempted</td>
				  <td class="value"><xsl:value-of select="//stat/attempted"/></td>
				</tr>
				<tr>
				  <td class="label">assembled</td>
				  <td class="value"><xsl:value-of select="//stat/assembled"/></td>
				</tr>
				<tr>
				  <td class="label">assembly success</td>
				  <td class="value"><xsl:value-of select="//stat/assembly-success"/></td>
				</tr>

			  </table>
			</td>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">average length</td>
				  <td class="value"><xsl:value-of select="//stat/length-average"/></td>
				</tr>
				<tr>
				  <td class="label">max. length</td>
				  <td class="value"><xsl:value-of select="//stat/length-maximum"/></td>
				</tr>
				<tr>
				  <td class="label">median length</td>
				  <td class="value"><xsl:value-of select="//stat/length-median"/></td>
				</tr>
				<tr>
				  <td class="label">minimum length</td>
				  <td class="value"><xsl:value-of select="//stat/length-minimum"/></td>
				</tr>
			  </table>
			</td>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">average base quality</td>
				  <td class="value"><xsl:value-of select="//stat/quality-base-average"/></td>
				</tr>
				<tr>
				  <td class="label">less than 20 bases<br/> per assembly quality</td>
				  <td class="value"><xsl:value-of select="//stat/quality-less-than-20-bases-per-assembly"/></td>
				</tr>
			  </table>
			</td>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
				  <col width="50%"/>
				  <col width="50%"/>
				</colgroup>
				<tr>
				  <td class="label">assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled"/></td>
				</tr>
				<tr>
				  <td class="label">average assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled-average"/></td>
				</tr>
				<tr>
				  <td class="label">maximum assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled-maximum"/></td>
				</tr>
				<tr>
				  <td class="label">median assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled-median"/></td>
				</tr>
				<tr>
				  <td class="label">miniumum assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled-minimum"/></td>
				</tr>
				<tr>
				  <td class="label">successful assembled reads</td>
				  <td class="value"><xsl:value-of select="//stat/reads-assembled-success"/></td>
				</tr>
			  </table>
			</td>
          </tr>
        </table>
        <h2 class="report_section">Quality Histogram</h2>
          <table width="100%" cellpadding="5" cellspacing="0">
            <tr>
              <td valign="middle"><img src="/static/report_resources/apipe_dashboard/images/axis_label_v_quality.png" width="25" height="68"/></td>
              <td align="center" valign="middle" width="100%">
                <div id="placeholder" class="graph_placeholder"></div>
              </td>
              <td valign="middle"><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>            
            <tr>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
              <td valign="middle" align="center"><img src="/static/report_resources/apipe_dashboard/images/axis_label_h_length.png" width="68" height="26"/></td>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>
            <tr>
              <td colspan="3">
                <div id="choices">
                  <p>
                    <strong>Show:</strong>
                  </p>
                </div>  
              </td>
            </tr>
          </table>

      </div>
    </div>
  </div>
</body>

</html>

  </xsl:template>
</xsl:stylesheet>
