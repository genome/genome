<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:template match="/">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
    <title><xsl:value-of select="//report-meta/name"/></title>
    <link href="layout.css" rel="stylesheet" type="text/css"></link>
    <link rel="shortcut icon" href="/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png"/>
    <link rel="stylesheet" href="/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen"/>
    <link rel="stylesheet" href="/res/old/report_resources/apipe_dashboard/css/tablesorter.css" type="text/css" media="screen" />
    <script src="/res/old/report_resources/db_snp_concordance/js/jquery.js" type="text/javascript"></script>
    <script src="/res/old/report_resources/db_snp_concordance/js/jquery.flot.js" type="text/javascript"></script>
    <script type="text/javascript" src="/res/old/report_resources/jquery/jquery.tablesorter.min.js"></script>
    <script type="text/javascript">
      $(document).ready(function() {
          $("#model_list").tablesorter({
          // sort on first column, ascending
          sortList: [[0,0]]
          });
      });
    </script>

    <script  type="text/javascript">
      $(function() {
      var datasets = {
      <xsl:for-each select="//dataset/model">
        "<xsl:value-of select="@name"/>": {
            label: "<xsl:value-of select="@name"/>",
            hoverable: true,
            clickable: true,
            shadowSize: 0,
            data:[<xsl:for-each select="build">[ <xsl:value-of select="gigabases"/>, <xsl:value-of select="goldsnp-concordance-filtered"/>, <xsl:value-of select="@id"/>, <xsl:value-of select="lanes"/> ], </xsl:for-each>]
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
                      var build_id = item.datapoint[2]
                      var label = item.datapoint[3];

                      // showTooltip(item.pageX, item.pageY, item.series.label + " at " + x + " = " + y);
                      showTooltip(item.pageX, item.pageY, "<strong>build id:</strong> " + build_id + "<br/><strong>lanes:</strong> " + label + "<br/><strong>goldSNP concordance:</strong> " + y + "%");
                  }
              }
              else {
                  $("#tooltip").remove();
                  previousPoint = null;
              }
          });

          $("#placeholder").bind("plotclick", function (event, pos, item) {
              if (item) {
                window.open('https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=' + item.datapoint[2]);
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
                      legend: { show: true, position: "se", margin: 20},
                      lines: { show: true },
                      points: { show: true },
                      xaxis: { min: 0, tickFormatter: function (v, axis) { return addCommas(v); } },
                      yaxis: { ticks: 11, min: 0, max: 100},
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
    <script type="text/javascript" src="/res/old/report_resources/jquery/boxy/src/javascripts/jquery.boxy.js"></script>
    <link rel="stylesheet" href="/res/old/report_resources/jquery/boxy/src/stylesheets/boxy.css" type="text/css" />
    <script type="text/javascript">
      <![CDATA[
               function event_popup(eventObject) {
                   // assemble event info into a table
                   var popup_content = '<table class="boxy_info" cellpadding="0" cellspacing="0" border="0"><tbody>';
                   for (prop in eventObject) {
                       if (prop != 'popup_title') {
                           popup_content += '<tr><td class="label">' + prop.replace(/_/g," ") + ':</td><td class="value">' + eventObject[prop] + '</td></tr>';
                       }
                   }

                   popup_content += '</tbody></table>';
                   // create popup
                   var popup = new Boxy(popup_content, {title:eventObject.popup_title, fixed:false});
                   popup.center();
               }
      ]]>
    </script>

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
		            <h1>Analysis Reports v0.2</h1>
		          </td>
		        </tr>
		      </table>
		    </div>
		    <div class="page_padding">


<h2 class="page_title"><xsl:value-of select="//report-meta/name"/></h2>
      <div class="page_padding">

        <h2 class="report_section" style="margin-bottom: 0;">Models Compared</h2>
        <table id="model_list" class="list tablesorter" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
          <colgroup>
            <col width="40%" />
            <col />
            <col />
            <col />
          </colgroup>
          <thead>
            <tr>
              <th>model name</th>
              <th>model ID</th>
              <th>processing profile</th>
              <th class="last">username</th>
            </tr>
          </thead>
          <tbody>
                <xsl:for-each select="//dataset/model">
                  <tr onmouseover="this.className = 'hover'" onmouseout="this.className=''">
                    <xsl:attribute name="onclick">
                      <xsl:text>javascript:document.location.href='</xsl:text>
                      <xsl:call-template name="object_link_href"/>
                      <xsl:text>'</xsl:text>
                    </xsl:attribute>
                    <td><xsl:value-of select="@name"/></td>
                    <td><xsl:value-of select="@id"/></td>
                    <td><xsl:value-of select="@processing-profile"/></td>
                    <td class="last"><xsl:value-of select="@username"/></td>
                  </tr>
                </xsl:for-each>
          </tbody>
        </table>

        <h2 class="report_section">Comparison Graph</h2>
          <table width="100%" cellpadding="5" cellspacing="0">
            <tr>
              <td valign="middle"><img src="/res/old/report_resources/apipe_dashboard/images/axis_label_v_snp_concordance.png" width="21" height="204" alt="SNP Concordance %" style="margin-right: 10px;"/></td>
              <td align="center" valign="middle" width="100%">
                <div id="placeholder" class="graph_placeholder"></div>
              </td>
              <td valign="middle"><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
            </tr>
            <tr>
              <td><xsl:text disable-output-escaping="yes">&amp;nbsp;</xsl:text></td>
              <td valign="middle" align="center"><img src="/res/old/report_resources/apipe_dashboard/images/axis_label_h_gigabases.png" width="104" height="26" alt="Gigabases" style="margin-top: 10px;"/></td>
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
        </div>
        <div id="ajax_status"/>
      </body>
    </html>

  </xsl:template>

</xsl:stylesheet>
