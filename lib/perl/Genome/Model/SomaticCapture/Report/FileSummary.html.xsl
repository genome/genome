<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
                xmlns:set="http://exslt.org/sets">

  <xsl:output method="html"/>
  <xsl:output encoding="utf-8"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:template match="/">

    <html>
      <head>
        <title>File Summary: <xsl:value-of select="//model-info/name"/>&#160;<xsl:value-of select="//model-info/id"/></title>
        <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png" />
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen" />
        <style type="text/css">
        .page_footer {
            background-color: #f2f3e0;
            border-top: 1px solid #8dc643;
            padding: 5px;
        } 
        </style>
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/jquery/jquery.js"></script>
        
        <!-- initialize data tables -->
        <!-- note: dataTables doesn't like to be applied to a table with no column headers (which will happen if we create a 'None found' table), so must be applied using $(document).ready on a per-table basis in the body of the page. -->
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.js"></script>
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.plugin.formatted-num.js"></script>
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/jquery/dataTables-1.5/media/css/gc_table.css" type="text/css" media="screen"></link>
      </head>

      <body>
        <div class="container">
          <div class="background">
            <div class="page_header">
              <table cellpadding="0" cellspacing="0" border="0">
                <tr>
                  <td>
                    <img src="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/gc_header_logo2.png" width="44" height="45" align="absmiddle" />
                  </td>
                  <td>
                    <h1><xsl:value-of select="//model-info/name"/>&#160;<xsl:value-of select="//model-info/id"/>&#160; File Summary</h1>
                  </td>
                </tr>
              </table>
            </div>
            <div class="page_padding">
            
            
            <!--  FILE SUMMARY TABLE -->
              <h2 class="report_section" style="margin-bottom: 0">File Summary</h2>
              <table id="models" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//files/file) > 0">
                    <thead>
                      <tr>
                        <th>file name</th>
                        <th>line count</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//files/file">
                        <tr>
                          <td>
                            <xsl:value-of select="@file-name"/>
                          </td>
                          <td>
                            <xsl:value-of select="@count"/>
                          </td>
                        </tr>
                      </xsl:for-each>
                    </tbody>
                  </xsl:when>
                  <xsl:otherwise>
                    <tr><td><span class="note">None found.</span></td></tr>
                  </xsl:otherwise>
                </xsl:choose>
              </table>
              <xsl:if test="count(//files/file) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#models').dataTable( {
                          "aaSorting": [[1, 'asc']],
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "bPaginate": false,
                          "aoColumns": [
                              null,
                              null,
                          ]
                      } );
                  } );
                </script>
              </xsl:if>
              <br clear="all"/>
            </div>
            <div class="page_footer">
              Generated <xsl:value-of select="//report-meta/date"/> for build <xsl:value-of select="//report-meta/generator-params/build-id"/>
            </div>
          </div>
        </div>
      </body>
    </html>

  </xsl:template>

</xsl:stylesheet>
