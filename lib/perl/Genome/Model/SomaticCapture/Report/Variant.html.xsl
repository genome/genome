<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
                xmlns:set="http://exslt.org/sets">

  <xsl:output method="html"/>
  <xsl:output encoding="utf-8"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:template match="/">

    <html>
      <head>
        <title>Cancer Summary: <xsl:value-of select="//individual/@name"/>&#160;<xsl:value-of select="//individual/@gender"/>&#160;<xsl:value-of select="//individual/@id"/></title>
        <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png" />
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen" />
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/jquery/jquery.js"></script>
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/cancer_card/css/zoom2.css" type="text/css" media="screen"></link>
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/cancer_card/css/cancer_card.css" type="text/css" media="screen"></link>

        <!-- initialize circos graph zoom -->
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/cancer_card/js/dom-drag.js"></script>
        <script type="text/javascript" src="https://imp.gsc.wustl.edu/res/old/report_resources/cancer_card/js/HotSpot2.js"></script>
        <script type="text/javascript">
          addEvent(window, 'load', function() {
              HotSpotController.init("zoomImage",300, '<xsl:value-of select="//individual/circos-images/@large"/>','ZTbutton'); });

              function addEvent(obj, evType, fn) {
                  if (obj.addEventListener) {
                      obj.addEventListener(evType, fn, false);
                      return true;
                  } else if (obj.attachEvent) {
                      var r = obj.attachEvent("on" + evType, fn);
                      return r;
                  } else {
                      return false;
                  }
              }
        </script>

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
                    <h1><xsl:value-of select="//individual/@name"/>&#160;<xsl:value-of select="//individual/@gender"/>&#160;<xsl:value-of select="//individual/@id"/>&#160;Cancer Summary</h1>
                  </td>
                </tr>
              </table>
            </div>
            <div class="page_padding">
              <table cellpadding="0" cellspacing="0" border="0" class="info_table_group">
                <tr>
                  <td>
                    <h3 class="group_header" style="margin-top: 0;">Clinical Data</h3>
                    <table border="0" cellpadding="0" cellspacing="0">
                      <tr>
                        <td>
                          <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                            <colgroup>
                              <col/>
                              <col width="100%"/>
                            </colgroup>
                            <tr><td class="label">Name:</td><td class="value"><xsl:value-of select="//individual/@name"/></td></tr>
                            <tr>
                              <td class="label">Gender:</td>
                              <td class="value">
                                <xsl:choose>
                                  <xsl:when test="string(//individual/@gender)">
                                    <xsl:value-of select="//individual/@gender"/>
                                  </xsl:when>
                                  <xsl:otherwise>
                                    Not Provided
                                  </xsl:otherwise>
                                </xsl:choose>
                              </td>
                            </tr>
                            <tr><td class="label">ID:</td><td class="value"><xsl:value-of select="//individual/@id"/></td></tr>

                          </table>
                        </td>
                        <td>
                          <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%" style="float: left;">
                            <colgroup>
                              <col/>
                              <col width="100%"/>
                            </colgroup>

                            <tr><td class="label">Year Diagnosed:</td><td class="value"><xsl:value-of select="//clinical-data/@diagnosis-year"/></td></tr>
                            <tr><td class="label">Diagnosed at Age:</td><td class="value"><xsl:value-of select="//clinical-data/@diagnosis-age"/></td></tr>
                            <xsl:choose>
                              <xsl:when test="//clinical-data/@alive = '1'">
                                <tr><td class="label">Survived:</td><td class="value">Yes</td></tr>
                              </xsl:when>
                              <xsl:otherwise>
                                <tr><td class="label">Survived:</td><td class="value">No</td></tr>
                                <tr><td class="label">Days Survived:</td><td class="value"><xsl:value-of select="//clinical-data/@days-survived"/></td></tr>
                              </xsl:otherwise>
                            </xsl:choose>
                          </table>
                        </td>
                        <td>
                          <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%" style="float: left;">
                            <colgroup>
                              <col/>
                              <col width="100%"/>
                            </colgroup>

                            <tr><td class="label">Treatment:</td><td class="value"><xsl:value-of select="//clinical-data/@treatment"/></td></tr>
                            <tr><td class="label">Outcome:</td><td class="value"><xsl:value-of select="//clinical-data/@outcome"/></td></tr>
                            <tr><td class="label">AMP:</td><td class="value"><xsl:value-of select="//clinical-data/@amp"/></td></tr>

                          </table>
                        </td>
                      </tr>
                    </table>
                  </td>
                  <td>
                    <h3 class="group_header" style="margin-top: 0;">Sequencing Stats</h3>
                    <table border="0" cellpadding="0" cellspacing="0" class="info_table" width="100%">
                      <colgroup>
                        <col/>
                        <col width="100%"/>
                      </colgroup>
                      <tr><td class="label">Normal Coverage:</td><td class="value"><xsl:value-of select="//samples/sample/models/model/@normal-haploid-coverage"/>X</td></tr>
                      <tr><td class="label">Tumor Coverage:</td><td class="value"><xsl:value-of select="//samples/sample/models/model/@tumor-haploid-coverage"/>X</td></tr>
                      <tr><td class="label">Purity Estimate:</td><td class="value"><xsl:value-of select="//samples/sample/models/model/@purity-estimate"/></td></tr>
                      <xsl:if test="string(//samples/sample/models/model/@ploidy-estimate)">
                        <tr><td class="label">Ploidy Estimate:</td><td class="value"><xsl:value-of select="//samples/sample/models/model/@ploidy-estimate"/></td></tr>
                      </xsl:if>
                    </table>
                  </td>
                </tr>
              </table>
              <hr style="margin-bottom: 0;"/>
              <h2 class="report_section" style="margin-bottom: 0; margin-top: 0">Somatic Genome Graph&#160;&#160;&#160;<a id="ZTbutton" href="javascript: void(0);" style="font-size: 85%; font-weight: normal;">[toggle zoom]</a></h2>
              <p id="ZTthumbnail">
                <img>
                  <xsl:attribute name="id">zoomImage</xsl:attribute>
                  <xsl:attribute name="src"><xsl:value-of select="//individual/circos-images/@small"/></xsl:attribute>
                  <xsl:attribute name="width">920</xsl:attribute>
                  <xsl:attribute name="height">920</xsl:attribute>
                </img>
              </p>
              <h2 class="report_section" style="margin-bottom: 0">Tier 1 SNVs</h2>
              <table id="tier_1_snps" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//variants/snps/snp) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th>chromosome</th>
                        <th class="last">start</th>
                        <th class="last">reference</th>
                        <th class="last">variant</th>
                        <th class="last">gene</th>
                        <th class="last">amino acid change</th>
                        <th class="last">trv type</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//variants/snps/snp">
                        <xsl:sort select="@validation-status" data-type="text" order="ascending"/>
                        <xsl:sort select="@chromosome" data-type="number" order="ascending"/>
                        <tr>
                          <td class="validation_status">
                            <!-- expand validation status codes -->
                            <xsl:variable select="@validation-status" name="validation-status"/>
                            <xsl:choose>
                              <xsl:when test="$validation-status = 'G'">
                                germline
                              </xsl:when>
                              <xsl:when test="$validation-status = 'WT'">
                                wildtype
                              </xsl:when>
                              <xsl:when test="$validation-status = 'O'">
                                other
                              </xsl:when>
                              <xsl:when test="$validation-status = 'LQ'">
                                low quality
                              </xsl:when>
                              <xsl:when test="$validation-status = 'S'">
                                somatic
                              </xsl:when>
                              <xsl:when test="$validation-status = 'V'">
                                variant
                              </xsl:when>
                              <xsl:when test="$validation-status = 'A'">
                                ambiguous
                              </xsl:when>
                              <xsl:when test="$validation-status = 'NC'">
                                no coverage
                              </xsl:when>
                              <xsl:when test="$validation-status = 'P'">
                                predicted
                              </xsl:when>
                              <xsl:otherwise>
                                <xsl:value-of select="@validation-status"/>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td><xsl:value-of select="@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="start" select="@start"/><xsl:value-of select="format-number($start, '#,##0')"/>
                          </td>
                          <td class="last"><xsl:value-of select="@reference-allele"/></td>
                          <td class="last"><xsl:value-of select="@variant-allele"/></td>
                          <td class="last"><xsl:value-of select="@gene"/></td>
                          <td class="last"><xsl:value-of select="@amino-acid-change"/></td>
                          <td class="last"><xsl:value-of select="@trv-type"/></td>
                        </tr>
                      </xsl:for-each>
                    </tbody>
                  </xsl:when>
                  <xsl:otherwise>
                    <tr><td><span class="note">None found.</span></td></tr>
                  </xsl:otherwise>
                </xsl:choose>
              </table>
              <xsl:if test="count(//variants/snps/snp) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#tier_1_snps').dataTable( {
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  null,
							  null,
							  null,
							  null,
							  null
						  ]
					  } );
                  } );
                </script>
              </xsl:if>

              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0">Tier 1 Insertions</h2>
              <table id="tier_1_insertions" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//variants/insertions/insertion) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th>chromosome</th>
                        <th class="last">start</th>
                        <th class="last">stop</th>
                        <th class="last">variant</th>
                        <th class="last">gene</th>
                        <th class="last">amino acid change</th>
                        <th class="last">trv type</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//variants/insertions/insertion">
                        <xsl:sort select="@validation-status" data-type="text" order="ascending"/>
                        <xsl:sort select="@chromosome" data-type="number" order="ascending"/>
                        <tr>
                          <td class="validation_status">
                            <!-- expand validation status codes -->
                            <xsl:variable select="@validation-status" name="validation-status"/>
                            <xsl:choose>
                              <xsl:when test="$validation-status = 'G'">
                                germline
                              </xsl:when>
                              <xsl:when test="$validation-status = 'WT'">
                                wildtype
                              </xsl:when>
                              <xsl:when test="$validation-status = 'O'">
                                other
                              </xsl:when>
                              <xsl:when test="$validation-status = 'LQ'">
                                low quality
                              </xsl:when>
                              <xsl:when test="$validation-status = 'S'">
                                somatic
                              </xsl:when>
                              <xsl:when test="$validation-status = 'V'">
                                variant
                              </xsl:when>
                              <xsl:when test="$validation-status = 'A'">
                                ambiguous
                              </xsl:when>
                              <xsl:when test="$validation-status = 'NC'">
                                no coverage
                              </xsl:when>
                              <xsl:when test="$validation-status = 'P'">
                                predicted
                              </xsl:when>
                              <xsl:otherwise>
                                <xsl:value-of select="@validation-status"/>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td><xsl:value-of select="@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="start" select="@start"/><xsl:value-of select="format-number($start, '#,##0')"/>
                          </td>
                          <td class="last">
                            <xsl:variable name="stop" select="@stop"/><xsl:value-of select="format-number($stop, '#,##0')"/>
                          </td>
                          <td class="last"><xsl:value-of select="@variant-allele"/></td>
                          <td class="last"><xsl:value-of select="@gene"/></td>
                          <td class="last"><xsl:value-of select="@amino-acid-change"/></td>
                          <td class="last"><xsl:value-of select="@trv-type"/></td>
                        </tr>
                      </xsl:for-each>
                    </tbody>
                  </xsl:when>
                  <xsl:otherwise>
                    <tr><td><span class="note">None found.</span></td></tr>
                  </xsl:otherwise>
                </xsl:choose>
              </table>
              <xsl:if test="count(//variants/insertions/insertion) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#tier_1_insertions').dataTable({
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
							  null,
							  null,
							  null,
							  null
						  ]
					  });
                  } );
                </script>
              </xsl:if>

              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0">Tier 1 Deletions</h2>
              <table id="tier_1_deletions" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//variants/deletions/deletion) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th>chromosome</th>
                        <th class="last">start</th>
                        <th class="last">stop</th>
                        <th class="last">reference</th>
                        <th class="last">gene</th>
                        <th class="last">amino acid change</th>
                        <th class="last">trv type</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//variants/deletions/deletion">
                        <xsl:sort select="@validation-status" data-type="text" order="ascending"/>
                        <xsl:sort select="@chromosome" data-type="number" order="ascending"/>
                        <tr>
                          <td class="validation_status">
                            <!-- expand validation status codes -->
                            <xsl:variable select="@validation-status" name="validation-status"/>
                            <xsl:choose>
                              <xsl:when test="$validation-status = 'G'">
                                germline
                              </xsl:when>
                              <xsl:when test="$validation-status = 'WT'">
                                wildtype
                              </xsl:when>
                              <xsl:when test="$validation-status = 'O'">
                                other
                              </xsl:when>
                              <xsl:when test="$validation-status = 'LQ'">
                                low quality
                              </xsl:when>
                              <xsl:when test="$validation-status = 'S'">
                                somatic
                              </xsl:when>
                              <xsl:when test="$validation-status = 'V'">
                                variant
                              </xsl:when>
                              <xsl:when test="$validation-status = 'A'">
                                ambiguous
                              </xsl:when>
                              <xsl:when test="$validation-status = 'NC'">
                                no coverage
                              </xsl:when>
                              <xsl:when test="$validation-status = 'P'">
                                predicted
                              </xsl:when>
                              <xsl:otherwise>
                                <xsl:value-of select="@validation-status"/>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td><xsl:value-of select="@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="start" select="@start"/><xsl:value-of select="format-number($start, '#,##0')"/>
                          </td>
                          <td class="last">
                            <xsl:variable name="stop" select="@stop"/><xsl:value-of select="format-number($stop, '#,##0')"/>
                          </td>
                          <td class="last"><xsl:value-of select="@reference-allele"/></td>
                          <td class="last"><xsl:value-of select="@gene"/></td>
                          <td class="last"><xsl:value-of select="@amino-acid-change"/></td>
                          <td class="last"><xsl:value-of select="@trv-type"/></td>
                        </tr>
                      </xsl:for-each>
                    </tbody>
                  </xsl:when>
                  <xsl:otherwise>
                    <tr><td><span class="note">None found.</span></td></tr>
                  </xsl:otherwise>
                </xsl:choose>
              </table>
              <xsl:if test="count(//variants/deletions/deletion) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#tier_1_deletions').dataTable({
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
							  null,
							  null,
							  null,
							  null
						  ]
					  });
                  } );
                </script>
              </xsl:if>

              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0;">Structural Variations (translocations)</h2>
              <table id="sv_translocations" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//structural-variants/translocations/translocation) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th class="last">chromosome</th>
                        <th class="last">position</th>
                        <th>&#160;</th>
                        <th class="last">chromosome</th>
                        <th class="last">position</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//structural-variants/translocations/translocation">
                        <xsl:sort select="@validation-status" data-type="text" order="ascending"/>
                        <xsl:sort select="start/@chromosome" data-type="number" order="ascending"/>
                        <tr>
                          <td class="validation_status">
                            <!-- expand validation status codes -->
                            <xsl:variable select="@validation-status" name="validation-status"/>
                            <xsl:choose>
                              <xsl:when test="$validation-status = 'G'">
                                germline
                              </xsl:when>
                              <xsl:when test="$validation-status = 'WT'">
                                wildtype
                              </xsl:when>
                              <xsl:when test="$validation-status = 'O'">
                                other
                              </xsl:when>
                              <xsl:when test="$validation-status = 'LQ'">
                                low quality
                              </xsl:when>
                              <xsl:when test="$validation-status = 'S'">
                                somatic
                              </xsl:when>
                              <xsl:when test="$validation-status = 'V'">
                                variant
                              </xsl:when>
                              <xsl:when test="$validation-status = 'A'">
                                ambiguous
                              </xsl:when>
                              <xsl:when test="$validation-status = 'NC'">
                                no coverage
                              </xsl:when>
                              <xsl:when test="$validation-status = 'P'">
                                predicted
                              </xsl:when>
                              <xsl:otherwise>
                                <xsl:value-of select="@validation-status"/>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td class="last"><xsl:value-of select="start/@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="start_position" select="start/@position"/><xsl:value-of select="format-number($start_position, '#,##0')"/>
                          </td>

                          <td class="last"><span style="font-size: 100%; font-weight: bold;">&#8594;</span></td>

                          <td class="last"><xsl:value-of select="stop/@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="stop_position" select="stop/@position"/><xsl:value-of select="format-number($stop_position, '#,##0')"/>
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
              <xsl:if test="count(//structural-variants/translocations/translocation) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#sv_translocations').dataTable({
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  null,
							  null,
							  { "sType": "formatted-num" },
						  ]
					  });
                  } );
                </script>
              </xsl:if>

              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0;">Structural Variations (insertions)</h2>
              <table id="sv_insertions" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                
                <xsl:choose>
                  <xsl:when test="count(//structural-variants/insertions/insertion) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th class="last">chromosome</th>
                        <th class="last">start</th>
                        <th class="last">stop</th>
                        <th class="last">size</th>
                      </tr>
                    </thead>
                    <tbody>
                    <xsl:for-each select="//structural-variants/insertions/insertion">
                      <xsl:sort select="@validation-status" data-type="number" order="ascending"/>
                      <xsl:sort select="start/@chromosome" data-type="number" order="ascending"/>
                      <tr>
                        <td class="validation_status">
                          <!-- expand validation status codes -->
                          <xsl:variable select="@validation-status" name="validation-status"/>
                          <xsl:choose>
                            <xsl:when test="$validation-status = 'G'">
                              germline
                            </xsl:when>
                            <xsl:when test="$validation-status = 'WT'">
                              wildtype
                            </xsl:when>
                            <xsl:when test="$validation-status = 'O'">
                              other
                            </xsl:when>
                            <xsl:when test="$validation-status = 'LQ'">
                              low quality
                            </xsl:when>
                            <xsl:when test="$validation-status = 'S'">
                              somatic
                            </xsl:when>
                            <xsl:when test="$validation-status = 'V'">
                              variant
                            </xsl:when>
                            <xsl:when test="$validation-status = 'A'">
                              ambiguous
                            </xsl:when>
                            <xsl:when test="$validation-status = 'NC'">
                              no coverage
                            </xsl:when>
                            <xsl:when test="$validation-status = 'P'">
                              predicted
                            </xsl:when>
                            <xsl:otherwise>
                              <xsl:value-of select="@validation-status"/>
                            </xsl:otherwise>
                          </xsl:choose>
                        </td>
                        <td class="last"><xsl:value-of select="start/@chromosome"/></td>
                        <td class="last">
                          <xsl:variable name="start_position" select="start/@position"/><xsl:value-of select="format-number($start_position, '#,##0')"/>
                        </td>
                        <td class="last">
                          <xsl:variable name="stop_position" select="stop/@position"/><xsl:value-of select="format-number($stop_position, '#,##0')"/>
                        </td>
                        <td class="last">
                          <xsl:variable name="size" select="@size"/><xsl:value-of select="format-number($size, '#,##0')"/>
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
              <xsl:if test="count(//structural-variants/insertions/insertion) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#sv_insertions').dataTable({
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
						  ]
					  });
                  } );
                </script>
              </xsl:if>

              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0;">Structural Variations (deletions)</h2>
              <table id="sv_deletions" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//structural-variants/deletions/deletion) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th class="last">chromosome</th>
                        <th class="last">start</th>
                        <th class="last">stop</th>
                        <th class="last">size</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//structural-variants/deletions/deletion">
                        <xsl:sort select="@validation-status" data-type="number" order="ascending"/>
                        <xsl:sort select="start/@chromosome" data-type="number" order="ascending"/>
                        <tr>
                          <td class="validation_status">
                            <!-- expand validation status codes -->
                            <xsl:variable select="@validation-status" name="validation-status"/>
                            <xsl:choose>
                              <xsl:when test="$validation-status = 'G'">
                                germline
                              </xsl:when>
                              <xsl:when test="$validation-status = 'WT'">
                                wildtype
                              </xsl:when>
                              <xsl:when test="$validation-status = 'O'">
                                other
                              </xsl:when>
                              <xsl:when test="$validation-status = 'LQ'">
                                low quality
                              </xsl:when>
                              <xsl:when test="$validation-status = 'S'">
                                somatic
                              </xsl:when>
                              <xsl:when test="$validation-status = 'V'">
                                variant
                              </xsl:when>
                              <xsl:when test="$validation-status = 'A'">
                                ambiguous
                              </xsl:when>
                              <xsl:when test="$validation-status = 'NC'">
                                no coverage
                              </xsl:when>
                              <xsl:when test="$validation-status = 'P'">
                                predicted
                              </xsl:when>
                              <xsl:otherwise>
                                <xsl:value-of select="@validation-status"/>
                              </xsl:otherwise>
                            </xsl:choose>
                          </td>
                          <td class="last"><xsl:value-of select="start/@chromosome"/></td>
                          <td class="last">
                            <xsl:variable name="start_position" select="start/@position"/><xsl:value-of select="format-number($start_position, '#,##0')"/>
                          </td>
                          <td class="last">
                            <xsl:variable name="stop_position" select="stop/@position"/><xsl:value-of select="format-number($stop_position, '#,##0')"/>
                          </td>
                          <td class="last">
                            <xsl:variable name="size" select="@size"/><xsl:value-of select="format-number($size, '#,##0')"/>
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
              <xsl:if test="count(//structural-variants/deletions/deletion) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#sv_deletions').dataTable( {
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
							  { "sType": "formatted-num" },
						  ]
					  } );
                  } );
                </script>
              </xsl:if>
              
              <br clear="all"/>

              <h2 class="report_section" style="margin-bottom: 0;">Structural Variations (inversions)</h2>
              <table id="sv_inversions" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                
                <xsl:choose>
                  <xsl:when test="count(//structural-variants/inversions/inversion) > 0">
                    <thead>
                      <tr>
                        <th>validation status</th>
                        <th class="last">chromosome</th>
                        <th class="last">start</th>
                        <th class="last">stop</th>
                        <th class="last">size</th>
                      </tr>
                    </thead>
                    <tbody>
                    <xsl:for-each select="//structural-variants/inversions/inversion">
                      <xsl:sort select="@validation-status" data-type="number" order="ascending"/>
                      <xsl:sort select="start/@chromosome" data-type="number" order="ascending"/>
                      <tr>
                        <td class="validation_status">
                          <!-- expand validation status codes -->
                          <xsl:variable select="@validation-status" name="validation-status"/>
                          <xsl:choose>
                            <xsl:when test="$validation-status = 'G'">
                              germline
                            </xsl:when>
                            <xsl:when test="$validation-status = 'WT'">
                              wildtype
                            </xsl:when>
                            <xsl:when test="$validation-status = 'O'">
                              other
                            </xsl:when>
                            <xsl:when test="$validation-status = 'LQ'">
                              low quality
                            </xsl:when>
                            <xsl:when test="$validation-status = 'S'">
                              somatic
                            </xsl:when>
                            <xsl:when test="$validation-status = 'V'">
                              variant
                            </xsl:when>
                            <xsl:when test="$validation-status = 'A'">
                              ambiguous
                            </xsl:when>
                            <xsl:when test="$validation-status = 'NC'">
                              no coverage
                            </xsl:when>
                            <xsl:when test="$validation-status = 'P'">
                              predicted
                            </xsl:when>
                            <xsl:otherwise>
                              <xsl:value-of select="@validation-status"/>
                            </xsl:otherwise>
                          </xsl:choose>
                        </td>
                        <td class="last"><xsl:value-of select="start/@chromosome"/></td>
                        <td class="last">
                          <xsl:variable name="start_position" select="start/@position"/><xsl:value-of select="format-number($start_position, '#,##0')"/>
                        </td>
                        <td class="last">
                          <xsl:variable name="stop_position" select="stop/@position"/><xsl:value-of select="format-number($stop_position, '#,##0')"/>
                        </td>
                        <td class="last">
                          <xsl:variable name="size" select="@size"/><xsl:value-of select="format-number($size, '#,##0')"/>
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
              <xsl:if test="count(//structural-variants/inversions/inversion) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#sv_inversions').dataTable({
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "aoColumns": [
                              null,
                              null,
                              { "sType": "formatted-num" },
                              { "sType": "formatted-num" },
                              { "sType": "formatted-num" },
                          ]
                      });
                  } );
                </script>
              </xsl:if>

              <br clear="all" />
            </div>
          </div>
        </div>
      </body>
    </html>

  </xsl:template>

</xsl:stylesheet>
