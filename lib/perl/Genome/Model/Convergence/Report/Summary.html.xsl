<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
                xmlns:set="http://exslt.org/sets">

  <xsl:output method="html"/>
  <xsl:output encoding="utf-8"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:template match="/">

    <html>
      <head>
        <title>Summary: <xsl:value-of select="//model-info/name"/>&#160;<xsl:value-of select="//model-info/id"/></title>
        <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png" />
        <link rel="stylesheet" href="https://imp.gsc.wustl.edu/res/old/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen" />
        <style type="text/css">
        .page_footer {
            background-color: #f2f3e0;
            border-top: 1px solid #8dc643;
            padding: 5px;
        }
        /* This report is extra-wide */
        div.container {
            width: 1020px !important;
        }
        div.background {
            width: 1020px !important;
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
                    <h1><xsl:value-of select="//model-info/name"/>&#160;<xsl:value-of select="//model-info/id"/>&#160;Summary</h1>
                  </td>
                </tr>
              </table>
            </div>
            <div class="page_padding">
            
            
            <!--  SUMMARY TABLE -->
              <h2 class="report_section" style="margin-bottom: 0">Models</h2>
              <table id="models" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//members/member) > 0">
                    <thead>
                      <tr>
                        <th>model id</th>
                        <th>name</th>
                        <th>build id</th>
                        <th>date completed</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//members/member">
                        <tr>
                          <td>
                            <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/status.html?id=<xsl:value-of select="@model-id"/></xsl:attribute><xsl:value-of select="@model-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@name"/>
                          </td>
                          <td>
                             <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=<xsl:value-of select="@build-id"/></xsl:attribute><xsl:value-of select="@build-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@completed"/>
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
              <xsl:if test="count(//members/member) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#models').dataTable( {
                          "aaSorting": [[1, 'asc']],
						  "bAutoWidth": false,
						  "bStateSave": true,
						  "aoColumns": [
							  null,
							  null,
							  null,
							  null,
						  ]
					  } );
                  } );
                </script>
              </xsl:if>
              <br clear="all"/>
              
              <!-- ALIGNMENT METRICS -->
              <h2 class="report_section" style="margin-bottom: 0">Alignment Metrics</h2>
              <table id="metrics" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//metrics/build) > 0">
                    <thead>
                      <tr>
                        <th>build id</th>
                        <th>name</th>
                        <th>lanes</th>
                        <th>haploid coverage</th>
                        <th>input base count (kb)</th>
                        <th>unfiltered SNP calls</th>
                        <th>unfiltered Het %</th>
                        <th>unfiltered dbSNP concordance</th>
                        <th>filtered SNP calls</th>
                        <th>filtered Het %</th>
                        <th>filtered dbSNP concordance</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//metrics/build">
                        <tr>
                          <td>
                           <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=<xsl:value-of select="@build-id"/></xsl:attribute><xsl:value-of select="@build-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@name"/>
                          </td>
                          <td>
                            <xsl:value-of select="@lanes"/>
                          </td>
                          <td>
                            <xsl:value-of select="@haploid-coverage"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(@instrument-data-total-kb,'#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="@unfiltered-snp-calls"/>
                          </td>
                          <td>
                            <xsl:value-of select="@unfiltered-diploid-heterozygous-percentage"/>
                          </td>
                          <td>
                            <xsl:value-of select="@unfiltered-dbsnp-concordance"/>
                          </td>
                          <td>
                            <xsl:value-of select="@filtered-snp-calls"/>
                          </td>
                          <td>
                            <xsl:value-of select="@filtered-diploid-heterozygous-percentage"/>
                          </td>
                          <td>
                            <xsl:value-of select="@filtered-dbsnp-concordance"/>
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
              <xsl:if test="count(//metrics/build) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#metrics').dataTable( {
                          "aaSorting": [[1, 'asc']],
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "aoColumns": [
                              null,
                              null,
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last", "sType": "formatted-num" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last", "sType": "formatted-num" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                          ]
                      } );
                  } );
                </script>
              </xsl:if>
              <br clear="all"/>

              <!-- ALIGNMENT FLAGSTAT READ METRICS -->
              <h2 class="report_section" style="margin-bottom: 0">Alignment Flagstat Read Metrics</h2>
              <table id="metrics" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//metrics/build) > 0">
                    <thead>
                      <tr>
                        <th>build id</th>
                        <th>name</th>
                        <th>total reads</th>
                        <th>failing qc</th>
                        <th>marked duplicates</th>
                        <th>mapped</th>
                        <th>mapped %</th>
                        <th>paired in sequencing</th>
                        <th>marked as read1</th>
                        <th>marked as read2</th>
                        <th>mapped in proper pairs</th>
                        <th>mapped in proper pairs %</th>
                        <th>mapped in pair</th>
                        <th>mapped as singleton</th>
                        <th>mapped as singleton %</th>
                        <th>mapped in interchromosomal pairs</th>
                        <th>hq reads mapped in interchromosomal pairs</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//metrics/build">
                        <tr>
                          <td>
                           <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=<xsl:value-of select="@build-id"/></xsl:attribute><xsl:value-of select="@build-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@name"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-total_reads"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_marked_failing_qc"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_marked_duplicates"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_percentage"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_paired_in_sequencing"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_marked_as_read1"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_marked_as_read2"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_in_proper_pairs"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_in_proper_pairs_percentage"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_in_pair"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_as_singleton"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_as_singleton_percentage"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-reads_mapped_in_interchromosomal_pairs"/>
                          </td>
                          <td>
                            <xsl:value-of select="@flagstat-hq_reads_mapped_in_interchromosomal_pairs"/>
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
              <xsl:if test="count(//metrics/build) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#metrics').dataTable( {
                          "aaSorting": [[1, 'asc']],
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "aoColumns": [
                              null,
                              null,
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                          ]
                      } );
                  } );
                </script>
              </xsl:if>
              <br clear="all"/>

              <!-- SOMATIC SNP METRICS -->
              <h2 class="report_section" style="margin-bottom: 0">Somatic SNV Metrics</h2>
              <table id="snp-metrics" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//somatic-stats/build) > 0">
                    <thead>
                      <tr>
                        <th>build id</th>
                        <th>name</th>
                        <th>sniper</th>
						<th>filter</th>
						<th>ceu yri</th>
						<th>dbsnp</th>
						<th>loh</th>
						<th>loh fail</th>
						<th>annotate</th>
						<th>ucsc</th>
						<th>ucsc unanno.</th>
						<th>tier 1 (hc)</th>
						<th>tier 2 (hc)</th>
						<th>tier 3 (hc)</th>
						<th>tier 4 (hc)</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//somatic-stats/build">
                        <tr>
                          <td>
                            <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=<xsl:value-of select="@build-id"/></xsl:attribute><xsl:value-of select="@build-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@name"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='sniper_snp_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='snp_filter_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='filter_ceu_yri_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='dbsnp_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='loh_output_file']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='loh_fail_output_file']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='annotate_output_snp']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='ucsc_output_snp']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='ucsc_unannotated_output_snp']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_1_snp_file']/@count, '#,##0')"/> (<xsl:value-of select="format-number(files/file[@file-name='tier_1_snp_high_confidence_file']/@count, '#,##0')"/>)
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_2_snp_file']/@count, '#,##0')"/> (<xsl:value-of select="format-number(files/file[@file-name='tier_2_snp_high_confidence_file']/@count, '#,##0')"/>)
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_3_snp_file']/@count, '#,##0')"/> (<xsl:value-of select="format-number(files/file[@file-name='tier_3_snp_high_confidence_file']/@count, '#,##0')"/>)
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_4_snp_file']/@count, '#,##0')"/> (<xsl:value-of select="format-number(files/file[@file-name='tier_4_snp_high_confidence_file']/@count, '#,##0')"/>)
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
              <xsl:if test="count(//somatic-stats/build) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#snp-metrics').dataTable( {
                          "aaSorting": [[1, 'asc']],
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "aoColumns": [
                              null,
                              null,
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                          ]
                      } );
                  } );
                </script>
              </xsl:if>
              <br clear="all"/>
              
              <!-- SOMATIC INDEL METRICS -->
              <h2 class="report_section" style="margin-bottom: 0">Somatic Indel Metrics</h2>
              <table id="indel-metrics" class="list display" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;">
                <xsl:choose>
                  <xsl:when test="count(//somatic-stats/build) > 0">
                    <thead>
                      <tr>
                        <th>build id</th>
                        <th>name</th>
                        <th>sniper</th>
                        <th>lib filter single</th>
                        <th>lib filter multi</th>
                        <th>annotate</th>
                        <th>ucsc</th>
                        <th>ucsc unanno.</th>
                        <th>tier 1</th>
                        <th>tier 2</th>
                        <th>tier 3</th>
                        <th>tier 4</th>
                      </tr>
                    </thead>
                    <tbody>
                      <xsl:for-each select="//somatic-stats/build">
                        <tr>
                          <td>
                            <a><xsl:attribute name="href">https://imp.gsc.wustl.edu/view/Genome/Model/Build/status.html?id=<xsl:value-of select="@build-id"/></xsl:attribute><xsl:value-of select="@build-id"/></a>
                          </td>
                          <td>
                            <xsl:value-of select="@name"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='sniper_indel_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='indel_lib_filter_single_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='indel_lib_filter_multi_output']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='annotate_output_indel']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='ucsc_output_indel']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='ucsc_unannotated_output_indel']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_1_indel_file']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_2_indel_file']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_3_indel_file']/@count, '#,##0')"/>
                          </td>
                          <td>
                            <xsl:value-of select="format-number(files/file[@file-name='tier_4_indel_file']/@count, '#,##0')"/>
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
              <xsl:if test="count(//somatic-stats/build) > 0">
                <script type="text/javascript" charset="utf-8">
                  $(document).ready( function() {
                      $('#indel-metrics').dataTable( {
                          "aaSorting": [[1, 'asc']],
                          "bAutoWidth": false,
                          "bStateSave": true,
                          "aoColumns": [
                              null,
                              null,
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
                              { "sClass": "last" },
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
