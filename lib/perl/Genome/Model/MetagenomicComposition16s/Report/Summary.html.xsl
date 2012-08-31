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
    </head>
  <body>
    <div class="container"><div class="background">
      <h1 class="page_title"><xsl:value-of select="//report-meta/description"/></h1>
      <div class="page_padding">
        
      <h2 class="report_section">Model and Build Information</h2>
      <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
      <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
          <col width="30%"/>
          <col width="70%"/>
				</colgroup>
				<tr>
          <td class="label" width="25%">Build Id</td>
				  <td class="value"><xsl:value-of select="//model-info/build-id"/></td>
				</tr>
				<tr>
          <td class="label">Data Directory</td>
          <td class="value"><a><xsl:attribute name="href">https://gscweb.gsc.wustl.edu/<xsl:value-of select="//model-info/data-directory"/></xsl:attribute><xsl:text>View Directory</xsl:text></a></td>
				</tr>
				<tr>
          <td class="label" width="25%">Model Id</td>
				  <td class="value"><xsl:value-of select="//model-info/id"/></td>
				</tr>
				<tr>
          <td class="label" width="25%">Name</td>
				  <td class="value"><xsl:value-of select="//model-info/name"/></td>
				</tr>
				<tr>
          <td class="label">Subject</td>
				  <td class="value"><xsl:value-of select="//model-info/subject-name"/></td>
				</tr>
				<tr>
          <td class="label">Processing Profile</td>
				  <td class="value"><xsl:value-of select="//model-info/processing-profile-name"/></td>
				</tr>
				<tr>
          <td class="label">Seq Platform</td>
				  <td class="value"><xsl:value-of select="//model-info/sequencing-platform"/></td>
				</tr>
				<tr>
          <td class="label">Seq Center</td>
				  <td class="value"><xsl:value-of select="//model-info/sequencing-center"/></td>
				</tr>
				<tr>
          <td class="label">Min Amplicon Size</td>
				  <td class="value"><xsl:value-of select="//model-info/amplicon-size"/></td>
				</tr>
				<tr>
          <td class="label">Exclude Contminated Amps</td>
          <td class="value"><xsl:value-of select="//model-info/exclude-contaminated-amplicons"/></td>
				</tr>
				<tr>
          <td class="label">Use Latest Read Iteration</td>
          <td class="value"><xsl:value-of select="//model-info/only-use-latest-iteration-of-reads"/></td>
				</tr>
				<tr>
          <td class="label">Assembler</td>
				  <td class="value"><xsl:value-of select="//model-info/assembler"/></td>
				</tr>
				<tr>
          <td class="label">Assembler Params</td>
          <td class="value"><xsl:value-of select="//model-info/assembler-params"/></td>
				</tr>
				<tr>
          <td class="label">Trimmer</td>
				  <td class="value"><xsl:value-of select="//model-info/trimmer"/></td>
				</tr>
				<tr>
          <td class="label">Trimmer Params</td>
          <td class="value"><xsl:value-of select="//model-info/trimmer-params"/></td>
				</tr>
				<tr>
          <td class="label">Classifier</td>
				  <td class="value"><xsl:value-of select="//model-info/classifier"/></td>
				</tr>
				<tr>
          <td class="label">Classifier Params</td>
          <td class="value"><xsl:value-of select="//model-info/classifier-params"/></td>
				</tr>
			  </table>
			</td>
      </tr>
      </table>

      <h2 class="report_section">Amplicon Statistics</h2>
      <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
      <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
          <col width="30%"/>
          <col width="70%"/>
				</colgroup>
				<tr>
				  <td class="label">Attempted</td>
				  <td class="value"><xsl:value-of select="//stat/amplicons-attempted"/></td>
				</tr>
				<tr>
          <td class="label">Processed</td>
          <td class="value"><xsl:value-of select="//stat/amplicons-processed"/></td>
				</tr>
				<tr>
          <xsl:variable name="processed-success" select="//stat/amplicons-processed-success"/>
          <td class="label">Processed Success</td>
          <td class="value"><xsl:value-of select="$processed-success * 100"/>%</td>
				</tr>
				<tr>
          <td class="label">Classified</td>
				  <td class="value"><xsl:value-of select="//stat/amplicons-classified"/></td>
				</tr>
				<tr>
          <xsl:variable name="classified-success" select="//stat/amplicons-classified-success"/>
          <td class="label">Classified Success</td>
          <td class="value"><xsl:value-of select="$classified-success * 100"/>%</td>
				</tr>
				<tr>
          <td class="label">Classification Error</td>
          <td class="value"><xsl:value-of select="//stat/amplicons-classification-error"/></td>
				</tr>
				<tr>
          <td class="label">Average Length w/o Primer</td>
				  <td class="value"><xsl:value-of select="//stat/length-average"/></td>
				</tr>
				<tr>
          <td class="label">Max Length w/o Primer</td>
				  <td class="value"><xsl:value-of select="//stat/length-maximum"/></td>
				</tr>
				<tr>
          <td class="label">Median Length w/o Primer</td>
				  <td class="value"><xsl:value-of select="//stat/length-median"/></td>
				</tr>
				<tr>
          <td class="label">Min Length w/o Primer</td>
				  <td class="value"><xsl:value-of select="//stat/length-minimum"/></td>
				</tr>
			  </table>
			</td>
      </tr>
      </table>

      <h2 class="report_section">Report Info</h2>
      <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
      <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
          <col width="30%"/>
          <col width="70%"/>
				</colgroup>
				<tr>
				  <td class="label">Date Generated</td>
          <td class="value"><xsl:value-of select="//report-meta/date"/></td>
				</tr>
				<tr>
          <td class="label">Generator</td>
          <td class="value"><xsl:value-of select="//report-meta/generator"/></td>
				</tr>
			  </table>
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
