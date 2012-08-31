<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
  <xsl:variable name="reads_attempted" select="//aspect[@name='reads_attempted']/value"/>
  <xsl:variable name="reads_processed" select="//aspect[@name='reads_processed']/value"/>
  <xsl:variable name="reads_assembled" select="//aspect[@name='reads_assembled']/value"/>
  <xsl:variable name='reads_not_assembled' select="//aspect[@name='reads_not_assembled']/value"/>
  <xsl:template match="/">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
    <title><xsl:value-of select="//report-meta/description"/></title> 
    <link href="layout.css" rel="stylesheet" type="text/css"></link>
    <link rel="shortcut icon" href="https://imp.gsc.wustl.edu/static/report_resources/apipe_dashboard/images/gc_favicon.png" type="image/png"/>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/static/report_resources/apipe_dashboard/css/master.css" type="text/css" media="screen"/>
    </head>
  <body>
    <div class="container"><div class="background">
      <h1 class="page_title">Summary Report</h1>
      <div class="page_padding">
        
      <h2 class="report_section">Assembler Information</h2>
      <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
      <tr>
      <td>
        <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
        <colgroup>
          <col width="30%"/>
          <col width="70%"/>
        </colgroup>
				<tr>
          <td class="label">Subject</td>
          <td class="value"><xsl:value-of select="//aspect[@name='subject_name']/value"/></td>
				</tr>
				<tr>
          <td class="label">Assembler Name</td>
          <td class="value"><xsl:value-of select="//aspect[@name='assembler_name']/value"/></td>
				</tr>
				<tr>
          <td class="label">Assembler Version</td>
          <td class="value"><xsl:value-of select="//aspect[@name='assembler_version']/value"/></td>
				</tr>
				<tr>
          <td class="label">Assembler Params</td>
          <td class="value"><xsl:value-of select="//aspect[@name='assembler_params']/value"/></td>
				</tr>
				<tr>
          <td class="label">Assembler Kmer</td>
          <td class="value"><xsl:value-of select="//aspect[@name='assembler_kmer']/value"/></td>
				</tr>
				<tr>
          <td class="label">Read Coverage Used</td>
          <td class="value"><xsl:value-of select="//aspect[@name='coverage']/value"/></td>
				</tr>
				<tr>
          <td class="label">Genome Size Used</td>
          <td class="value"><xsl:value-of select="//aspect[@name='genome_size']/value"/></td>
				</tr>
				<tr>
	        <td class="label">Average Insert Size</td>
          <td class="value"><xsl:value-of select="//aspect[@name='insert_size']/value"/></td>
				</tr>
			  </table>
			</td>
      </tr>
      </table>
      <h2 class="report_section">Assembly Metrics</h2>
      <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table_group">
      <tr>
			<td>
			  <table width="100%" cellspacing="0" cellpadding="0" border="0" class="info_table">
				<colgroup>
          <col width="30%"/>
          <col width="70%"/>
				</colgroup>
				<tr>
          <td class="label">Reads Attempted</td>
          <td class="value"><xsl:value-of select="$reads_attempted"/></td>
				</tr>
				<tr>
          <td class="label">Reads Processed</td>
          <td class="value"><xsl:value-of select="$reads_processed"/></td>
				</tr>
				<tr>
          <td class="label">Reads Processed Success</td>
          <td class="value"><xsl:value-of select="format-number(($reads_processed div $reads_attempted), '#.0%')"/></td>
				</tr>
				<tr>
          <td class="label">Reads Assembled</td>
          <td class="value"><xsl:value-of select="$reads_assembled"/></td>
				</tr>
	      <tr>
	        <td class="label">Reads Assembled Success</td>
          <td class="value"><xsl:value-of select="format-number(($reads_assembled div $reads_processed), '#.0%')"/></td>
	      </tr>
	      <tr>
	        <td class="label">Reads Not Assembled (Chaff rate)</td>
          <td class="value"><xsl:value-of select="//aspect[@name='reads']/value * 100"/><xsl:text>&#37;</xsl:text></td>
	      </tr>
			  <tr>
          <td class="label">Average Read Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='reads_processed_average_length']/value"/></td>
			  </tr>
			  <tr>
          <td class="label">Assembly Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='assembly_length']/value"/></td>
			  </tr>
			  <tr>
          <td class="label">Contigs</td>
          <td class="value"><xsl:value-of select="//aspect[@name='contigs_count']/value"/></td>
			  </tr>
			  <tr>
          <td class="label">Average Contig Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='contigs_average_length']/value"/></td>
			  </tr>
	      <tr>
          <td class="label">Average Major Contig Length (>=<xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/>)</td>
          <td class="value"><xsl:value-of select="//aspect[@name='contigs_major_average_length']/value"/></td>
	      </tr>
				<tr>
          <td class="label">n50 Contig Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='contigs_n50_length']/value"/></td>
				</tr>
	      <tr>
          <td class="label">n50 Major Contig Length (>=<xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/>)</td>
          <td class="value"><xsl:value-of select="//aspect[@name='contigs_major_n50_length']/value"/></td>
	      </tr>
				<tr>
          <td class="label">Supercontigs</td>
          <td class="value"><xsl:value-of select="//aspect[@name='supercontigs_count']/value"/></td>
				</tr>
				<tr>
          <td class="label">Average Supercontig Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='supercontigs_average_length']/value"/></td>
				</tr>
	      <tr>
           <td class="label">Average Major Supercontig Length (>=<xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/>)</td>
           <td class="value"><xsl:value-of select="//aspect[@name='supercontigs_major_average_length']/value"/></td>
	      </tr>
				<tr>
          <td class="label">n50 Supercontig Length</td>
          <td class="value"><xsl:value-of select="//aspect[@name='supercontigs_n50_length']/value"/></td>
				</tr>
	      <tr>
          <td class="label">n50 Major Supercontig Length (>=300)</td>
          <td class="value"><xsl:value-of select="//aspect[@name='supercontigs_major_n50_length']/value"/></td>
	      </tr>

	  <xsl:if test="//metric/major-contig-length='500'">
	    <tr>
              <td class="label">n50 Major Supercontig Length (>=500)</td>
              <td class="value"><xsl:value-of select="//metric/n50-supercontig-length-gt-500"/></td>
	    </tr>
	  </xsl:if>
	  
	  <xsl:if test="//metric/read-depths-ge-5x='NA'">
	    <tr>
	      <td class="label">Read Coverage (>= 5x)</td>
	      <td class="value"><xsl:value-of select="//metric/read-depths-ge-5x"/></td>
	    </tr>
	  </xsl:if>

	  <xsl:if test="//metric/read-depths-ge-5x!='NA'">
	    <tr>
	      <td class="label">Read Coverage (>= 5x)</td>
	      <td class="value"><xsl:value-of select="//metric/read-depths-ge-5x"/><xsl:text>&#37;</xsl:text></td>
	    </tr>
	  </xsl:if>
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
