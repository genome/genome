<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  
  <xsl:template name="genome_model_somatic" match="object[./types[./isa[@type='Genome::Model::Somatic']]]">
    <xsl:variable name="build_directory_url">
      <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='last_complete_build']/object/aspect[@name='data_directory']/value)" />
    </xsl:variable>
    
    <xsl:call-template name="genome_model">
	  <xsl:with-param name="build_directory_url">
	    <xsl:value-of select="$build_directory_url"/>
	  </xsl:with-param>
	  <xsl:with-param name="summary_report_url">
	    <xsl:value-of select="$build_directory_url"/><xsl:text>/cancer_report.html</xsl:text>
	  </xsl:with-param>
    </xsl:call-template>
  </xsl:template>

</xsl:stylesheet>
