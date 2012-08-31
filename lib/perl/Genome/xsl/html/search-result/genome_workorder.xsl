<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_workorder" match="object[./types[./isa[@type='Genome::WorkOrder']]]">
    <div class="search_result">
      <div class="result_icon genome_workorder_32">
        <br/>
      </div>
      <div class="result">
        <h3>Work Order:
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext">
            <xsl:value-of select="aspect[@name='name']/value"/>
          </xsl:with-param>
        </xsl:call-template>
        </h3>
        <p class="resource_buttons">
          <a class="mini btn"><xsl:attribute name="href">https://gscweb.gsc.wustl.edu/wiki/<xsl:value-of select="aspect[@name='name']"/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>wiki page</a>

          <a class="mini btn"><xsl:attribute name="href">http://linus222:8090/view/genome/search/query/status.html?query=<xsl:value-of select="aspect[@name='barcode']/value"/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>search for barcode <xsl:value-of select="aspect[@name='barcode']/value"/></a>

        </p>
        <p class="result_summary">
          <strong>Pipeline: </strong> <xsl:value-of select="aspect[@name='pipeline']/value"/><xsl:text>; </xsl:text>
          <strong>Project: </strong> <xsl:value-of select="aspect[@name='project_name']/value"/><xsl:text>; </xsl:text>
          <strong>Description: </strong> <xsl:value-of select="aspect[@name='description']/value"/><xsl:text> </xsl:text>
        </p>
      </div>
    </div> <!-- end search_result -->
  </xsl:template>

</xsl:stylesheet>
