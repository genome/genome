<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="static_html_holder" match="object[@type='Genome::View::Static::Html']">
    <xsl:call-template name="control_bar_view"/>

    <div class="content rounded shadow" style="padding-top: 0;">
      <div class="header rounded-top gradient-grey">
        <div class="container">
          <div><xsl:attribute name="class">title span-24 last <xsl:value-of select="'app_analysis_search_32'"/></xsl:attribute>
          <h1>
            <xsl:value-of select="./display_name"/>
          </h1>
          </div>
        </div>
      </div>
      <div class="container">
        <xsl:copy-of select='./body_html'/>
      </div>
    </div>

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>
  </xsl:template>

</xsl:stylesheet>
