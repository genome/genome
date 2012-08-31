<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_modelgroup" match="object[./types[./isa[@type='Genome::ModelGroup']]]">
    <div class="search_result">
      <div class="result_icon genome_modelgroup_32">
        <br/>
      </div>
      <div class="result">
        <h3>Model Group: <xsl:call-template name="object_link"/></h3>
        <p class="resource_buttons">
          <xsl:for-each select="aspect[@name='convergence_model']/object">
            <xsl:call-template name="object_link_button">
              <xsl:with-param name="linktext" select="'convergence model'" />
              <xsl:with-param name="icon" select="'sm-icon-extlink'" />
            </xsl:call-template>
            <xsl:text> </xsl:text>
            <xsl:choose>
              <xsl:when test="aspect[@name='last_complete_build']/object">
                <xsl:for-each select="aspect[@name='last_complete_build']/object">
                  <xsl:call-template name="object_link_button">
                    <xsl:with-param name="linktext" select="'last succeeded build'" />
                    <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                  </xsl:call-template>
                  <xsl:text> </xsl:text>
                  <xsl:variable name="build_directory_url">
                    <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='data_directory']/value)" />
                  </xsl:variable>
                  <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$build_directory_url'/><xsl:text>/reports/Summary/report.html</xsl:text></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>summary report</a>
                </xsl:for-each>
              </xsl:when>
              <xsl:otherwise>
                [No succeeded builds.]
              </xsl:otherwise>
            </xsl:choose>
          </xsl:for-each>
        </p>
        <p class="result_summary">
          <strong>Members: </strong>
          <xsl:for-each select="aspect[@name='models']/object">
            <xsl:value-of select="aspect[@name='name']/value"/>
            <xsl:if test="position() != last()"><xsl:text>, </xsl:text></xsl:if>
          </xsl:for-each>
        </p>
      </div>
    </div> <!-- end search_result -->
  </xsl:template>

</xsl:stylesheet>