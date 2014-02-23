<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_model" match="object[./types[./isa[@type='Genome::Model']]]">
    <!-- These parameters are only valid when there is a last_complete_build; they are here for overriding by subclasses -->
    <xsl:param name="build_directory_url">
      <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='last_complete_build']/object/aspect[@name='data_directory']/value)" />
    </xsl:param>
    <xsl:param name="summary_report_url">
      <xsl:value-of select="$build_directory_url"/><xsl:text>/reports/Summary/report.html</xsl:text>
    </xsl:param>

    <div class="search_result">
      <div class="result_icon genome_model_32">
        <br/>
      </div>
      <div class="result">
        <h3>Model: <xsl:call-template name="object_link" /></h3>
        <p class="resource_buttons">
          <xsl:choose>
            <xsl:when test="aspect[@name='last_complete_build']/object">
              <xsl:for-each select="aspect[@name='last_complete_build']/object">
                <xsl:call-template name="object_link_button">
                  <xsl:with-param name="linktext" select="'last succeeded build'" />
                  <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                </xsl:call-template>

                <xsl:text> </xsl:text>

                <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$build_directory_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>data directory</a>

                <xsl:text> </xsl:text>

                <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$summary_report_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>summary report</a>
              </xsl:for-each>
            </xsl:when>
            <xsl:otherwise>
              This model contains no succeeded builds.
            </xsl:otherwise>
          </xsl:choose>
        </p>
        <p class="result_summary">
          <strong>ID: </strong><xsl:value-of select="aspect[@name='genome_model_id']/value"/>
          <xsl:text>; </xsl:text>

          <xsl:variable name="model_type">
                <xsl:call-template name="substring-after-last">
                <xsl:with-param name="input" select="@type"/>
                <xsl:with-param name="substr" select="':'"/>
                </xsl:call-template>
          </xsl:variable>

          <strong>Created by: </strong><xsl:value-of select="aspect[@name='created_by']/value"/>
          <xsl:text>; </xsl:text>
          <strong>Run as: </strong><xsl:value-of select="aspect[@name='run_as']/value"/>
          <xsl:text>; </xsl:text>
          <strong>Type: </strong><xsl:value-of select="$model_type"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
