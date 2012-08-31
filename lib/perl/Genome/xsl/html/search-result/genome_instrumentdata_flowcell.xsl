<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_instrumentdata_flowcell" match="object[./types[./isa[@type='Genome::InstrumentData::FlowCell']]]">
    <div class="search_result">
      <div class="result_icon genome_instrumentdata_flowcell_32">
        <br/>
      </div>
      <div class="result">
        <h3>Flow Cell:
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext">
            <xsl:value-of select="aspect[@name='flow_cell_id']/value"/>
          </xsl:with-param>
        </xsl:call-template>
        </h3>
        <p class="resource_buttons">
          <a class="mini btn"><xsl:attribute name="href">/solexa/equipment/flowcell?flow_cell_id=<xsl:value-of select="normalize-space(aspect[@name='flow_cell_id']/value)"/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>production reports</a>
          <xsl:choose>
            <xsl:when test="aspect[@name='lanes']/object">
              <xsl:call-template name="object_link_button">
                <xsl:with-param name="linktext" select="'analysis reports'"/>
                <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
              </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
              <span style="color: #aaa;">(no analysis reports available)</span>
            </xsl:otherwise>
          </xsl:choose>

        </p>
        <p class="result_summary">
          <strong>Machine Name: </strong><xsl:value-of select="aspect[@name='machine_name']/value"/><xsl:text>; </xsl:text>
          <strong>Run Type: </strong><xsl:value-of select="aspect[@name='run_name']/value"/> <xsl:value-of select="aspect[@name='run_type']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
