<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_instrumentdata_flowcell_table">
    <xsl:comment>template: status/genome_instrumentdata.xsl; name="genome_model_link_table"</xsl:comment>
    <div class="generic_lister">
      <div class="box_header span-24 last rounded-top">
        <div class="box_title"><h3 class="genome_instrumentdata_flowcell_16 span-24 last">Flowcell Lanes</h3></div>
      </div>
      <div class="box_content rounded-bottom span-24 last">
        <table class="lister">
          <thead>
            <tr>
              <th>direction</th>
              <th>model</th>
              <th><br/></th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="aspect[@name='to_models']/object | aspect[@name='to_builds']/object">
              <xsl:call-template name="genome_model_link_table_row">
                <xsl:with-param name="type">to</xsl:with-param>
              </xsl:call-template>
            </xsl:for-each>
            <xsl:for-each select="aspect[@name='from_models']/object | aspect[@name='from_builds']/object">
              <xsl:call-template name="genome_model_link_table_row">
                <xsl:with-param name="type">from</xsl:with-param>
              </xsl:call-template>
            </xsl:for-each>
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>


  <xsl:template name="genome_instrumentdata_flowcell_row">
    <xsl:param name="type" select="''" />

    <xsl:comment>template: status/genome_instrumentdata.xsl; name="genome_model_link_table_table_row"</xsl:comment>

    <tr>
      <td><xsl:value-of select="$type"/></td>
      <td>
        <xsl:value-of select="@type"/><xsl:text>: </xsl:text>
        <xsl:choose>
          <xsl:when test="aspect[@name='name']/value">
            <xsl:value-of select="aspect[@name='name']/value"/> (#<xsl:value-of select="@id"/>)
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="@id"/>
          </xsl:otherwise>
        </xsl:choose>
      </td>
      <td class="buttons">
        <xsl:call-template name="object_link_button_tiny">
          <xsl:with-param name="icon" select="'sm-icon-extlink'" />
        </xsl:call-template>
      </td>
    </tr>
  </xsl:template>

</xsl:stylesheet>