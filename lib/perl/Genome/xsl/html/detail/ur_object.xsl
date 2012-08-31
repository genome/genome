<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rest="urn:rest">
  <xsl:template name="ur_object_header" match="/object">
    <xsl:comment>template: status/ur_object.xsl:ur_object_header  match: /object</xsl:comment>
    <div class="object">
      <div class="header_object">
        <div class="display_name">
          <xsl:call-template name="header_display_name"/>
        </div>
      </div>
      <div class="content">
        <xsl:if test="count(aspect) > 0">
          <div class="aspects">
            <table class="aspects" cellpadding="0" cellspacing="0" border="0">
              <tbody>
                <xsl:for-each select="aspect">
                  <tr>
                    <td class="name">
                      <strong><xsl:value-of select="@name"/></strong>
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="normalize-space(.)">
                          <xsl:apply-templates/>
                        </xsl:when>
                        <xsl:otherwise>
                          <p>--</p>
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>
                </xsl:for-each>
              </tbody>
            </table>
          </div>
        </xsl:if>
      </div>
    </div>

  </xsl:template>

  <xsl:template name="header_display_name">
    <xsl:comment>template: status/ur_object.xsl:header_display_name match: header_display_name</xsl:comment>
    <span style="font-weight: bold;">
      <h1>
        <xsl:value-of select="@type"/><span class="id"> (<xsl:value-of select="display_name"/>)</span>
      </h1>
    </span>
  </xsl:template>

  <xsl:template name="ur_object" match="object">
    <xsl:comment>template: status/ur_object.xsl:ur_object match: object</xsl:comment>
    <p>
      <xsl:apply-templates select="display_name"/>
    </p>
    <xsl:if test="count(aspect) > 0">
      <div class="aspects">
        <table class="aspects" cellpadding="0" cellspacing="0" border="0">
          <tbody>
            <xsl:for-each select="aspect">
              <tr>
                <td class="name">
                  <strong><xsl:value-of select="@name"/></strong>
                </td>
                <td class="value"><xsl:apply-templates/></td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>
      </div>
    </xsl:if>
  </xsl:template>

  <xsl:template match="display_name">
    <xsl:comment>
      template: status/ur_object.xsl:display_name match: display_name
    </xsl:comment>
    <xsl:variable name="typeLink">
      <xsl:value-of select="rest:typetourl(../@type)" />
    </xsl:variable>
    <span>
      <span class="display_name"><xsl:value-of select="../@type"/></span><span class="id"> (<a>
      <xsl:attribute name="href">
        <xsl:value-of select="$rest"/>
        <xsl:text>/</xsl:text>
        <xsl:value-of select="$typeLink"/>
        <xsl:text>/</xsl:text>
        <xsl:value-of select="$currentPerspective"/>
        <xsl:text>.</xsl:text>
        <xsl:value-of select="$currentToolkit"/>
        <xsl:text>?id=</xsl:text>
        <xsl:value-of select="../@id"/>
      </xsl:attribute>
      <xsl:value-of select="."/>
      </a>)</span>
    </span>
  </xsl:template>

  <xsl:template match="exception">
    <xsl:comment>
      template: status/ur_object.xsl:exception match: exception
    </xsl:comment>
    <p class="exception">
      Exception <span class="trigger">[toggle view]</span>
    </p>
    <div class="toggle_container">
      <p><xsl:value-of select="."/></p>
    </div>
  </xsl:template>

  <xsl:template match="value">
    <xsl:comment>
      match: value
    </xsl:comment>
    <p><xsl:value-of select="."/></p>
  </xsl:template>

  <xsl:template match="perldata/scalar">
    <xsl:comment>
      match: perldata/scalar
    </xsl:comment>
    <p><xsl:value-of select="."/></p>
  </xsl:template>

  <xsl:template match="perldata/scalarref">
    <xsl:comment>
      match: perldata/scalarref
    </xsl:comment>
    <p><xsl:value-of select="."/></p>
  </xsl:template>

  <xsl:template match="perldata/arrayref">
    <xsl:comment>
      match: perldata/arrayref
    </xsl:comment>
    <p><xsl:value-of select="@blessed_package"/>=ARRAY(<xsl:value-of select="@memory_address"/>)</p>
  </xsl:template>

  <xsl:template match="perldata/hashref">
    <xsl:comment>
      match: perldata/hashref
    </xsl:comment>
    <p><xsl:value-of select="@blessed_package"/>=HASH(<xsl:value-of select="@memory_address"/>) <span class="trigger">[toggle view]</span></p>
    <div class="toggle_container">
      <table class="hash">
        <tbody>
          <xsl:for-each select="item">
            <tr>
              <td class="name"><xsl:value-of select="@key"/></td>
              <td class="value"><xsl:value-of select="."/></td>
            </tr>
          </xsl:for-each>
        </tbody>
      </table>
    </div>
  </xsl:template>

</xsl:stylesheet>
