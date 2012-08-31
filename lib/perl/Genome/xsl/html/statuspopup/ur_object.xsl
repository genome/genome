<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rest="urn:rest">

  <xsl:template name="ur_object" match="object">
    <xsl:comment>
      name: ur_object  match: object
    </xsl:comment>
    <xsl:if test="count(aspect) > 0">
      <div class="aspects">
        <table class="aspects" cellpadding="0" cellspacing="0" border="0">
          <colgroup>
           <col width="40%"/>
           <col/>
          </colgroup>
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
      match: display_name
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
      match: exception
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
    <xsl:value-of select="."/>
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
