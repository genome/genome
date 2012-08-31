<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rest="urn:rest">

  <xsl:strip-space elements="*"/>

  <xsl:template name="object_link_href">
    <xsl:param name="rest" select="$rest" />
    <xsl:param name="type" select="./@type"/>
    <xsl:param name="key" select="'id'"/>
    <xsl:param name="keys" select="''"/>
    <xsl:param name="id" select="./@id"/>
    <xsl:param name="perspective" select="'status'"/>
    <xsl:param name="toolkit" select="'html'"/>

    <xsl:value-of select="$rest"/>
    <xsl:text>/</xsl:text>
    <xsl:value-of select="rest:typetourl($type)"/>
    <xsl:text>/</xsl:text>
    <xsl:value-of select="$perspective"/>
    <xsl:text>.</xsl:text>
    <xsl:value-of select="$toolkit"/>
    <xsl:text>?</xsl:text>

    <!-- Can pass in a xml node with children having key/value pairs for URL generation.
    Each child with a key attribute will be used.
    The url will be generated with the key value as key and content of the tag as value
    Example:
        <keys>
          <pair1 key='status'>Completed</key>
          <pair2 key='user'>guest</key>
        </keys>
    Turns into:
        ?status=Completed&user=guest
    -->
    <xsl:choose>
      <xsl:when test='$keys'>
        <xsl:for-each select="$keys/*">
          <xsl:choose>
            <xsl:when test="@key">
              <xsl:value-of select="@key"/>
              <xsl:text>=</xsl:text>
              <xsl:value-of select="."/>
              <xsl:text>&amp;</xsl:text>
            </xsl:when>
          </xsl:choose>
        </xsl:for-each>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$key"/>
        <xsl:text>=</xsl:text>
        <xsl:value-of select="$id"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="object_link">
    <xsl:param name="type" select="./@type"/>
    <xsl:param name="id" select="./@id"/>
    <xsl:param name="perspective" select="'status'"/>
    <xsl:param name="toolkit" select="'html'"/>
    <xsl:param name="linktext" select="./aspect[@name='name']/value"/>

    <a>
      <xsl:attribute name="href">
        <xsl:call-template name="object_link_href">
          <xsl:with-param name="type" select="$type"/>
          <xsl:with-param name="id" select="$id"/>
          <xsl:with-param name="perspective" select="$perspective"/>
          <xsl:with-param name="toolkit" select="$toolkit"/>
        </xsl:call-template>
      </xsl:attribute>
      <xsl:copy-of select="$linktext"/>
    </a>
  </xsl:template>

  <xsl:template name="string-replace-all">
    <xsl:param name="text" />
    <xsl:param name="replace" />
    <xsl:param name="by" />
    <xsl:choose>
      <xsl:when test="contains($text, $replace)">
        <xsl:value-of select="substring-before($text,$replace)" />
        <xsl:value-of select="$by" />
        <xsl:call-template name="string-replace-all">
          <xsl:with-param name="text"
                          select="substring-after($text,$replace)" />
          <xsl:with-param name="replace" select="$replace" />
          <xsl:with-param name="by" select="$by" />
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$text" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- creates a button with a jQueryUI icon -->
  <xsl:template name="object_link_button">
    <xsl:param name="type" select="./@type"/>
    <xsl:param name="id" select="./@id"/>
    <xsl:param name="key" select="'id'"/>
    <xsl:param name="keys" select="''"/>
    <xsl:param name="perspective" select="'status'"/>
    <xsl:param name="toolkit" select="'html'"/>
    <xsl:param name="linktext" select="./aspect[@name='name']/value"/>
    <xsl:param name="icon"/>

    <xsl:comment>template: common_includes/links.xsl; name: object_link_button</xsl:comment>

    <a class="mini btn">
      <xsl:attribute name="href">
        <xsl:call-template name="object_link_href">
          <xsl:with-param name="type" select="$type"/>
          <xsl:with-param name="key" select="$key"/>
          <xsl:with-param name="keys" select="$keys"/>
          <xsl:with-param name="id" select="$id"/>
          <xsl:with-param name="perspective" select="$perspective"/>
          <xsl:with-param name="toolkit" select="$toolkit"/>
        </xsl:call-template>
      </xsl:attribute>
      <xsl:value-of select="$linktext"/>
    </a>

  </xsl:template>

  <!-- creates a tiny button with no label -->
  <xsl:template name="object_link_button_tiny">
    <xsl:param name="type" select="./@type"/>
    <xsl:param name="id" select="./@id"/>
    <xsl:param name="perspective" select="'status'"/>
    <xsl:param name="toolkit" select="'html'"/>
    <xsl:param name="icon"/>

    <xsl:variable name="button_href">
      <xsl:call-template name="object_link_href">
        <xsl:with-param name="type" select="$type"/>
        <xsl:with-param name="id" select="$id"/>
        <xsl:with-param name="perspective" select="$perspective"/>
        <xsl:with-param name="toolkit" select="$toolkit"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:comment>template: status/ur_object.xsl:object_link_button_tiny</xsl:comment>

    <a class="mini-icon btn">
      <xsl:attribute name="href">
        <xsl:value-of select="$button_href"/>
      </xsl:attribute>
      <span><xsl:attribute name="class"><xsl:text>sm-icon </xsl:text><xsl:value-of select="$icon"/></xsl:attribute><br/></span>
    </a>

  </xsl:template>

</xsl:stylesheet>
