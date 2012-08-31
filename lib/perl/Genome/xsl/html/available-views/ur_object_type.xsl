<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rest="urn:rest">

  <xsl:template name="ur_object_type" match="perspectives">
    <xsl:comment>
      name: ur_object_type  match: perspectives
    </xsl:comment>
    <h2>Available Perspectives</h2>
    <ul>
      <xsl:for-each select="perspective[./toolkit/@name='html']">
        <xsl:if test="@name != 'static'">
          <li>
            <a>
              <xsl:attribute name="href">
                <xsl:call-template name="string-replace-all">
                  <xsl:with-param name="text">
                    <xsl:call-template name="object_link_href">
                       <xsl:with-param name="type" select="/perspectives/@type" />
                       <xsl:with-param name="perspective" select="./@name" />
                       <xsl:with-param name="id" select="''" />
                    </xsl:call-template>
                  </xsl:with-param>
                  <xsl:with-param name="replace" select="'?id='" />
                  <xsl:with-param name="by" select="''" />
                </xsl:call-template>
              </xsl:attribute>
              <xsl:value-of select="./@name" />
            </a>
          </li>
        </xsl:if>
      </xsl:for-each>
    </ul>
  </xsl:template>
</xsl:stylesheet>
