<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output
      method="xml"
      omit-xml-declaration="yes"
      doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"
      doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"
      indent="yes"/>

  <xsl:strip-space elements="*"/>


  <xsl:template match="/">
    <xsl:comment>template: /html/status/root.xsl; match="/" </xsl:comment>

    <xsl:variable name="title">
      <xsl:choose>
        <xsl:when test="/object">
          <xsl:value-of select="/object/label_name" />: <xsl:value-of select="/object/display_name" /> - Status
        </xsl:when>
        <xsl:otherwise>
          Status
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <html>
      <head>
        <xsl:call-template name="html_head_page">
          <xsl:with-param name="title" select="$title"/>
        </xsl:call-template>
      </head>

      <body>
        <div class="page">

          <xsl:apply-templates/>

        </div>
      </body>
    </html>

  </xsl:template>

</xsl:stylesheet>
