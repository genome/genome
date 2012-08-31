<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:strip-space elements="*"/>

  <xsl:template match="/">
    <xsl:comment>template: /html/statuspopup/root.xsl match="/"</xsl:comment>

    <xsl:apply-templates/>
  </xsl:template>

</xsl:stylesheet>
