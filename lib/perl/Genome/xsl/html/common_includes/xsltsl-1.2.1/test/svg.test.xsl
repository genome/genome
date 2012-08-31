<?xml version="1.0"?>
<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:svg="http://xsltsl.org/svg"
>
  <xsl:include href="../svg.xsl"/>

  <doc:article>
    <doc:title>SVG Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the SVG module.</doc:para>
  </doc:article>

  <xsl:template name="svg">

    <xsl:if test="$debug='true'">
      <xsl:message>SVG tests starting</xsl:message>
    </xsl:if>

    <xsl:if test="$debug='true'">
      <xsl:message>Test svg:aqua-button-defs template</xsl:message>
    </xsl:if>

    <!-- Don't actually do any testing - the proof is in the pudding! -->

  </xsl:template>

</xsl:stylesheet>
