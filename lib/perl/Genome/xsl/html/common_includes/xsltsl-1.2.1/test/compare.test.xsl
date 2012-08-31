<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:cmp="http://xsltsl.org/cmp"
  exclude-result-prefixes="doc cmp">

  <xsl:include href="../cmp.xsl"/>

  <doc:article>
    <doc:title>Comparison Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the comparison module.</doc:para>
  </doc:article>

  <xsl:template name="compare">

    <xsl:if test="$debug='true'">
      <xsl:message>Comparison tests starting</xsl:message>
    </xsl:if>

    <xsl:if test="$debug='true'">
      <xsl:message>Test cmp:diff template</xsl:message>
    </xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">cmp:diff test (empty nodesets)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="cmp:diff">
          <xsl:with-param name='ns1' select='/..'/>
          <xsl:with-param name='ns2' select='/..'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
