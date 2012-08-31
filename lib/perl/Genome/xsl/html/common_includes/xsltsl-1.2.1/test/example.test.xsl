<?xml version="1.0"?>
<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:eg="http://xsltsl.org/example"
>
  <xsl:include href="../example.xsl"/>

  <doc:article>
    <doc:title>Example Module Test Suite</doc:title>
    <doc:para>This stylesheet is a template for a module test.</doc:para>
  </doc:article>

  <xsl:template name="example">

    <xsl:if test="$debug='true'">
      <xsl:message>Example tests starting</xsl:message>
    </xsl:if>

    <xsl:if test="$debug='true'">
      <xsl:message>Test eg:example template</xsl:message>
    </xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">eg:example test (default value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="eg:example"/>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
