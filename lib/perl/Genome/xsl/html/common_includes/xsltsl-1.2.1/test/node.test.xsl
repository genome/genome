<?xml version="1.0"?>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:node="http://xsltsl.org/node"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
>
  <xsl:include href="../node.xsl"/>

  <doc:article>
    <doc:title>Node Module Test Suite</doc:title>
    <doc:para>Tests for the node stylesheet module.</doc:para>
  </doc:article>

  <xsl:template name="node">

    <xsl:if test="$debug='true'">
      <xsl:message>Node tests starting</xsl:message>
    </xsl:if>

    <xsl:if test="$debug='true'">
      <xsl:message>Test node:xpath template</xsl:message>
    </xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">node:xpath test (root node)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="node:xpath">
	  <xsl:with-param name="node" select="/"/>
	</xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">/</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">node:xpath test (element node)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="node:xpath">
	  <xsl:with-param name="node" select="/*"/>
	</xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">/tests[1]</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">node:xpath test (element node, three levels deep)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="node:xpath">
          <xsl:with-param name="node" select="/tests/module[@name='node']/testdata/element[1]"/>
	</xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">/tests[1]/module[3]/testdata[1]/element[1]</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">node:xpath test (element node, multiple siblings)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="node:xpath">
	  <xsl:with-param name="node" select="/*/module[2]"/>
	</xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">/tests[1]/module[2]</xsl:with-param>
    </xsl:call-template>

<!--
    <xsl:call-template name="test">
      <xsl:with-param name="description">node:xpath test (text node)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="node:xpath">
	  <xsl:with-param name="node" select="//module[@name = 'node']//text()[1]"/>
	</xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">/tests/module[3]/testdata/element[1]/text/text()</xsl:with-param>
    </xsl:call-template>
-->

  </xsl:template>

</xsl:stylesheet>
