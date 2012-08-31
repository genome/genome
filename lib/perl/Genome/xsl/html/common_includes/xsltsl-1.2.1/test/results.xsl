<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
>
  <xsl:template match="results">
    <html>
      <head>
        <title>XSLTSL Test Results</title>
        <link rel="stylesheet" type="text/css" href="results.css" />
      </head>
      <body>
        <p>
          <b>Modules: </b>
          <xsl:for-each select="module">
            <a href="#{@name}"><xsl:value-of select="@name" /></a>
            <xsl:if test="position() &lt; last()">, </xsl:if>
          </xsl:for-each>
          <br />
          <b>Summary: </b>
          <xsl:value-of select="count(module/test[@pass = 'true'])" />
          <xsl:text> tests passed. </xsl:text>
          <xsl:value-of select="count(module/test[@pass = 'false'])" />
          <xsl:text> tests failed.</xsl:text>
        </p>
        <xsl:if test='@processor'>
          <p>
            <xsl:text>Tests run using XSLT processor "</xsl:text>
            <xsl:value-of select='@processor'/>
            <xsl:text>".</xsl:text>
          </p>
        </xsl:if>
        <xsl:apply-templates />
      </body>
    </html>
  </xsl:template>
  
  <xsl:template match="module">
    <a name="{@name}" />
    <table class="module">
      <tr>
        <th colspan="3">
          <xsl:text>Results for </xsl:text>
          <font size="+1"><xsl:value-of select="@name" /></font>
          <xsl:text> module</xsl:text>
        </th>
      </tr>
      <xsl:apply-templates />
    </table>
  </xsl:template>
  
  <xsl:template match="test">
    <tr>
      <xsl:attribute name="class">
        <xsl:choose>
          <xsl:when test="@pass = 'true'">pass</xsl:when>
          <xsl:otherwise>fail</xsl:otherwise>
        </xsl:choose>
      </xsl:attribute>
      <td><xsl:value-of select="description" /></td>
      <td><xsl:value-of select="expected" /></td>
      <td><xsl:value-of select="result" /></td>
    </tr>
  </xsl:template>
  
</xsl:stylesheet>
