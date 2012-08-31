<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version='1.0'
                xmlns:docbook="http://www.oasis-open.org/docbook/xml/4.0"
                xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
                xmlns:xt="http://www.jclark.com/xt"
                xmlns:saxon="http://icl.com/saxon"
                xmlns:lxslt="http://xml.apache.org/xslt"
                xmlns:xalanredirect="org.apache.xalan.xslt.extensions.Redirect"
                exclude-result-prefixes="doc docbook xsl"
                extension-element-prefixes="xt saxon xalanredirect lxslt">

<xsl:output
     method="xml"
  indent='yes'
/>

<xsl:preserve-space elements="xsl:variable"/>
<xsl:strip-space elements="xsl:stylesheet"/>

<!-- ********************************************************************
     $Id: xsl2jref.xsl,v 1.5 2004/10/08 07:25:03 balls Exp $
     ********************************************************************

     This file is part of the XSL DocBook Stylesheet distribution.
     See ../README or http://nwalsh.com/docbook/xsl/ for copyright
     and other information.

     ******************************************************************** -->

<!-- ==================================================================== -->

<xsl:param name="output-file" select="''"/>

<!-- ==================================================================== -->

<xsl:param name="includestylesheets" select="1"/>

<!-- ==================================================================== -->

<xsl:template match="/">
  <xsl:choose>
    <xsl:when test="$output-file = ''">
      <xsl:message terminate='yes'>
        <xsl:text>You must set the output-file parameter!</xsl:text>
      </xsl:message>
    </xsl:when>
    <xsl:when test="/xsl:stylesheet/doc:*">
      <xt:document method="xml" href="{$output-file}">
        <xsl:apply-templates/>
        <xsl:fallback>
          <xalanredirect:write file="{$output-file}">
            <xsl:apply-templates/>
            <xsl:fallback>
              <saxon:output method="xml" href="{$output-file}">
                <xsl:apply-templates/>
                <xsl:fallback>
                  <xsl:apply-templates/>
                </xsl:fallback>
              </saxon:output>
            </xsl:fallback>
          </xalanredirect:write>
        </xsl:fallback>
      </xt:document>
    </xsl:when>
    <xsl:otherwise>
      <!-- nop -->
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="*">
  <xsl:variable name="block-element" select="name(.) = 'para'
                                             or name(.) = 'sidebar'
                                             or name(.) = 'variablelist'
                                             or name(.) = 'itemizedlist'"/>
  <xsl:element namespace="" name="{name(.)}">
    <xsl:apply-templates select="@*" mode="copy-attr"/>
    <xsl:apply-templates/>
  </xsl:element>
</xsl:template>

<xsl:template match="@*" mode="copy-attr">
  <xsl:attribute name="{name(.)}">
    <xsl:value-of select="."/>
  </xsl:attribute>
</xsl:template>

<xsl:template match="xsl:*"></xsl:template>

<xsl:template match="xsl:stylesheet">
  <xsl:choose>
    <xsl:when test="doc:reference">
      <reference>
        <xsl:apply-templates/>
      </reference>
    </xsl:when>
    <xsl:when test="doc:book">
      <book>
        <xsl:apply-templates/>
      </book>
    </xsl:when>
    <xsl:when test="doc:chapter">
      <chapter>
        <xsl:apply-templates/>
      </chapter>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="xsl:include">
  <xsl:if test="$includestylesheets = '1'">
    <xsl:apply-templates select="document(@href)/*"/>
  </xsl:if>
</xsl:template>

<xsl:template match="@*" mode="copy">
  <xsl:text> </xsl:text>
  <xsl:value-of select="name(.)"/>
  <xsl:text>="</xsl:text>
  <xsl:value-of select="."/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="*" mode="copy">
  <xsl:variable name="content">
    <xsl:apply-templates mode="copy"/>
  </xsl:variable>

  <xsl:text>&lt;</xsl:text>
  <xsl:value-of select="name(.)"/>
  <xsl:apply-templates select="@*" mode="copy"/>
  <xsl:choose>
    <xsl:when test="$content = ''">
      <xsl:text>/</xsl:text>
      <xsl:text>&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&gt;</xsl:text>
      <xsl:copy-of select="$content"/>
      <xsl:text>&lt;/</xsl:text>
      <xsl:value-of select="name(.)"/>
      <xsl:text>&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="comment()" mode="copy">
  <xsl:text>&lt;!--</xsl:text>
  <xsl:value-of select="."/>
  <xsl:text>--&gt;</xsl:text>
</xsl:template>

<xsl:template match="processing-instruction()" mode="copy">
  <xsl:text>&lt;?</xsl:text>
  <xsl:value-of select="name(.)"/>
  <xsl:text> </xsl:text>
  <xsl:value-of select="."/>
  <xsl:text>?&gt;</xsl:text>
</xsl:template>

<!-- ==================================================================== -->

<xsl:template match="doc:reference|doc:book|doc:chapter">
  <!-- only process the children; doc:reference logically wraps the entire
       stylesheet, even if it can't syntactically do so. -->
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="doc:variable">
  <xsl:variable name="name" select="@name"/>

  <refentry id="var.{$name}">
    <refnamediv>
      <refname><xsl:value-of select="$name"/></refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <refsynopsisdiv>
      <synopsis>
        <xsl:apply-templates select="../xsl:variable[@name=$name]"
                             mode="copy-template"/>
      </synopsis>
    </refsynopsisdiv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<xsl:template match="doc:param">
  <xsl:variable name="name" select="@name"/>

  <refentry id="param.{$name}">
    <refnamediv>
      <refname><xsl:value-of select="$name"/></refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <refsynopsisdiv>
      <synopsis>
        <xsl:apply-templates select="../xsl:param[@name=$name]"
                             mode="copy-template"/>
      </synopsis>
    </refsynopsisdiv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<xsl:template match="doc:template[@name]">
  <xsl:variable name="name" select="@name"/>

  <refentry id="template.{$name}">
    <refnamediv>
      <refname><xsl:value-of select="$name"/></refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <refsynopsisdiv>
      <synopsis>
        <xsl:apply-templates select="../xsl:template[@name=$name]"
                             mode="copy-template"/>
      </synopsis>
    </refsynopsisdiv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<xsl:template match="doc:template[@match]">
  <xsl:variable name="match" select="@match"/>
  <xsl:variable name="mode" select="@mode"/>

  <refentry>
    <refnamediv>
      <refname>
        <xsl:value-of select="$match"/>
        <xsl:if test="@mode">
          <xsl:text> (in </xsl:text>
          <xsl:value-of select="$mode"/>
          <xsl:text> mode)</xsl:text>
        </xsl:if>
      </refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <refsynopsisdiv>
      <synopsis>
        <xsl:choose>
          <xsl:when test="@mode">
            <xsl:apply-templates select="../xsl:template[@match=$match and @mode=$mode]"
                                 mode="copy-template"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:apply-templates select="../xsl:template[@match=$match]"
                                 mode="copy-template"/>
          </xsl:otherwise>
        </xsl:choose>
      </synopsis>
    </refsynopsisdiv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<xsl:template match="doc:mode">
  <xsl:variable name="name" select="@mode"/>

  <refentry id="mode.{$name}">
    <refnamediv>
      <refname><xsl:value-of select="$name"/> mode</refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<xsl:template match="doc:attribute-set">
  <xsl:variable name="name" select="@mode"/>

  <refentry id="attrset.{$name}">
    <refnamediv>
      <refname><xsl:value-of select="$name"/> mode</refname>
      <xsl:apply-templates select="refpurpose"/>
    </refnamediv>
    <xsl:apply-templates select="*[name(.)!='refpurpose']"/>
  </refentry>
</xsl:template>

<!-- ==================================================================== -->

<xsl:template match="*|text()|comment()|processing-instruction()"
              mode="copy-template">
  <!-- suppress -->
</xsl:template>

<xsl:template match="xsl:param" mode="copy-template">
  <xsl:apply-templates select="." mode="copy"/>
</xsl:template>

<xsl:template match="xsl:template" mode="copy-template">
  <xsl:variable name="content">
    <xsl:apply-templates mode="copy-template"/>
  </xsl:variable>

  <xsl:text>&lt;</xsl:text>
  <xsl:value-of select="name(.)"/>
  <xsl:apply-templates select="@*" mode="copy"/>
  <xsl:choose>
    <xsl:when test="$content = ''">
      <xsl:text>/</xsl:text>
      <xsl:text>&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&gt;</xsl:text>
      <xsl:copy-of select="$content"/>
      <xsl:text>  ...</xsl:text>
      <xsl:text>&lt;/</xsl:text>
      <xsl:value-of select="name(.)"/>
      <xsl:text>&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- ==================================================================== -->

</xsl:stylesheet>
