<?xml version="1.0"?>
<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:markup="http://xsltsl.org/markup"
>
  <xsl:include href="../markup.xsl"/>

  <doc:article>
    <doc:title>Markup Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the Markup module.</doc:para>
  </doc:article>

  <xsl:template name="markup">

    <xsl:if test="$debug='true'"><xsl:message>Markup tests starting</xsl:message></xsl:if>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:xml-declaration template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:xml-declaration test (default value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:xml-declaration"/>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<?xml version="1.0"?>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:xml-declaration test (standalone value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:xml-declaration">
          <xsl:with-param name="standalone">yes</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<?xml version="1.0" standalone="yes"?>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:xml-declaration test (encoding value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:xml-declaration">
          <xsl:with-param name="encoding">utf-16</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<?xml version="1.0" encoding="utf-16"?>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:xml-declaration test (standalone and encoding values)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:xml-declaration">
          <xsl:with-param name="encoding">utf-8</xsl:with-param>
          <xsl:with-param name="standalone">no</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<?xml version="1.0" standalone="no" encoding="utf-8"?>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:doctype-declaration template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (doc element)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (external DTD subset, system identifier)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
          <xsl:with-param name="systemid">foo.dtd</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo SYSTEM "foo.dtd">
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (external DTD subset, public identifier)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
          <xsl:with-param name="publicid">-//Foo//Foo schema//EN</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo PUBLIC "-//Foo//Foo schema//EN">
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (external DTD subset, both identifiers)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
          <xsl:with-param name="systemid" select="'foobar.dtd'"/>
          <xsl:with-param name="publicid">-//Foo//Foo schema//EN</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo PUBLIC "-//Foo//Foo schema//EN" "foobar.dtd">
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (internal DTD subset)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
          <xsl:with-param name="internaldtd"><![CDATA[
<!ENTITY % foo SYSTEM "foo.mod">
%foo;
]]></xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo [
<!ENTITY % foo SYSTEM "foo.mod">
%foo;
]>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:doctype-declaration test (internal and external DTD subsets)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:doctype-declaration">
          <xsl:with-param name="docel">Foo</xsl:with-param>
          <xsl:with-param name="systemid" select="'foobar.dtd'"/>
          <xsl:with-param name="internaldtd"><![CDATA[
<!ENTITY % foo SYSTEM "foo.mod">
%foo;
]]></xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!DOCTYPE Foo SYSTEM "foobar.dtd" [
<!ENTITY % foo SYSTEM "foo.mod">
%foo;
]>
]]></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:element-declaration template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:element-declaration test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:element-declaration">
          <xsl:with-param name="type">Foo</xsl:with-param>
          <xsl:with-param name="content-spec" select="'EMPTY'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ELEMENT Foo EMPTY>]]></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:attlist-declaration template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:attlist-declaration test (explicit attr defns)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:attlist-declaration">
          <xsl:with-param name="type">Foo</xsl:with-param>
          <xsl:with-param name="attr-defns">
bar CDATA #IMPLIED</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ATTLIST Foo 
bar CDATA #IMPLIED>]]></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:attlist-declaration test (using attr defn)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:attlist-declaration">
          <xsl:with-param name="type">Foo</xsl:with-param>
          <xsl:with-param name="attr-defns">
            <xsl:call-template name="markup:attribute-definition">
              <xsl:with-param name="name" select="'bar'"/>
              <xsl:with-param name="type" select="'CDATA'"/>
              <xsl:with-param name="default" select="'#FIXED &quot;foobar&quot;'"/>
            </xsl:call-template>
          </xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ATTLIST Foo  bar CDATA #FIXED "foobar">]]></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:entity-declaration template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:entity-declaration test (internal general entity as text)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:entity-declaration">
          <xsl:with-param name="name">atest</xsl:with-param>
          <xsl:with-param name="text">this is a test</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ENTITY atest "this is a test">]]></xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:entity-declaration test (internal general entity with special characters)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:entity-declaration">
          <xsl:with-param name="name">atest</xsl:with-param>
          <xsl:with-param name="nodes" select="//module[@name='markup']/testdata/*"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ENTITY atest "&lt;sample>this is a test&lt;/sample>">]]></xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:entity-declaration test (external entity)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:entity-declaration">
          <xsl:with-param name="name">atest</xsl:with-param>
          <xsl:with-param name="systemid">external.xml</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<!ENTITY atest SYSTEM "external.xml">]]></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:quote-value template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:quote-value test (no value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:quote-value"/>
      </xsl:with-param>
      <xsl:with-param name="expect">""</xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:quote-value test (value has no quotes)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:quote-value">
          <xsl:with-param name="value">simple value</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">"simple value"</xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:quote-value test (value has double quotes)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:quote-value">
          <xsl:with-param name="value">double " quote</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">'double " quote'</xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:quote-value test (value has single quote)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:quote-value">
          <xsl:with-param name="value">single ' quote</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">"single ' quote"</xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:quote-value test (value has both quotes)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:quote-value">
          <xsl:with-param name="value">both "' quotes</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">"both &amp;quot;' quotes"</xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test markup:as-xml template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:as-xml test (no value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:as-xml"/>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:as-xml test (text)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:as-xml">
          <xsl:with-param name='nodes' select='//module[@name="markup"]/testdata/*/text()'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">this is a test</xsl:with-param>
    </xsl:call-template>
    <xsl:call-template name="test">
      <xsl:with-param name="description">markup:as-xml test (element)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="markup:as-xml">
          <xsl:with-param name='nodes' select='//module[@name="markup"]/testdata/*'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"><![CDATA[<sample>this is a test</sample>]]></xsl:with-param>
    </xsl:call-template>
    <!-- TODO: tests for attributes, comments, PIs -->

    <!-- TODO: tests for notation declarations, CDATA sections, entity references -->

  </xsl:template>

</xsl:stylesheet>
