<?xml version="1.0"?>

<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:uri="http://xsltsl.org/uri"
>
  <xsl:include href="../uri.xsl"/>

  <doc:article>
    <doc:title>URI Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the uri stylesheet module.</doc:para>
  </doc:article>

  <xsl:template name="uri">

    <xsl:if test="$debug='true'"><xsl:message>URI tests starting</xsl:message></xsl:if>

    <xsl:if test="$debug='true'"><xsl:message>Test uri:is-absolute-uri template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:is-absolute-uri test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:is-absolute-uri">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'true'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:is-absolute-uri template (with relative URI)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:is-absolute-uri test (with relative URI)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:is-absolute-uri">
          <xsl:with-param name="uri" select="'path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="''"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-scheme template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-scheme test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-scheme">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'scheme'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-scheme template (with no scheme)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-scheme test (with no scheme)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-scheme">
          <xsl:with-param name="uri" select="'path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="''"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-authority template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-authority test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-authority">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'authority'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-authority template (with no scheme)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-authority test (with no scheme)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-authority">
          <xsl:with-param name="uri" select="'//authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'authority'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-authority template (with no path, query, or fragment)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-authority test (with no path, query, or fragment)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-authority">
          <xsl:with-param name="uri" select="'scheme://authority'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'authority'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'/path'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template (with path only)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test (with path only)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'path'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'path'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template (with no authority)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test (with no authority)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'scheme:/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'/path'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template (with no scheme or authority)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test (with no scheme or authority)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'path'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template (with no authority or '/')</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test (with no authority or '/')</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'scheme:path'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'path'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-path template (with no path)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-path test (with no path)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-path">
          <xsl:with-param name="uri" select="'scheme://authority?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="''"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-query template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-query test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-query">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'query'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-query template (with no query)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-query test (with no query)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-query">
          <xsl:with-param name="uri" select="'scheme://authority/path#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="''"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-query template (with no fragment)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-query test (with no fragment)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-query">
          <xsl:with-param name="uri" select="'scheme://authority/path?query'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'query'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-fragment template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-fragment test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-fragment">
          <xsl:with-param name="uri" select="'scheme://authority/path?query#fragment'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'fragment'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:get-uri-fragment template (with no fragment)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:get-uri-fragment test (with no fragment)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:get-uri-fragment">
          <xsl:with-param name="uri" select="'scheme://authority/path?query'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="''"/>
    </xsl:call-template>


    <!-- Normal Examples from RFC 2396 Appendix C.1. -->


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g:h)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g:h)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g:h'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'g:h'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (./g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (./g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'./g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g/)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g/)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g/'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (/g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (/g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'/g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (//g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (//g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'//g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (?y)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (?y)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'?y'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/?y'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g?y)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g?y)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g?y'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g?y'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (#s)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (#s)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'#s'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/d;p#s'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (#s with document URI)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (#s with document URI)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="reference" select="'#s'"/>
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="document" select="'http://foo/bar?baz#quux'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://foo/bar#s'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g#s)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g#s)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g#s'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g#s'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g?y#s)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g?y#s)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g?y#s'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g?y#s'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (;x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (;x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="';x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/;x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g;x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g;x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g;x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g;x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g;x?y#s)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g;x?y#s)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g;x?y#s'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g;x?y#s'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (.)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (.)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'.'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (./)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (./)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'./'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (..)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (..)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'..'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../..)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../..)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../..'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../../)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../../)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../../'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/g'"/>
    </xsl:call-template>


    <!-- Abnormal Examples from RFC 2396 Appendix C.2. -->


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template ("")</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test ("")</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="''"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/d;p'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template ("" with document URI)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test ("" with document URI)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="reference" select="''"/>
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="document" select="'http://foo/bar?baz#quux'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://foo/bar'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../../../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../../../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../../../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/../g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (../../../../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (../../../../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'../../../../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/../../g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (/./g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (/./g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'/./g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/./g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (/../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (/../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'/../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/../g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g.)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g.)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g.'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g.'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (.g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (.g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'.g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/.g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g..)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g..)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g..'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g..'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (..g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (..g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'..g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/..g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (./../g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (./../g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'./../g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/g'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (./g/.)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (./g/.)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'./g/.'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g/'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g/./h)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g/./h)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g/./h'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g/h'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g/../h)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g/../h)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g/../h'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/h'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g;x=1/./y)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g;x=1/./y)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g;x=1/./y'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g;x=1/y'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g;x=1/../y)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g;x=1/../y)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g;x=1/../y'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/y'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g?y/./x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g?y/./x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g?y/./x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g?y/./x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g?y/../x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g?y/../x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g?y/../x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g?y/../x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g#s/./x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g#s/./x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g#s/./x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g#s/./x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (g#s/../x)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (g#s/../x)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'g#s/../x'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http://a/b/c/g#s/../x'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test uri:resolve-uri template (http:g)</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">uri:resolve-uri test (http:g)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="uri:resolve-uri">
          <xsl:with-param name="base" select="'http://a/b/c/d;p?q'"/>
          <xsl:with-param name="reference" select="'http:g'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'http:g'"/>
    </xsl:call-template>


  </xsl:template>

</xsl:stylesheet>
