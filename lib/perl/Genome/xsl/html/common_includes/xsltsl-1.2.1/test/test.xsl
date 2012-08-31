<?xml version="1.0"?>
<axsl:stylesheet xmlns:doc="http://xsltsl.org/xsl/documentation/1.0" xmlns:axsl="http://www.w3.org/1999/XSL/Transform" exclude-result-prefixes="doc" version="1.0">
  <axsl:param name="debug" select="'false'"/>
<!--====================================-->
<!--=                                  =-->
<!--=   DO NOT EDIT THIS STYLESHEET    =-->
<!--=                                  =-->
<!--= This stylesheet is generated     =-->
<!--= using the gentest.xsl stylesheet =-->
<!--=                                  =-->
<!--====================================-->
  <doc:chapter xmlns="">
    <chapterinfo xmlns="">
      <author>
        <surname>Ball</surname>
        <firstname>Steve</firstname>
      </author>
      <copyright>
        <year>2001</year>
        <holder>Steve Ball</holder>
      </copyright>
    </chapterinfo>
    <title xmlns="">XSLT Standard Library Test Suite</title>
    <section xmlns="">
      <title>About The Test Suite</title>
      <para>This stylesheet tests the stylesheet modules contained in the XSLT Standard Library</para>
    </section>
    <section xmlns="">
      <title>Preparing The Tests</title>
      <para>The tests are performed using the stylesheet
          <filename>test.xsl</filename>.  Normally this stylesheet is supplied
          with the XSLTSL distribution, so you shouldn't need to do
          anything.</para>
      <para>However, if you have added a new stylesheet module to
          the library then the <filename>test.xsl</filename> stylesheet must be
          regenerated.  This is done using the <filename>gentest.xsl</filename>
          stylesheet.  To do this, process the <filename>test.xml</filename>
          file with the stylesheet <filename>gentest.xsl</filename>
          using your favourite XSLT processor.</para>
    </section>
    <section xmlns="">
      <title>Running The Tests</title>
      <para>Process the <filename>test.xml</filename> file with
          the stylesheet <filename>test.xsl</filename> using your
          favourite XSLT processor.  The result is an XML document.
          You can process this with
          the <filename>results.xsl</filename> stylesheet to get formatted html
          output. A modern browser combined with the
          <filename>results.css</filename> file gives you green lines for the
          tests that passed and ugly red ones for the ones that didn't. When the
          tests fail they also print a text message to the console for
          quick verification of the correctness.</para>
    </section>
    <section xmlns="">
      <title>Debugging mode</title>
      <para>When a test seems to be stuck in a loop or if you want
          to verify the tests in some other way, start the
          <filename>test.xsl</filename> stylesheet with the parameter
          <userinput>debug='true'</userinput> to get more info printed
          to the console.
        </para>
    </section>
  </doc:chapter>
  <axsl:include href="string.test.xsl"/>
  <axsl:include href="date-time.test.xsl"/>
  <axsl:include href="node.test.xsl"/>
  <axsl:include href="uri.test.xsl"/>
  <axsl:include href="compare.test.xsl"/>
  <axsl:include href="markup.test.xsl"/>
  <axsl:include href="math.test.xsl"/>
  <axsl:include href="svg.test.xsl"/>
  <axsl:output method="xml" encoding="utf-8"/>
  <axsl:template match="tests">
    <axsl:processing-instruction name="xml-stylesheet">type="text/xsl" href="results.xsl"</axsl:processing-instruction>
    <results>
      <module name="string">
        <axsl:call-template name="string"/>
      </module>
      <module name="date-time">
        <axsl:call-template name="date-time"/>
      </module>
      <module name="node">
        <axsl:call-template name="node"/>
      </module>
      <module name="uri">
        <axsl:call-template name="uri"/>
      </module>
      <module name="compare">
        <axsl:call-template name="compare"/>
      </module>
      <module name="markup">
        <axsl:call-template name="markup"/>
      </module>
      <module name="math">
        <axsl:call-template name="math"/>
      </module>
      <module name="svg">
        <axsl:call-template name="svg"/>
      </module>
    </results>
  </axsl:template>
  <axsl:template name="test">
    <axsl:param name="expect"/>
    <axsl:param name="result"/>
    <axsl:param name="description">(no description)</axsl:param>
    <axsl:param name="passed">
      <axsl:value-of select="$expect = $result"/>
    </axsl:param>
    <axsl:if test="$passed='false'">
      <axsl:message>FAILED: <axsl:value-of select="$description"/></axsl:message>
    </axsl:if>
    <test>
      <axsl:attribute name="pass">
        <axsl:value-of select="$passed"/>
      </axsl:attribute>
      <description>
        <axsl:value-of select="$description"/>
      </description>
      <expected>
        <axsl:value-of select="$expect"/>
      </expected>
      <result>
        <axsl:value-of select="$result"/>
      </result>
    </test>
  </axsl:template>
</axsl:stylesheet>
