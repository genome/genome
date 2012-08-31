<?xml version="1.0"?>
<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:math="http://xsltsl.org/math"
>
  <xsl:include href="../math.xsl"/>

  <doc:article>
    <doc:title>Math Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the math module.</doc:para>
  </doc:article>

  <xsl:template name="math">

    <xsl:if test='$debug = "true"'><xsl:message>Math module tests starting</xsl:message></xsl:if>

    <xsl:if test='$debug = "true"'><xsl:message>Test math:power template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (exponent 0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='2'/>
          <xsl:with-param name='power' select='0'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">1</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (power 1)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='2'/>
          <xsl:with-param name='power' select='1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">2</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (arbitrary power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='3'/>
          <xsl:with-param name='power' select='4'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">81</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (negative power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='2'/>
          <xsl:with-param name='power' select='-1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">0.5</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (negative power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='5'/>
          <xsl:with-param name='power' select='-2'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">0.04</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (fractional power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='2'/>
          <xsl:with-param name='power' select='0.5'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (invalid number)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base'>wrong</xsl:with-param>
          <xsl:with-param name='power' select='2'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (invalid number, zero power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base'>wrong</xsl:with-param>
          <xsl:with-param name='power' select='0'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:power test (invalid power)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:power">
          <xsl:with-param name='base' select='3'/>
          <xsl:with-param name='power' select='"wrong"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:if test='$debug = "true"'><xsl:message>Test math:power template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:cvt-hex-decimal test (empty value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:cvt-hex-decimal">
          <xsl:with-param name='value'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:cvt-hex-decimal test (invalid value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:cvt-hex-decimal">
          <xsl:with-param name='value' select='"fooey"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select='"NaN"'/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:cvt-hex-decimal test (single digit value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:cvt-hex-decimal">
          <xsl:with-param name='value' select='"3"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select='3'/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:cvt-hex-decimal test (high single digit value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:cvt-hex-decimal">
          <xsl:with-param name='value' select='"c"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select='12'/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:cvt-hex-decimal test (two digit value)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:cvt-hex-decimal">
          <xsl:with-param name='value' select='"ff"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select='255'/>
    </xsl:call-template>

    <xsl:if test='$debug = "true"'><xsl:message>Test math:ordinal</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='0'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">0th</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (less than 0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='-1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (1)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">1st</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (2)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='2'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">2nd</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (3)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='3'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">3rd</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (4)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='4'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">4th</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (less than 10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='8'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">8th</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='10'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">10th</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (less than 20)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='12'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">12th</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (21)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='21'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">21st</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal test (1002)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal">
          <xsl:with-param name='number' select='1002'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">1002nd</xsl:with-param>
    </xsl:call-template>

    <xsl:if test='$debug = "true"'><xsl:message>Test math:ordinal-as-word</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='0'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">zeroth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (less than 0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='-1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (1)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">first</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (2)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='2'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">second</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (3)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='3'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">third</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (4)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='4'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">fourth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (less than 10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='8'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">eighth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='10'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">tenth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (less than 20)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='12'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twelveth</xsl:with-param>
    </xsl:call-template>


    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (teens)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='15'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">fifteenth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (21)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='21'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twenty first</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (hundreds)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='334'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">three hundred and thirty fourth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (thousands - 1002)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1002'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one thousand and second</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (thousands - 1409)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1409'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one thousand four hundred and ninth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (thousands - 4957)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='4957'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four thousand nine hundred and fifty seventh</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (tens of thousands - 24817)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='24817'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twenty four thousand eight hundred and seventeenth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (hundreds of thousands - 481700)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='481700'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four hundred and eighty one thousand seven hundredth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (hundreds of thousands - 481703)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='481703'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four hundred and eighty one thousand seven hundred and third</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1000000'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one millionth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1000100)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1000100'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million and one hundredth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1000004)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1000004'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million and fourth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1000013)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1000013'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million and thirteenth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1000990)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1000990'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million nine hundred and ninetieth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1003000)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1003000'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million and three thousandth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1045768)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1045768'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million forty five thousand seven hundred and sixty eighth</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 1240102)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='1240102'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million two hundred and forty thousand one hundred and second</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:ordinal-as-word test (millions - 20240102)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:ordinal-as-word">
          <xsl:with-param name='number' select='20240102'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twenty million two hundred and forty thousand one hundred and second</xsl:with-param>
    </xsl:call-template>

    <xsl:if test='$debug = "true"'><xsl:message>Test math:number-as-word</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='0'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">zero</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (less than 0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='-1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">minus one</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (1)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (2)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='2'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">two</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (3)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='3'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">three</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (4)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='4'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (less than 10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='8'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">eight</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (10)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='10'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">ten</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (less than 20)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='12'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twelve</xsl:with-param>
    </xsl:call-template>


    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (teens)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='15'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">fifteen</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (21)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='21'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twenty one</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (hundreds)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='334'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">three hundred and thirty four</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (thousands)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1002'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one thousand and two</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (thousands)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1409'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one thousand four hundred and nine</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (thousands)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='4957'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four thousand nine hundred and fifty seven</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (tens of thousands)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='24817'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">twenty four thousand eight hundred and seventeen</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (hundreds of thousands - 481700)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='481700'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four hundred and eighty one thousand seven hundred</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (hundreds of thousands - 481705)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='481705'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">four hundred and eighty one thousand seven hundred and five</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (millions)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1000000'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (millions)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1000100'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million one hundred</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">math:number-as-word test (millions)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="math:number-as-word">
          <xsl:with-param name='number' select='1240102'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">one million two hundred and forty thousand one hundred and two</xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
