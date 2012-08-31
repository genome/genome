<?xml version="1.0"?>

<xsl:stylesheet
  version="1.0"
  exclude-result-prefixes="doc"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:doc="http://xsltsl.org/xsl/documentation/1.0"
  xmlns:dt="http://xsltsl.org/date-time"
>
  <xsl:include href="../date-time.xsl"/>

  <doc:article>
    <doc:title>Date/Time Module Test Suite</doc:title>
    <doc:para>This stylesheet tests the date-time stylesheet module.</doc:para>
  </doc:article>

  <xsl:template name="date-time">

    <xsl:if test="$debug='true'"><xsl:message>Date/Time tests starting</xsl:message></xsl:if>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time test (with default format)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="year" select="2001"/>
          <xsl:with-param name="month" select="3"/>
          <xsl:with-param name="day" select="29"/>
          <xsl:with-param name="hour" select="19"/>
          <xsl:with-param name="minute" select="26"/>
          <xsl:with-param name="second" select="17"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">2001-03-29T19:26:17</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time test (xsd datetime, with default format)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name='xsd-date-time'>2001-03-29T19:26:17</xsl:with-param>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'2001-03-29T19:26:17'"/>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="year" select="1999"/>
          <xsl:with-param name="month" select="4"/>
          <xsl:with-param name="day" select="1"/>
          <xsl:with-param name="hour" select="12"/>
          <xsl:with-param name="minute" select="0"/>
          <xsl:with-param name="second" select="0"/>
          <xsl:with-param name="format" select="'%A, %b %d, %Y'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'Thursday, Apr 01, 1999'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="year" select="1"/>
          <xsl:with-param name="month" select="2"/>
          <xsl:with-param name="day" select="3"/>
          <xsl:with-param name="hour" select="4"/>
          <xsl:with-param name="minute" select="5"/>
          <xsl:with-param name="second" select="6"/>
          <xsl:with-param name="format" select="'%a %A %b %B %d %H %I %m %M %p %S %w %Y %%'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'Sat Saturday Feb February 03 04 04 02 05 am 06 6 0001 %'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%I, hour=0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="0"/>
          <xsl:with-param name="format" select="'%I'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'12'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%I, hour=9)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="9"/>
          <xsl:with-param name="format" select="'%I'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'09'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%I, hour=12)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="12"/>
          <xsl:with-param name="format" select="'%I'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'12'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%I, hour=21)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="21"/>
          <xsl:with-param name="format" select="'%I'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'09'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%I, hour=23)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="23"/>
          <xsl:with-param name="format" select="'%I'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'11'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%i, hour=0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="0"/>
          <xsl:with-param name="format" select="'%i'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'12'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%i, hour=9)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="9"/>
          <xsl:with-param name="format" select="'%i'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'9'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%i, hour=12)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="12"/>
          <xsl:with-param name="format" select="'%i'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'12'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%i, hour=21)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="21"/>
          <xsl:with-param name="format" select="'%i'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'9'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%i, hour=23)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="23"/>
          <xsl:with-param name="format" select="'%i'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'11'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%p, hour=0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="0"/>
          <xsl:with-param name="format" select="'%p'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'am'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%p, hour=12)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="12"/>
          <xsl:with-param name="format" select="'%p'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'pm'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%P, hour=0)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="0"/>
          <xsl:with-param name="format" select="'%P'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'AM'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time (%P, hour=12)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="hour" select="12"/>
          <xsl:with-param name="format" select="'%P'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'PM'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-day-of-the-week template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-day-of-the-week test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-day-of-the-week">
          <xsl:with-param name="year" select="1974"/>
          <xsl:with-param name="month" select="5"/>
          <xsl:with-param name="day" select="27"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="1"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-julian-day template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-julian-day test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-julian-day">
          <xsl:with-param name="year" select="2000"/>
          <xsl:with-param name="month" select="1"/>
          <xsl:with-param name="day" select="1"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="2451545"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-julian-day template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-julian-day test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-julian-day">
          <xsl:with-param name="julian-day" select="2451545"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'2000-01-01'"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="2000"/>
          <xsl:with-param name="month" select="1"/>
          <xsl:with-param name="day" select="1"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="52"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="2001"/>
          <xsl:with-param name="month" select="1"/>
          <xsl:with-param name="day" select="1"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="1"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="2001"/>
          <xsl:with-param name="month" select="1"/>
          <xsl:with-param name="day" select="8"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="2"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="1996"/>
          <xsl:with-param name="month" select="12"/>
          <xsl:with-param name="day" select="31"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="1"/>
    </xsl:call-template>


    <!-- Most years have 52 weeks, but years that start
    on a Thursday and leap years that start on a Wednesday
    have 53 weeks. -->

    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="1998"/>
          <xsl:with-param name="month" select="12"/>
          <xsl:with-param name="day" select="31"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="53"/>
    </xsl:call-template>


    <xsl:if test="$debug='true'"><xsl:message>Test dt:calculate-week-number template</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:calculate-week-number test</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:calculate-week-number">
          <xsl:with-param name="year" select="1992"/>
          <xsl:with-param name="month" select="12"/>
          <xsl:with-param name="day" select="31"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="53"/>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:format-date-time with prepended zeroes</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:format-date-time test (with default format)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:format-date-time">
          <xsl:with-param name="year" select="2001"/>
          <xsl:with-param name="month" select="'01'"/>
          <xsl:with-param name="day" select="'01'"/>
          <xsl:with-param name="hour" select="'01'"/>
          <xsl:with-param name="minute" select="'01'"/>
          <xsl:with-param name="second" select="'01'"/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'2001-01-01T01:01:01'"/>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-year</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-year test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-year">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">2001</xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-month-number</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-month-number test (full)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-month-number">
          <xsl:with-param name="month" select='"january"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">1</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-month-number test (abbrev)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-month-number">
          <xsl:with-param name="month" select='"Jun"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">6</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-year test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-year">
          <xsl:with-param name="xsd-date-time" select='"32001-04-01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">32001</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-year test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-year">
          <xsl:with-param name="xsd-date-time" select='"03:02:01-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-month</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-month test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-month">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">04</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-month test (date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-month">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">04</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-month test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-month">
          <xsl:with-param name="xsd-date-time" select='"32001-04-01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">04</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-month test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-month">
          <xsl:with-param name="xsd-date-time" select='"03:02:01-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-day</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-day test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-day">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T03:02:01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">05</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-day test (date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-day">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">05</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-day test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-day">
          <xsl:with-param name="xsd-date-time" select='"32001-04-05"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">05</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-day test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-day">
          <xsl:with-param name="xsd-date-time" select='"03:02:01-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-hour</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-hour test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-hour">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">06</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-hour test (date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-hour">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-hour test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-hour">
          <xsl:with-param name="xsd-date-time" select='"32001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">06</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-hour test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-hour">
          <xsl:with-param name="xsd-date-time" select='"06:07:08.1234-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">06</xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-minute</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-minute test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-minute">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">07</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-minute test (date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-minute">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-minute test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-minute">
          <xsl:with-param name="xsd-date-time" select='"32001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">07</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-minute test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-minute">
          <xsl:with-param name="xsd-date-time" select='"06:07:08.1234-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">07</xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-second</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (datetime)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (datetime w/- UTC timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234Z"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (datetime w/- + timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234+08:30"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (datetime w/- - timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05T06:07:08.1234-05:45"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"2001-04-05"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"></xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (far future date)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"32001-04-05T06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (time)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"06:07:08.1234"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (time w/- timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"06:07:08.1234-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08.1234</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (time, no fraction)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"06:07:08"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08</xsl:with-param>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-second test (time, no fraction w/- timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-second">
          <xsl:with-param name="xsd-date-time" select='"06:07:08-01:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect">08</xsl:with-param>
    </xsl:call-template>

    <xsl:if test="$debug='true'"><xsl:message>Test dt:get-xsd-datetime-timezone</xsl:message></xsl:if>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (datetime, no timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (date, no timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (time, no timezone)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"03:02:01"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (datetime, UTC)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01Z"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'Z'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (date, UTC)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01Z"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'Z'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (time, UTC)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"03:02:01Z"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'Z'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (datetime, +)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01+10:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'+10:00'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (date, +)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01+10:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'+10:00'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (time, +)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"03:02:01+10:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'+10:00'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (datetime, -)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01T03:02:01-02:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'-02:00'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (date, -)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"2001-04-01-02:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'-02:00'"/>
    </xsl:call-template>

    <xsl:call-template name="test">
      <xsl:with-param name="description">dt:get-xsd-datetime-timezone test (time, -)</xsl:with-param>
      <xsl:with-param name="result">
        <xsl:call-template name="dt:get-xsd-datetime-timezone">
          <xsl:with-param name="xsd-date-time" select='"03:02:01-02:00"'/>
        </xsl:call-template>
      </xsl:with-param>
      <xsl:with-param name="expect" select="'-02:00'"/>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
