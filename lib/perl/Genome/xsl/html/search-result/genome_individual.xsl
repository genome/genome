<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_individual" match="object[./types[./isa[@type='Genome::Individual']]]">
    <div class="search_result">
      <div class="result_icon genome_individual_32">
        <br/>
      </div>
      <div class="result">
        <h3>Individual:
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext">
            <xsl:choose>
              <xsl:when test="normalize-space(aspect[@name='common_name']/value)">
                <xsl:value-of select="aspect[@name='common_name']/value"/>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="@id"/>
              </xsl:otherwise>
            </xsl:choose>
          </xsl:with-param>
        </xsl:call-template>
        </h3>
        <p class="result_summary">
          <xsl:choose>
            <xsl:when test="string(aspect[@name='name']/value)">
              <strong>Name: </strong><xsl:value-of select="aspect[@name='name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <strong>Name: </strong> --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <xsl:choose>
            <xsl:when test="string(aspect[@name='gender']/value)">
              <strong>Gender: </strong><xsl:value-of select="aspect[@name='gender']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <strong>Gender: </strong> --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <xsl:choose>
            <xsl:when test="string(aspect[@name='species_name']/value)">
              <strong>Species: </strong><xsl:value-of select="aspect[@name='species_name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <strong>Species: </strong> --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <xsl:choose>
            <xsl:when test="string(aspect[@name='description']/value)">
              <strong>Description: </strong><xsl:value-of select="aspect[@name='description']/value"/>
            </xsl:when>
            <xsl:otherwise>
              <strong>Description: </strong> --
            </xsl:otherwise>
          </xsl:choose>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>