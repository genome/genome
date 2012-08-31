<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_populationgroup" match="object[./types[./isa[@type='Genome::PopulationGroup']]]">
    <div class="search_result">
      <div class="result_icon genome_populationgroup_32">
        <br/>
      </div>
      <div class="result">
        <h3>Population Group:
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext">
            <xsl:choose>
              <xsl:when test="normalize-space(aspect[@name='common_name']/value)">
                <xsl:value-of select="aspect[@name='common_name']/value"/>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="aspect[@name='name']/value"/>
              </xsl:otherwise>
            </xsl:choose>
          </xsl:with-param>
        </xsl:call-template>
        </h3>

        <p class="result_summary">
          <xsl:choose>
            <xsl:when test="string(aspect[@name='description']/value)">
              <strong>Description: </strong><xsl:value-of select="aspect[@name='description']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <strong>Description: </strong> --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>Members: </strong>
          <xsl:choose>
            <xsl:when test="count(aspect[@name='members']/object) &gt; 0">
              <xsl:for-each select="aspect[@name='members']/object">
                <xsl:choose>
                  <xsl:when test="string(normalize-space(aspect[@name='common_name']/value))">
                    <xsl:value-of select="aspect[@name='common_name']/value"/>
                  </xsl:when>
                  <xsl:when test="string(aspect[@name='name']/value)">
                    <xsl:value-of select="aspect[@name='name']/value"/>
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="display_name"/>
                  </xsl:otherwise>
                </xsl:choose>
                <xsl:if test="position() != last()"><xsl:text>; </xsl:text></xsl:if>
              </xsl:for-each>
            </xsl:when>
            <xsl:otherwise>
              No members found in this Population Group.
            </xsl:otherwise>
          </xsl:choose>

        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
