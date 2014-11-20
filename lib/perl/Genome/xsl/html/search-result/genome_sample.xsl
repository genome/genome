<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_sample" match="object[./types[./isa[@type='Genome::Sample']]]">
    <div class="search_result">
      <div class="result_icon genome_sample_32">
        <br/>
      </div>
      <div class="result">
        <h3>Sample: <xsl:call-template name="object_link"/></h3>
        <p class="result_summary">
          <strong>Extraction Type: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='extraction_type']/value)">
              <xsl:value-of select="aspect[@name='extraction_type']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
               --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>Common Name: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='common_name']/value)">
              <xsl:value-of select="aspect[@name='common_name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
               --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>Extraction Name: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='extraction_name']/value)">
              <xsl:value-of select="aspect[@name='extraction_name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>Tissue Description: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='tissue_desc']/value)">
              <xsl:value-of select="aspect[@name='tissue_desc']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
               --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>Organ Name: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='organ_name']/value)">
              <xsl:value-of select="aspect[@name='organ_name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
               --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>