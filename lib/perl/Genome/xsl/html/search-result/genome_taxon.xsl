<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_taxon" match="object[./types[./isa[@type='Genome::Taxon']]]">
    <div class="search_result">
      <div class="result_icon genome_taxon_32">
        <br/>
      </div>
      <div class="result">
        <h3>Taxon:
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext" select="aspect[@name='species_name']/value"/>
        </xsl:call-template>
        </h3>

        <p class="result_summary">
          <strong>Latin Species Name: </strong>
          <xsl:choose>
            <xsl:when test="string(aspect[@name='species_latin_name']/value)">
              <xsl:value-of select="aspect[@name='species_latin_name']/value"/>
              <xsl:text>; </xsl:text>
            </xsl:when>
            <xsl:otherwise>
              --
              <xsl:text> ; </xsl:text>
            </xsl:otherwise>
          </xsl:choose>

          <strong>NCBI Taxon ID: </strong>
          <xsl:choose>
            <xsl:when test="string(normalize-space(aspect[@name='ncbi_taxon_id']/value))">
              <xsl:value-of select="normalize-space(aspect[@name='ncbi_taxon_id']/value)"/>
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
