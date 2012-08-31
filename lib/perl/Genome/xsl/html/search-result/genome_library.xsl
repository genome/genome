<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_library" match="object[./types[./isa[@type='Genome::Library']]]">
    <div class="search_result">
      <div class="result_icon genome_library_32">
        <br/>
      </div>
      <div class="result">
        <h3>Library: <xsl:call-template name="object_link"/></h3>
        <p class="result_summary">
          <strong>Taxon: </strong><xsl:value-of select="aspect[@name='taxon']/object/aspect[@name='species_name']/value"/>
          <xsl:text>; </xsl:text>
          <strong>Sample: </strong><xsl:value-of select="aspect[@name='sample']/object/aspect[@name='name']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->
  </xsl:template>

</xsl:stylesheet>
