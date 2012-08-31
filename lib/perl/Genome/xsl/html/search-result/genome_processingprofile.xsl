<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_processingprofile" match="object[./types[./isa[@type='Genome::ProcessingProfile']]]">
    <div class="search_result">
      <div class="result_icon genome_processingprofile_32">
        <br/>
      </div>
      <div class="result">
        <h3>Processing Profile: <xsl:call-template name="object_link"/></h3>
        <p class="result_summary">
          <strong>Type: </strong><xsl:value-of select="aspect[@name='type_name']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
