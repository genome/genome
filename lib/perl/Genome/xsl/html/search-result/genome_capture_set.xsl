<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_capture_set" match="object[./types[./isa[@type='Genome::Capture::Set']]]">
    <div class="search_result">
      <div class="result_icon genome_capture_set_32">
        <br/>
      </div>
      <div class="result">
        <h3>Capture Set: <xsl:call-template name="object_link" /></h3>

        <p class="result_summary">
          <strong>Status: </strong><xsl:value-of select="aspect[@name='status']/value"/> --
          <strong>Description: </strong><xsl:value-of select="aspect[@name='description']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
