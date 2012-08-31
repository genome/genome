<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_project" match="object[./types[./isa[@type='Genome::Project']]]">
    <div class="search_result">
      <div class="result_icon genome_project_32">
        <br/>
      </div>
      <div class="result">
        <h3>Project: <xsl:call-template name="object_link"/></h3>
        <p class="result_summary">
          <strong>Project Type: </strong> <xsl:value-of select="aspect[@name='project_type']/value"/>
          <xsl:text>; </xsl:text>

          <strong>Status: </strong> <xsl:value-of select="aspect[@name='status']/value"/>
          <xsl:text>; </xsl:text>

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
