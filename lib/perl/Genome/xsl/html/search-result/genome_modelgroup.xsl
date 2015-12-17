<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_modelgroup" match="object[./types[./isa[@type='Genome::ModelGroup']]]">
    <div class="search_result">
      <div class="result_icon genome_modelgroup_32">
        <br/>
      </div>
      <div class="result">
        <h3>Model Group: <xsl:call-template name="object_link"/></h3>
        <p class="result_summary">
          <strong>Members: </strong>
          <xsl:for-each select="aspect[@name='models']/object">
            <xsl:value-of select="aspect[@name='name']/value"/>
            <xsl:if test="position() != last()"><xsl:text>, </xsl:text></xsl:if>
          </xsl:for-each>
        </p>
      </div>
    </div> <!-- end search_result -->
  </xsl:template>

</xsl:stylesheet>
