<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_wiki_document" match="object[./types[./isa[@type='Genome::Wiki::Document']]]">
    <div class="search_result">
      <div class="result_icon genome_wiki_document_32">
        <br/>
      </div>
      <div class="result">
        <h3>Wiki: <a><xsl:attribute name="href">https://gscweb.gsc.wustl.edu/wiki/<xsl:value-of select="aspect[@name='title']"/></xsl:attribute><xsl:value-of select="aspect[@name='title']/value"/></a></h3>
        <p class="result_summary">
          <strong>Author: </strong><xsl:value-of select="aspect[@name='user']/value"/>
          <xsl:text>; </xsl:text>
          <strong>Timestamp: </strong><xsl:value-of select="aspect[@name='timestamp']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
