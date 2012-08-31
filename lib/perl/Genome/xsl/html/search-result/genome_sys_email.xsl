<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_sys_email" match="object[./types[./isa[@type='Genome::Sys::Email']]]">
    <xsl:variable name="message_url_base">
      <xsl:value-of select="normalize-space(aspect[@name='mail_server_path']/value)"/>
      <xsl:text>/</xsl:text>
      <xsl:value-of select="normalize-space(aspect[@name='list_name']/value)"/>
      <xsl:text>/</xsl:text>
      <xsl:value-of select="normalize-space(aspect[@name='month']/value)"/>
      <xsl:text>/</xsl:text>
    </xsl:variable>

    <div class="search_result">
      <div class="result_icon genome_sys_email_32">
        <br/>
      </div>
      <div class="result">
        <h3>Email:
        <a>
          <xsl:attribute name="href">
            <xsl:value-of select="$message_url_base"/>
            <xsl:value-of select="normalize-space(aspect[@name='message_id']/value)" />
            <xsl:text>.html</xsl:text>
          </xsl:attribute>
          <xsl:value-of select="aspect[@name='subject']/value"/>
        </a>
        </h3>
        <p class="resource_buttons">
          <a class="mini btn">
            <xsl:attribute name="href">
              <xsl:value-of select="normalize-space(aspect[@name='mail_list_path']/value)"/>
              <xsl:text>/</xsl:text>
              <xsl:value-of select="normalize-space(aspect[@name='list_name']/value)"/>
            </xsl:attribute>
            <span class="sm-icon sm-icon-extlink"><br/></span>
            <xsl:value-of select="normalize-space(aspect[@name='list_name']/value)"/> listinfo
          </a>

          <a class="mini btn">
            <xsl:attribute name="href">
              <xsl:value-of select="$message_url_base"/>
              <xsl:text>date.html</xsl:text>
            </xsl:attribute>
            <span class="sm-icon sm-icon-extlink"><br/></span>
            <xsl:value-of select="normalize-space(aspect[@name='list_name']/value)"/><xsl:text> </xsl:text>
            <xsl:value-of select="normalize-space(aspect[@name='month']/value)"/> archives
          </a>

        </p>
        <p class="result_summary">
          <xsl:value-of select="aspect[@name='blurb']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
