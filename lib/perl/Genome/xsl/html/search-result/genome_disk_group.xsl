<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_disk_group" match="object[./types[./isa[@type='Genome::Disk::Group']]]">
    <xsl:variable name="href">
      <xsl:text>/disk/allocationmonitor/listvolumes?disk_group_name=</xsl:text><xsl:value-of select="normalize-space(aspect[@name='disk_group_name']/value)" />
    </xsl:variable>

    <div class="search_result">
      <div class="result_icon genome_disk_group_32">
        <br/>
      </div>
      <div class="result">
        <h3>Disk Group: <a><xsl:attribute name="href"><xsl:value-of select="$href"/></xsl:attribute><xsl:value-of select="aspect[@name='disk_group_name']/value" /></a></h3>
        <p class="result_summary">
          <strong>User Name: </strong><xsl:value-of select="aspect[@name='user_name']/value"/>
          <xsl:text>; </xsl:text>
          <strong>Group Name: </strong><xsl:value-of select="aspect[@name='group_name']/value"/>
        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
