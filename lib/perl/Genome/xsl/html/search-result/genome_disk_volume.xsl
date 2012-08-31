<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_disk_volume" match="object[./types[./isa[@type='Genome::Disk::Volume']]]">
    <xsl:variable name="href">
        <xsl:text>/disk/allocationmonitor/listallocations?dv_id=</xsl:text><xsl:value-of select='@id'/>
    </xsl:variable>

    <div class="search_result">
      <div class="result_icon genome_disk_volume_32">
        <br/>
      </div>
      <div class="result">
        <h3>Disk Volume: <a><xsl:attribute name="href"><xsl:value-of select="$href"/></xsl:attribute><xsl:value-of select="aspect[@name='mount_path']/value" /></a>
        </h3>

        <p class="result_summary">
          <strong>Disk Group: </strong><xsl:value-of select="aspect[@name='disk_group_names']/value"/>

        </p>
      </div>
    </div> <!-- end search_result -->

  </xsl:template>

</xsl:stylesheet>
