<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--  <xsl:template name="genome_sys_google_spreadsheet" match="object[./types[./isa[@type='Genome::Sys::Google::Spreadsheet']]]">
<xsl:template name="genome_sys_google_spreadsheet" match="doc[./field[@name='class']]='Genome::Sys::Google::Spreadsheet'">
-->

<xsl:template name="genome_sys_google_spreadsheet" match="/doc[field[@name='class']='Genome::Sys::Google::Spreadsheet']">
<xsl:variable name="fullContent" select="/doc/field[@name='content']"/>

    <div class="search_result">
      <div class="result_icon genome_sys_google_spreadsheet_32">
        <br/>
      </div>
      <div class="result">
        <h3>Spreadsheet: <xsl:value-of select="/doc/field[@name='title']"/>
        </h3>
        <p class="resource_buttons">
          <a class="mini btn"><xsl:attribute name="href">https://somewhere/<xsl:value-of select="/doc/field[@name='object_id']"/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>Google doc</a>

        </p>
        <p class="result_summary">
            <strong>modified: </strong> <xsl:value-of select="/doc/field[@name='timestamp']"/><xsl:text>; </xsl:text>
            <strong>preview: </strong> <xsl:value-of select='substring($fullContent, 0, 200)'/><xsl:text>; </xsl:text>
        </p>
      </div>
    </div> <!-- end search_result -->

</xsl:template>

</xsl:stylesheet>
