<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- initializes the dataTable plugin for model set views -->
  <xsl:template name="genome_model_build_set_table_init" match="object[./types[./isa[@type='Genome::Model::Build']]]" mode="set_table_init">
    <xsl:comment>template: status/genome_model_build.xsl match: object[./types[./isa[@type='Genome::Model::Build']]] mode: set_table_init</xsl:comment>
    <script type="text/javascript" charset="utf-8" src="/res/js/pkg/TableTools/media/ZeroClipboard/ZeroClipboard.js"></script>
    <script type="text/javascript" charset="utf-8" src="/res/js/pkg/TableTools/media/js/TableTools.js"></script>
    <style type="text/css">
      @import "/res/js/pkg/TableTools/media/css/TableTools.css";
    </style>
    <script type="text/javascript">
      <xsl:text disable-output-escaping="yes">
        <![CDATA[
                 $(document).ready(function() {
                  TableToolsInit.sSwfPath = "/res/js/pkg/TableTools/media/swf/ZeroClipboard.swf";
                 window.setTable = $('#set').dataTable({
                 "sDom": 'T<"clear">lfrtip',
                 "sScrollX": "100%",
                 "sScrollInner": "110%",
                 "bJQueryUI": true,
                 "sPaginationType": "full_numbers",
                 "bStateSave": true,
                 "iDisplayLength": 25
                 })
                 });
        ]]>
      </xsl:text>
    </script>
  </xsl:template>

  <!-- describes the columns for model build set metric views -->
  <xsl:template name="genome_model_build_set_header" match="aspect[@name='members']" mode="set_header">
  <xsl:comment>template: status/genome_model_build.xsl match: aspect[@name='members'] mode: set_header</xsl:comment>
    <tr>
      <th>
        build ID
      </th>
      <xsl:for-each select="object[aspect[@name='status'] = 'Succeeded'][1]">
        <xsl:for-each select="aspect[@name='processing_profile']/object/aspect[@name='params']/object" >
          <xsl:sort select="aspect[@name='name']/value"/>
          <th>
            <xsl:value-of select="aspect[@name='name']/value" />
          </th>
        </xsl:for-each>
        <th>
          <xsl:value-of select="aspect[@name='inputs']/object/aspect[@name='name']/value" />
        </th>
        <xsl:for-each select="aspect[@name='metrics']/object">
          <xsl:sort select="aspect[@name='name']/value"/>
          <th>
            <xsl:value-of select="aspect[@name='name']/value" />
          </th>
        </xsl:for-each>
      </xsl:for-each>
    </tr>
  </xsl:template>

  <!-- describes the row for model set views -->
  <xsl:template name="genome_model_build_set_row" match="aspect[@name='members']" mode="set_row">
    <xsl:for-each select="object[aspect[@name='status'] = 'Succeeded']">
      <tr>
        <td>
          <xsl:value-of select="@id" />
        </td>
        <xsl:for-each select="aspect[@name='processing_profile']/object/aspect[@name='params']/object" >
          <xsl:sort select="aspect[@name='name']/value"/>
          <td>
            <xsl:value-of select="aspect[@name='value']/value" />
          </td>
        </xsl:for-each>
        <td>
          <xsl:value-of select="aspect[@name='inputs']/object/aspect[@name='value']/object/display_name" />
        </td>
        <xsl:for-each select="aspect[@name='metrics']/object">
          <xsl:sort select="aspect[@name='name']/value"/>
          <td>
            <xsl:value-of select="aspect[@name='value']/value" />
          </td>
        </xsl:for-each>
      </tr>
    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet>
