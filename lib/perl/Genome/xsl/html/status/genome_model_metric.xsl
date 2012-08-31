<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- initializes the dataTable plugin for model set views -->
  <xsl:template name="genome_model_metric_set_table_init" match="object[./types[./isa[@type='Genome::Model::Metric']]]" mode="set_table_init">
    <xsl:comment>template: status/genome_model_metric.xsl match: object[./types[./isa[@type='Genome::Model::Metric']]] mode: set_table_init</xsl:comment>
    <script type="text/javascript">
      <xsl:text disable-output-escaping="yes">
        <![CDATA[
                 $(document).ready(
                 window.setTable = $('#set').dataTable({
                 "sScrollX": "100%",
                 "sScrollInner": "110%",
                 "bJQueryUI": true,
                 "sPaginationType": "full_numbers",
                 "bStateSave": true,
                 "iDisplayLength": 25
                 })
                 );
        ]]>
      </xsl:text>
    </script>
  </xsl:template>

  <!-- describes the columns for model set views -->
  <xsl:template name="genome_model_metric_set_header" match="object[./types[./isa[@type='Genome::Model::Metric']]]" mode="set_header">
    <xsl:comment>template: status/genome_model_metric.xsl match: object[./types[./isa[@type='Genome::Model::Metric']]] mode: set_header</xsl:comment>
    <tr>
      <th>
        name
      </th>
      <th>
        value
      </th>
    </tr>
  </xsl:template>

  <!-- describes the row for model set views -->
  <xsl:template name="genome_model_metric_set_row" match="object[./types[./isa[@type='Genome::Model::Metric']]]" mode="set_row">
    <xsl:comment>template: status/genome_model_metric.xsl match: object[./types[./isa[@type='Genome::Model::Metric']]] mode: set_row</xsl:comment>

    <tr>
      <td>
        <xsl:value-of select="aspect[@name='name']"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='value']"/>
      </td>
    </tr>
  </xsl:template>

</xsl:stylesheet>
