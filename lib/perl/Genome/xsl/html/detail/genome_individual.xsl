<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_individual" match="object[./types[./isa[@type='Genome::Individual']]]">
    <xsl:comment>template: /html/detail/genome_individual.xsl  match: object[./types[./isa[@type='Genome::Individual']]]</xsl:comment>

    <script type="text/javascript" src="/res/js/pkg/dataTables/media/js/jquery.dataTables.js"></script>

    <script type="text/javascript" language="javascript">

    </script>

    <script type="text/javascript" language="javascript">
      window.individual_id = <xsl:value-of select="/object/@id"/>;
    </script>

    <script type="text/javascript" src="/res/js/app/detail/genome_individual.js"></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Individual'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_individual_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div class="box_header span-24 last rounded-top">
          <div class="box_title"><h3 class="nontyped span-24 last">Grouped by nomenclature</h3></div>
        </div>
        <div class="box_content rounded-bottom span-24 last">

          <!-- details for this individual -->
          <table id="individuals" class="lister datatable">
            <thead>
              <tr>
                <th>Nomenclature</th>
                <th>Nomenclature ID</th>
                <th>Attribute</th>
                <th>Value</th>
              </tr>
            </thead>
          </table>
        </div>
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>




