<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">


  <xsl:template name="genome_nomenclature_set" match="object[./types[./isa[@type='UR::Object::Set']]]">

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="set_header">
      <xsl:with-param name="display_name" select="'Nomenclatures'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">

        <hr class="space" style="height: 10px; margin: 0;"/>

        <div class="span-24 last">
          <table id="myTable" width="100%" cellpadding="0" cellspacing="0" border="0" class="dataTable">
            <thead>
             <th>Nomenclature Name</th>
             <th>Template</th>
             <th>Upload Samples</th>
            </thead>
            <tbody>
              <xsl:for-each select="aspect[@name='members']/object">
              <tr>
              <td>
                <a><xsl:attribute name="href">/view/genome/nomenclature/set/create.html#id=<xsl:value-of select='@id'/></xsl:attribute>
                <xsl:value-of select="display_name"/></a>
              </td>
              <td>
                <a><xsl:attribute name="href">/view/genome/nomenclature/detail.csv?id=<xsl:value-of select='@id'/></xsl:attribute>CSV</a>
              </td>
              <td>
                <a><xsl:attribute name="href">/view/genome/subject/set/create.html#nomenclature_name=<xsl:value-of select='display_name'/></xsl:attribute>
                Upload Samples</a>
              </td>
              </tr>
              </xsl:for-each>
            </tbody>
          </table>
        </div>
      </div> <!-- end container -->
    </div> <!-- end content -->

    <style type="text/css">
            .toolbar {
                float: left;
                margin-top: 10px;
                margin-left: 30px;
            }
            div .dataTables_length {
                width: inherit;
            }
    </style>
  <script type="text/javascript">

  <xsl:text disable-output-escaping="yes">
    <![CDATA[
                 $(document).ready(function(){
                    $('#myTable').dataTable({
                        "sScrollX": "100%",
                        "sScrollInner": "110%",
                        "bJQueryUI": true,
                        "sPaginationType": "full_numbers",
                        "bStateSave": true,
                        "iDisplayLength": 25,
                        "sDom": '<"H"l<"toolbar">f>t<"F"ip>' 
                    });

                    $("div.toolbar").html('<a href="/view/genome/nomenclature/set/create.html">create a new nomenclature</a>');
                 });
    ]]>
    </xsl:text>

    </script>

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
