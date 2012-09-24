<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">


  <xsl:template name="genome_model_build" match="object[./types[./isa[@type='Genome::Model::Build']]]">

    <xsl:call-template name="control_bar_view"/>

    <xsl:variable name="build_id" select="/object/@id"/>

    <xsl:call-template name="set_header">
      <xsl:with-param name="display_name" select="concat('Allocations for build: ',$build_id)" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">

        <hr class="space" style="height: 10px; margin: 0;"/>
    
        <span style="visibility:hidden" id="build_id"><xsl:value-of select="/object/@id"/></span>
        <div class="span-24 last">
          <table id="myTable" width="100%" cellpadding="0" cellspacing="0" border="0" class="dataTable">
            <thead>
             <th>Kilobytes</th>
             <th>Is Archived?</th>
             <th>Archivable</th>
             <th>Disk Group</th>
             <th>Path</th>
             <th>ID</th>
            </thead>
            <tbody>
              <xsl:for-each select="aspect[@name='members']/object">
              <tr>
              <td><xsl:value-of select="aspect[@name='user_id']"/></td>
              <td><xsl:value-of select="aspect[@name='status']"/></td>
              <td><xsl:value-of select="aspect[@name='time_submitted']"/></td>
              <td><xsl:value-of select="aspect[@name='time_finished']"/></td>
              <td><xsl:value-of select="aspect[@name='time_finished']"/></td>
              </tr>
              </xsl:for-each>
            </tbody>
          </table>
        </div>
      </div> <!-- end container -->
    </div> <!-- end content -->

  <script type="text/javascript">
        $(document).ready(function() {

            var build_id = $("#build_id").html();
            debugger;
            var url = "/viewajax/genome/model/build/data-table.json?id=" + build_id;

            $('#myTable').dataTable( {
                 "sAjaxSource": url,
                 "sScrollX": "100%",
                 "sScrollInner": "110%",
                 "bJQueryUI": true,
                 "sPaginationType": "full_numbers",
                 "bStateSave": true,
                 "iDisplayLength": 25
            } );
        } );
    </script>

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
