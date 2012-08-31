<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:str="http://xsltsl.org/string">

  <xsl:template name="genome_disk_volume" match="object[./types[./isa[@type='Genome::Disk::Volume']]]">
    <xsl:comment>template: status/genome_disk_volume.xsl match: object[./types[./isa[@type='Genome::Disk::Volume']]]</xsl:comment>
    <script type="text/javascript" src="/res/js/pkg/protovis.js"></script>
    <script type="text/javascript" src="/res/js/app/status/genome_disk_volume_treemap.js"></script>
    <script type="text/javascript" src="/res/js/pkg/json2.js"></script>
    <script type="text/javascript" src="/res/js/pkg/jquery.tipsy.js"></script>
    <script type="text/javascript" src="/res/js/pkg/protovis.tipsy.js"></script>
    <link href="/res/css/tipsy.css" type="text/css" rel="stylesheet"/>

    <xsl:comment>allocation data JSON</xsl:comment>
    <script type="text/javascript">

      <!--      <xsl:for-each select="//aspect[@name='allocations']">

var allocation_data = {
<xsl:for-each select="object">
{
"owner_id": "<xsl:value-of select="aspect[@name='owner_id']/value"/>",
"owner_class_name": "<xsl:call-template name="str:substring-after-last"><xsl:with-param name="text"> <xsl:value-of select="aspect[@name='owner_class_name']/value"/></xsl:with-param><xsl:with-param name="chars"><xsl:text>Genome::Model::</xsl:text></xsl:with-param></xsl:call-template>",

"kilobytes_requested": <xsl:value-of select="aspect[@name='kilobytes_requested']/value"/>,
<xsl:if test="aspect[@name='build']">
"build_id": "<xsl:value-of select="aspect[@name='build']/object/aspect[@name='build_id']/value"/>",
"model_id": "<xsl:value-of select="aspect[@name='build']/object/aspect[@name='model_id']/value"/>",
"build_status": "<xsl:value-of select="aspect[@name='build']/object/aspect[@name='status']/value"/>",
</xsl:if>
},
</xsl:for-each>
};
</xsl:for-each>
      -->
      <xsl:for-each select="//aspect[@name='allocations']">

        var allocations = [
        <xsl:for-each select="object">
          <xsl:variable name="short_class_name">
            <xsl:call-template name="str:substring-after-last">
              <xsl:with-param name="text"> <xsl:value-of select="aspect[@name='owner_class_name']/value"/> </xsl:with-param>
              <xsl:with-param name="chars"> <xsl:text>Genome::Model::</xsl:text> </xsl:with-param>
            </xsl:call-template>
          </xsl:variable>
          {
          "id": "<xsl:value-of select="@id"/>",
          "owner_id": "<xsl:value-of select="aspect[@name='owner_id']/value"/>",
          "absolute_path": "<xsl:value-of select="aspect[@name='absolute_path']/value"/>",
          "display_name": "<xsl:value-of select="display_name"/>",
          "label_name": "<xsl:value-of select="label_name"/>",
          "owner_class_name": "<xsl:value-of select="aspect[@name='owner_class_name']/value"/>",
          "absolute_path": "<xsl:value-of select="aspect[@name='absolute_path']/value"/>",
          "owner_class_name": "<xsl:value-of select="$short_class_name"/>",
          "kilobytes_requested": "<xsl:value-of select="aspect[@name='kilobytes_requested']/value"/>",
          },
        </xsl:for-each>
        ];
      </xsl:for-each>


      <xsl:for-each select="//aspect[@name='allocations']">

        var allocation_build_info = {
        <xsl:for-each select="object">
          <xsl:choose>
            <xsl:when test="aspect[@name='build']">
              <xsl:for-each select="aspect[@name='build']/object">
                <xsl:variable name="build_url">
                  <xsl:call-template name="object_link_href">
                    <xsl:with-param name="type" select="@type"/>
                    <xsl:with-param name="id" select="aspect[@name='build_id']/value"/>
                    <xsl:with-param name="perspective" select="'status'"/>
                    <xsl:with-param name="key" select="'id'"/>
                    <xsl:with-param name="toolkit" select="'html'"/>
                  </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="model_url">
                  <xsl:call-template name="object_link_href">
                    <xsl:with-param name="type" select="@type"/>
                    <xsl:with-param name="id" select="aspect[@name='model_id']/value"/>
                    <xsl:with-param name="perspective" select="'status'"/>
                    <xsl:with-param name="key" select="'id'"/>
                    <xsl:with-param name="toolkit" select="'html'"/>
                  </xsl:call-template>
                </xsl:variable>

                "<xsl:value-of select="../../@id"/>": {
                "model": '<a><xsl:attribute name="href"><xsl:value-of select="$model_url"/></xsl:attribute><xsl:attribute name="class">mini btn</xsl:attribute><xsl:value-of select="aspect[@name='model_id']/value"/></a>',
                "build": '<a><xsl:attribute name="href"><xsl:value-of select="$build_url"/></xsl:attribute><xsl:attribute name="class">mini btn</xsl:attribute><xsl:value-of select="aspect[@name='build_id']/value"/></a>',
                "build name": '<xsl:value-of select="display_name"/>',
                "build class": "<xsl:value-of select="@type"/>",
                "status": "<xsl:value-of select="aspect[@name='status']/value"/>",
                "run by": "<xsl:value-of select="aspect[@name='run_by']/value"/>",
                "date scheduled": "<xsl:value-of select="aspect[@name='date_scheduled']/value"/>",
                "date completed": "<xsl:value-of select="aspect[@name='date_completed']/value"/>",
                },
              </xsl:for-each>
            </xsl:when>
            <xsl:otherwise>
              "<xsl:value-of select="@id"/>": {
              "build_id": "--",
              "model_id": "--",
              "status": "--",
              "run_by": "--",
              "date_scheduled": "--",
              "date_completed": "--",
              },
            </xsl:otherwise>
          </xsl:choose>
        </xsl:for-each>
        };
      </xsl:for-each>

      var total_allocation = {
      "total_kb": <xsl:value-of select="aspect[@name='total_kb']/value"/>,
      "unallocated_kb": <xsl:value-of select="aspect[@name='unallocated_kb']/value"/>
      }

    </script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Disk Volume:'" />
      <xsl:with-param name="display_name" select="aspect[@name='mount_path']/value" />
      <xsl:with-param name="icon" select="'genome_disk_volume_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div class="box rounded">
          <div style="width: 25%; margin-right; 2%; display: inline-block">
            <table border="0" cellpadding="0" cellspacing="0" class="name-value" style="margin:0;">
              <tr>
                <td class="name">disk group:</td>
                <td class="value"><xsl:value-of select="aspect[@name='disk_group_names']/value"/></td>
              </tr>

              <tr>
                <td class="name">status:</td>
                <td class="value"><xsl:value-of select="aspect[@name='disk_status']/value"/></td>
              </tr>

              <tr>
                <td class="name">can allocate:</td>
                <td class="value"><xsl:value-of select="aspect[@name='can_allocate']/value"/></td>
              </tr>

            </table>
          </div>

          <div style="width: 50%; display: inline-block;">
            <table border="0" cellpadding="0" cellspacing="0" class="name-value" style="margin:0;">
              <tr>
                <td class="name">unallocated (kb):</td>
                <td class="value"><xsl:value-of select="format-number(aspect[@name='unallocated_kb']/value, '#,##0')"/></td>
              </tr>
              <tr>
                <td class="name">total (kb):</td>
                <td class="value"><xsl:value-of select="format-number(aspect[@name='total_kb']/value, '#,##0')"/></td>
              </tr>
              <tr>
                <td class="name">&#160;</td>
                <td class="value">&#160;</td>
              </tr>

            </table>

          </div>
        </div>

        <div class="generic_lister">
          <div class="box_header span-24 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-24 last">Disk Allocation Treemap</h3></div>
          </div>
          <div class="box_content rounded-bottom span-24 last">
            <div style="float: left; height: 600px; margin-bottom: 10px;border-bottom: 1px solid #cdcdcd;">
              <script type="text/javascript">
                render_treemap(allocations, 950, 600, allocation_build_info, total_allocation);
              </script>
            </div>
          </div>
        </div>

        <xsl:call-template name="genome_disk_volume_table"></xsl:call-template>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="genome_disk_volume_table">
    <xsl:comment>template: status/genome_disk_volume.xsl name: genome_disk_volume_table</xsl:comment>
    <table id="set" class="dataTable">
      <thead>
        <tr>
          <th>owner</th>
          <th>owner class</th>
          <th>build status</th>
          <th>requested</th>
          <th>absolute path</th>
        </tr>
      </thead>
      <tbody>

        <xsl:for-each select="/object/aspect[@name='allocations']/object">
          <xsl:sort select="aspect[@name='owner_id']/value" data-type="number" order="ascending"/>
          <xsl:call-template name="genome_disk_volume_table_row"/>
        </xsl:for-each>
      </tbody>
    </table>

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

  <xsl:template name="genome_disk_volume_table_row">
    <xsl:comment>template: status/genome_disk_volume.xsl name: genome_disk_volume_table_row</xsl:comment>
    <tr>
      <td>
        <xsl:choose>
          <xsl:when test="aspect[@name='build']">
            <xsl:for-each select="aspect[@name='build']/object">
              <xsl:call-template name="object_link_button">
                <xsl:with-param name="linktext" select="@id" />
                <xsl:with-param name="icon" select="'sm-icon-extlink'" />
              </xsl:call-template>
            </xsl:for-each>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="aspect[@name='owner_id']/value"/>
          </xsl:otherwise>
        </xsl:choose>

      </td>
      <td>
        <xsl:call-template name="str:substring-after-last">
          <xsl:with-param name="text"> <xsl:value-of select="aspect[@name='owner_class_name']/value"/> </xsl:with-param>
          <xsl:with-param name="chars"> <xsl:text>Genome::Model::</xsl:text> </xsl:with-param>
        </xsl:call-template>
      </td>
      <td>
        <xsl:choose>
          <xsl:when test="aspect[@name='build']">
            <xsl:variable name="e_status" select="aspect[@name='build']/object/aspect[@name='status']/value"/>
            <xsl:variable name="lc_e_status" select="translate($e_status,'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz')"/>
            <xsl:attribute name="class"><xsl:text disable-output-escaping="yes">status </xsl:text><xsl:value-of select="$lc_e_status"/></xsl:attribute>
            <xsl:value-of select="$lc_e_status"/>
          </xsl:when>
          <xsl:otherwise>
            --
          </xsl:otherwise>
        </xsl:choose>
      </td>

      <td>
        <xsl:variable name="kb_req" select="aspect[@name='kilobytes_requested']/value"/>
        <xsl:value-of select="format-number($kb_req, '#,##0')"/>
      </td>

      <td>
        <xsl:variable name="absolute_path" select="aspect[@name='absolute_path']/value"/>
        <xsl:value-of select="$absolute_path"/>
        <!-- <a> -->
        <!--   <xsl:attribute name="href"> -->
        <!--     <xsl:value-of select="$absolute_path"/> -->
        <!--   </xsl:attribute> -->
        <!--   <xsl:value-of select="substring($absolute_path,1,30)"/>... -->
        <!-- </a> -->
      </td>
    </tr>
  </xsl:template>

</xsl:stylesheet>
