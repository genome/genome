<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- full page display for a project -->
  <xsl:template name="genome_project" match="object[./types[./isa[@type='Genome::Project']]]">
    <xsl:comment>template: status/genome_project.xsl match: object[./types[./isa[@type='Genome::Project']]]</xsl:comment>
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Project'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_project_32'" />
    </xsl:call-template>

    <script type="text/javascript" src="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.js"></script>
    <script type="text/javascript" src="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.plugin.formatted-num.js"></script>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/css/gc_table.css" type="text/css" media="screen"></link>

    <script type='text/javascript' src='/res/js/app/genome_project.js'></script>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this project -->

        <div id="dialog-confirm" title="Removing selected items">
            <p><span class="ui-icon ui-icon-alert" style="float:left; margin:0 7px 20px 0;"></span>
            You are about to remove the selected items from this project (the items will
            not be deleted).
            </p>
        </div>

          <div class="container">
            <table id="parts_table" width="100%" cellspacing="0" cellpadding="0" border="0" style="margin-top: 0;" class="list display">
                <xsl:for-each select="aspect[@name='parts']/object">
                    <tr>
                        <td>
                            <input class="partCheckbox" type="checkbox">
                            <xsl:attribute name="id"><xsl:value-of select="@id"/></xsl:attribute>
                            </input>
                        </td>
                        <td><xsl:value-of select="aspect[@name='entity_class_name_pretty']/value"/></td>
                        <td>
                                <xsl:for-each select="aspect[@name='entity']/object">
                                    <xsl:call-template name="object_link">
                                        <xsl:with-param name="linktext">
                                            <xsl:value-of select="display_name"/>
                                        </xsl:with-param>
                                    </xsl:call-template>
                                </xsl:for-each>
                        </td>
                    </tr>
                </xsl:for-each>
            </table>
           </div>

        </div> <!-- end .objects -->

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
