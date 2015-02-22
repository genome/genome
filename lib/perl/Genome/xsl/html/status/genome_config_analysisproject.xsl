<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_config_analysisproject" match="object[./types[./isa[@type='Genome::Config::AnalysisProject']]]">
    <xsl:comment>template: status/genome_config_analysisproject.xsl match: object[./types[./isa[@type='Genome::Config::AnalysisProject']]]</xsl:comment>

    <script type='text/javascript' src='/res/js/app/status/genome_model-group.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Analysis Project:'" />
      <xsl:with-param name="display_name" select="./aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_project_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">
          <xsl:call-template name="genome_config_analysisproject_attributes_box"/>
        </div>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="genome_config_analysisproject_attributes_box" match="object[./types[./isa[@type='Genome::Config::AnalysisProject']]]" mode="attributes_box">
<div class="span_12_box_masonry">
  <div class="box_header span-12 last rounded-top">
    <div class="box_title"><h3 class="nontyped span-12 last">Analysis Project Attributes</h3></div>
    <div class="box_button">

    </div>
  </div>

  <div class="box_content rounded-bottom span-12 last">
    <table class="name-value">
      <tbody>
        <tr>
          <td class="name">ID:
          </td>
          <td class="value"><xsl:value-of select="@id"/>
          </td>
        </tr>

        <tr>
          <td class="name">Name:
          </td>
          <td class="value">
            <xsl:choose>
              <xsl:when test="string(normalize-space(aspect[@name='name']/value))">
                <xsl:value-of select="normalize-space(aspect[@name='name']/value)"/>
              </xsl:when>
              <xsl:otherwise>
                --
              </xsl:otherwise>
            </xsl:choose>
          </td>
        </tr>

        <tr>
          <td class="name">User Name:
          </td>
          <td class="value"><xsl:value-of select="aspect[@name='created_by']/value"/>
          </td>
        </tr>

        <tr>
          <td class="name">Created:
          </td>
          <td class="value"><xsl:value-of select="aspect[@name='created_at']/value"/>
          </td>
        </tr>

        <tr>
          <td class="name">Updated:
          </td>
          <td class="value"><xsl:value-of select="aspect[@name='updated_at']/value"/>
          </td>
        </tr>
      </tbody>
    </table>
  </div>
</div>
</xsl:template>


</xsl:stylesheet>
