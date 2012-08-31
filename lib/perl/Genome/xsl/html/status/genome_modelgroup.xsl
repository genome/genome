<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_modelgroup" match="object[./types[./isa[@type='Genome::ModelGroup']]]">
    <xsl:comment>template: status/genome_modelgroup.xsl match: object[./types[./isa[@type='Genome::ModelGroup']]]</xsl:comment>

    <script type='text/javascript' src='/res/js/app/status/genome_model-group.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Model Group:'" />
      <xsl:with-param name="display_name" select="./aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_modelgroup_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">
          <xsl:call-template name="genome_modelgroup_attributes_box"/>

          <xsl:for-each select="aspect[@name='convergence_model']/object">
            <xsl:call-template name="genome_model_convergence_box"/>
          </xsl:for-each>

        </div>

        <xsl:for-each select="aspect[@name='models']/object[./types[./isa[@type='Genome::Model']]]">
          <xsl:sort select="aspect[@name='is_default']" order="descending"/>
          <xsl:sort select="aspect[@name='name']"/>
          <xsl:call-template name="genome_model_builds_list_table"/>
        </xsl:for-each>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>


  <xsl:template name="genome_modelgroup_attributes_box" match="object[./types[./isa[@type='Genome::ModelGroup']]]" mode="attributes_box">
    <xsl:comment>template: genome_modelgroup.xsl:genome_model_attributes_box match: object[./types[./isa[@type='Genome::ModelGroup']]]  mode: box</xsl:comment>
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="nontyped span-7 last">Model Group Attributes</h3></div>
        <div class="box_button">

        </div>
      </div>

      <div class="box_content rounded-bottom span-8 last">
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
              <td class="value"><xsl:value-of select="aspect[@name='user_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Model Count:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='model_count']/value"/>
              </td>
            </tr>

              <tr>
                <td class="name">Coverage Report:
                </td>
                <td class="value">
                  <xsl:call-template name="object_link_button">
                    <xsl:with-param name="perspective" select="'coverage'" />
                    <xsl:with-param name="linktext" select="'coverage report'" />
                    <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                  </xsl:call-template>
                </td>
              </tr>
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>


</xsl:stylesheet>
