<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- full page display for a project -->
  <xsl:template name="genome_site_wugc_project" match="object[./types[./isa[@type='Genome::Site::TGI::Project']]]">
    <xsl:comment>template: status/genome_project.xsl match: object[./types[./isa[@type='Genome::Site::TGI::Project']]]</xsl:comment>
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Project:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_project_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this project -->
          <div class="span_8_box_masonry">
            <div class="box_header span-8 last rounded-top">
              <div class="box_title"><h3 class="nontyped span-7 last">Project Attributes</h3></div>
              <div class="box_button">

              </div>
            </div>

            <div class="box_content rounded-bottom span-8 last">
              <table class="name-value">
                <tbody>
                  <tr>
                    <td class="name">Name:</td>
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
                    <td class="name">Project Type:</td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='project_type']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='project_type']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Description:</td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='description']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='description']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">External Contact:</td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='external_contact']/object/aspect[@name='name']/value))">
                          <a>
                            <xsl:attribute name="href">mailto:<xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='email']/value)"/></xsl:attribute>
                            <xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='name']/value)"/>
                          </a>
                          (<xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='type']/value)"/>)
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Internal Contact:</td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='internal_contact']/object/aspect[@name='name']/value))">
                          <a>
                            <xsl:attribute name="href">mailto:<xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='email']/value)"/></xsl:attribute>
                            <xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='name']/value)"/>
                          </a>
                          (<xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='type']/value)"/>)
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                </tbody>
              </table>
            </div>
          </div>

          <xsl:for-each select="aspect[@name='samples']/object">
            <xsl:call-template name="genome_sample_box"/>
          </xsl:for-each>

        </div> <!-- end .objects -->


        <xsl:for-each select="aspect[@name='models']">
          <xsl:call-template name="genome_model_build_table"/>
        </xsl:for-each>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

      <!-- box element for project, intended for display in a jquery masonry layout -->
    <xsl:template name="genome_project_box">
      <xsl:comment>template: /html/status/genome_project.xsl; name="genome_project_box"</xsl:comment>
      <div class="span_8_box_masonry">
        <div class="box_header span-8 last rounded-top">
          <div class="box_title"><h3 class="genome_sample_16 span-7 last">Project</h3></div>
          <div class="box_button">
            <xsl:call-template name="object_link_button_tiny">
              <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
            </xsl:call-template>
          </div>
        </div>

        <div class="box_content rounded-bottom span-8 last">
          <table class="name-value">
            <tbody>
              <tr>
                <td class="name">Name:</td>
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
                <td class="name">Project Type:</td>
                <td class="value">
                  <xsl:choose>
                    <xsl:when test="string(normalize-space(aspect[@name='project_type']/value))">
                      <xsl:value-of select="normalize-space(aspect[@name='project_type']/value)"/>
                    </xsl:when>
                    <xsl:otherwise>
                      --
                    </xsl:otherwise>
                  </xsl:choose>
                </td>
              </tr>

              <tr>
                <td class="name">Description:</td>
                <td class="value">
                  <xsl:choose>
                    <xsl:when test="string(normalize-space(aspect[@name='description']/value))">
                      <xsl:value-of select="normalize-space(aspect[@name='description']/value)"/>
                    </xsl:when>
                    <xsl:otherwise>
                      --
                    </xsl:otherwise>
                  </xsl:choose>
                </td>
              </tr>

              <tr>
                <td class="name">External Contact:</td>
                <td class="value">
                  <xsl:choose>
                    <xsl:when test="string(normalize-space(aspect[@name='external_contact']/object/aspect[@name='name']/value))">
                      <a>
                        <xsl:attribute name="href">mailto:<xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='email']/value)"/></xsl:attribute>
                        <xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='name']/value)"/>
                      </a>
                      (<xsl:value-of select="normalize-space(aspect[@name='external_contact']/object/aspect[@name='type']/value)"/>)
                    </xsl:when>
                    <xsl:otherwise>
                      --
                    </xsl:otherwise>
                  </xsl:choose>
                </td>
              </tr>

              <tr>
                <td class="name">Internal Contact:</td>
                <td class="value">
                  <xsl:choose>
                    <xsl:when test="string(normalize-space(aspect[@name='internal_contact']/object/aspect[@name='name']/value))">
                      <a>
                        <xsl:attribute name="href">mailto:<xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='email']/value)"/></xsl:attribute>
                        <xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='name']/value)"/>
                      </a>
                      (<xsl:value-of select="normalize-space(aspect[@name='internal_contact']/object/aspect[@name='type']/value)"/>)
                    </xsl:when>
                    <xsl:otherwise>
                      --
                    </xsl:otherwise>
                  </xsl:choose>
                </td>
              </tr>

            </tbody>
          </table>
        </div>
      </div>
    </xsl:template>

</xsl:stylesheet>
