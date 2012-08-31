<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_softwareresult" match="object[./types[./isa[@type='Genome::SoftwareResult']]]">
    <xsl:comment>template: /html/status/genome_sample.xsl  match: object[./types[./isa[@type='Genome::SoftwareResult']]]</xsl:comment>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Software Result:'" />
      <xsl:with-param name="display_name" select="./display_name" />
      <xsl:with-param name="icon" select="'genome_result_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this sample -->
          <div class="span_12_box_masonry">
            <div class="box_header span-12 last rounded-top">
              <div class="box_title"><h3 class="nontyped span-7 last">Attributes</h3></div>
              <div class="box_button">

              </div>
            </div>

            <div class="box_content rounded-bottom span-12 last">
              <table class="name-value">
                <tbody>

                  <tr>
                    <td class="name">ID:
                    </td>
                    <td class="value">
                      <xsl:value-of select="normalize-space(@id)"/>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Subclass:
                    </td>
                    <td class="value">
                      <xsl:value-of select="normalize-space(aspect[@name='subclass_name']/value)"/>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Output Directory:
                    </td>
                    <td class="value">
                      <xsl:value-of select="normalize-space(aspect[@name='output_dir']/value)"/>
                    </td>
                    </tr>

                    <tr>
                    <td class="name">Test Name:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='test_name']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='test_name']/value)"/>
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

          <xsl:if test="count(aspect[@name='params']/object) > 0">
            <div class="span_12_box_masonry">
              <div class="box_header span-12 last rounded-top">
                <div class="box_title"><h3 class="nontyped span-7 last">Parameters</h3></div>
                <div class="box_button">
                </div>
              </div>

              <div class="box_content rounded-bottom span-12 last">
                <table class="name-value">
                  <tbody>
                    <xsl:for-each select="aspect[@name='params']/object">
                      <tr>
                        <td class="name">
                          <xsl:value-of select="aspect[@name='param_name']/value"/>
                        </td>
                        <td class="value">
                          <xsl:value-of select="aspect[@name='param_value']/value"/>
                        </td>
                      </tr>
                    </xsl:for-each>
                  </tbody>
                </table>
              </div>
            </div>
          </xsl:if>

          <xsl:if test="count(aspect[@name='inputs']/object) > 0">
            <div class="span_12_box_masonry">
              <div class="box_header span-12 last rounded-top">
                <div class="box_title"><h3 class="nontyped span-7 last">Inputs</h3></div>
                <div class="box_button">
                </div>
              </div>

              <div class="box_content rounded-bottom span-12 last">
                <table class="name-value">
                  <tbody>
                    <xsl:for-each select="aspect[@name='inputs']/object">
                      <tr>
                        <td class="name">
                          <xsl:value-of select="aspect[@name='input_name']/value"/>
                        </td>
                        <td class="value">
                          <xsl:value-of select="aspect[@name='input_value']/value"/>
                        </td>
                      </tr>
                    </xsl:for-each>
                  </tbody>
                </table>
              </div>
            </div>
          </xsl:if>

          <xsl:if test="count(aspect[@name='metrics']/object) > 0">
            <div class="span_12_box_masonry">
              <div class="box_header span-12 last rounded-top">
                <div class="box_title"><h3 class="nontyped span-7 last">Metrics</h3></div>
                <div class="box_button">
                </div>
              </div>

              <div class="box_content rounded-bottom span-12 last">
                <table class="name-value">
                  <tbody>
                    <xsl:for-each select="aspect[@name='metrics']/object">
                      <tr>
                        <td class="name">
                          <xsl:value-of select="aspect[@name='metric_name']/value"/>
                        </td>
                        <td class="value">
                          <xsl:value-of select="aspect[@name='metric_value']/value"/>
                        </td>
                      </tr>
                    </xsl:for-each>
                  </tbody>
                </table>
              </div>
            </div>
          </xsl:if>

          <xsl:if test="count(aspect[@name='users']/object) > 0">
            <div class="span_12_box_masonry">
              <div class="box_header span-12 last rounded-top">
                <div class="box_title"><h3 class="nontyped span-7 last">Users</h3></div>
                <div class="box_button">
                </div>
              </div>

              <div class="box_content rounded-bottom span-12 last">
                <table class="name-value">
                  <tbody>
                    <xsl:for-each select="aspect[@name='users']/object">
                      <tr>
                        <td class="name">
                          <xsl:value-of select="aspect[@name='label']/value"/>
                        </td>
                        <td class="value">
                          <xsl:for-each select="aspect[@name='user']/object">
                            <xsl:call-template name="object_link_button">
                              <xsl:with-param name="linktext" select="./display_name"/>
                              <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                            </xsl:call-template>
                          </xsl:for-each>
                        </td>
                      </tr>
                    </xsl:for-each>
                  </tbody>
                </table>
              </div>
            </div>
          </xsl:if>


        </div> <!-- end .masonry -->
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>
</xsl:stylesheet>
