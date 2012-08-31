<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- full page display for a processing profile -->

  <xsl:template name="genome_processingprofile" match="object[./types[./isa[@type='Genome::ProcessingProfile']]]">
    <xsl:comment>template: status/genome_taxon.xsl match: object[./types[./isa[@type='Genome::ProcessingProfile']]]</xsl:comment>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Processing Profile:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_processingprofile_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this taxon -->
          <div class="span_8_box_masonry">
            <div class="box_header span-12 last rounded-top">
              <div class="box_title"><h3 class="nontyped span-7 last">Processing Profile Attributes</h3></div>
              <div class="box_button">

              </div>
            </div>

            <div class="box_content rounded-bottom span-12 last">
              <table class="name-value">
                <tbody>
                  <tr>
                    <td class="name">Type:
                    </td>
                    <td class="value"><xsl:value-of select="aspect[@name='type_name']/value"/>

                    </td>
                  </tr>

                  <xsl:if test="normalize-space(aspect[@name='supersedes']/value)">
                    <tr>
                      <td class="name">Supercedes:
                      </td>
                      <td class="value"><xsl:value-of select="aspect[@name='supersedes']/value"/>
                      </td>
                    </tr>
                  </xsl:if>

                  <xsl:if test="count(aspect[@name='params']) > 0 ">
                    <xsl:for-each select="aspect[@name='params']/object">
                      <tr>
                        <td class="name"><xsl:value-of select="normalize-space(aspect[@name='name'])"/>:</td>
                        <td class="value"><xsl:value-of select="aspect[@name='_value_scalar_or_object']/value"/></td>
                      </tr>
                    </xsl:for-each>
                  </xsl:if>

                </tbody>
              </table>
            </div>
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

  <!-- box element for processing profile, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_processingprofile_box">

    <xsl:comment>template: status/genome_processingprofile.xsl:genome_processingprofile_box</xsl:comment>

    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_processingprofile_16 span-7 last">Processing Profile</h3></div>
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
              <td class="name">Name:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Type:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='type_name']/value"/>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>


</xsl:stylesheet>
