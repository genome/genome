<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_featurelist" match="object[./types[./isa[@type='Genome::FeatureList']]]">
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Feature List:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_capture_set_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

        <!-- details for this library -->
        <div class="span_8_box_masonry">
          <div class="box_header span-12 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-7 last">Feature List Attributes</h3></div>
            <div class="box_button">

            </div>
          </div>

          <div class="box_content rounded-bottom span-12 last">
            <table class="name-value">
              <tbody>
                <tr>
                  <td class="name">Name:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='name']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Source:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='source']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Format:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='format']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">MD5:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='file_content_hash']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">File Path:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='file_path']/value"/>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        </div><!-- end .objects -->
      </div> <!-- end .container -->
    </div> <!-- end .content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <!-- box element for library, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_featurelist_box">
    <xsl:comment>template: status/genome_featurelist.xsl:genome_featurelist_box</xsl:comment>
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_capture_set_16 span-7 last">Feature List</h3></div>
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
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>


</xsl:stylesheet>
