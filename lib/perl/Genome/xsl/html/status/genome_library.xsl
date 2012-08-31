<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_library" match="object[./types[./isa[@type='Genome::Library']]]">
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Library:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_library_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

        <!-- details for this library -->
        <div class="span_8_box_masonry">
          <div class="box_header span-8 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-7 last">Library Attributes</h3></div>
            <div class="box_button">

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

        <!-- taxon objects -->
        <xsl:for-each select="aspect[@name='taxon']/object">
          <xsl:call-template name="genome_taxon_box"/>
        </xsl:for-each>

        <!-- sample objects -->
        <xsl:for-each select="aspect[@name='sample']/object">
          <xsl:call-template name="genome_sample_box"/>
        </xsl:for-each>

        </div><!-- end .objects -->

        <xsl:for-each select="aspect[@name='models']/object[./types[./isa[@type='Genome::Model']]]">
          <xsl:sort select="aspect[@name='is_default']" order="descending"/>
          <xsl:sort select="aspect[@name='name']"/>
          <xsl:call-template name="genome_model_builds_list_table"/>
        </xsl:for-each>

      </div> <!-- end .container -->
    </div> <!-- end .content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <!-- box element for library, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_library_box">
    <xsl:comment>template: status/genome_library.xsl:genome_library_box</xsl:comment>
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_library_16 span-7 last">Library</h3></div>
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
