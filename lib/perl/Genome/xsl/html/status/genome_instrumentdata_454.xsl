<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- full page display for a instrumentdata flowcell -->
  <xsl:template name="genome_instrumentdata_454" match="object[./types[./isa[@type='Genome::InstrumentData::454']]]">
    <xsl:comment>template: /html/status/genome_instrumentdata_454.xsl match="object[./types[./isa[@type='Genome::InstrumentData::454']]]"</xsl:comment>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Instrument Data 454:'" />
      <xsl:with-param name="display_name" select="@id" />
      <xsl:with-param name="icon" select="'genome_instrumentdata_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <xsl:call-template name="genome_instrumentdata_solexa_box"/>

          <xsl:for-each select="aspect[@name='samples']/object">
            <xsl:call-template name="genome_sample_box"/>
          </xsl:for-each>

<!--          <xsl:for-each select="aspect[@name='taxon']/object">
            <xsl:call-template name="genome_taxon_box"/>
          </xsl:for-each>
-->

          <div class="span_8_box_masonry">
            <div class="box_header span-8 last rounded-top">
              <div class="box_title"><h3 class="genome_library_16 span-7 last">Library</h3></div>
              <div class="box_button">
                <xsl:for-each select="aspect[@name='library']/object">
                  <xsl:call-template name="object_link_button_tiny">
                    <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
                  </xsl:call-template>
                </xsl:for-each>
              </div>
            </div>

            <div class="box_content rounded-bottom span-8 last">
              <table class="name-value">
                <tbody>
                  <tr>
                    <td class="name">Name:
                    </td>
                    <td class="value"><xsl:value-of select="aspect[@name='library_name']/value"/>
                    </td>
                  </tr>
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

  <xsl:template name="genome_instrumentdata_solexa_box">

    <xsl:comment>template: genome_model.xsl:genome_model_attributes_box</xsl:comment>

    <!-- details for this solexa instrumentdata -->
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="nontyped span-7 last">Instrument Data Attributes</h3></div>
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
          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>

  <!-- box element for library, intended for display in a jquery masonry layout -->
  <!-- copied from genome_library.xml because the 454 report differs slightly -->
  <xsl:template name="genome_library_box_454">
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
              <td class="value"><xsl:value-of select="../aspect[@name='library_name']/value"/>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>


</xsl:stylesheet>
