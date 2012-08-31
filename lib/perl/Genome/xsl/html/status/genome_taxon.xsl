<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- full page display for a taxon -->
  <xsl:template name="genome_taxon" match="object[./types[./isa[@type='Genome::Taxon']]]">
    <xsl:comment>template: status/genome_taxon.xsl match: object[./types[./isa[@type='Genome::Taxon']]]</xsl:comment>
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Taxon:'" />
      <xsl:with-param name="display_name" select="aspect[@name='species_name']/value" />
      <xsl:with-param name="icon" select="'genome_taxon_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this taxon -->
          <div class="span_8_box_masonry">
            <div class="box_header span-8 last rounded-top">
              <div class="box_title"><h3 class="nontyped span-7 last">Taxon Attributes</h3></div>
              <div class="box_button">

              </div>
            </div>

            <div class="box_content rounded-bottom span-8 last">
              <table class="name-value">
                <tbody>
                  <tr>
                    <td class="name">Name:
                    </td>
                    <td class="value"><xsl:value-of select="aspect[@name='species_name']/value"/>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Latin Name:
                    </td>
                    <td class="value"><xsl:value-of select="aspect[@name='species_latin_name']/value"/>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">NCBI Taxon ID:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='ncbi_taxon_id']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='ncbi_taxon_id']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Estimated Genome Size:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='estimated_genome_size']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='estimated_genome_size']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Gram Stain Category:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='gram_stain_category']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='gram_stain_category']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Domain:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='domain']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='domain']/value)"/>
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

          <xsl:for-each select="aspect[@name='individuals']/object">
            <xsl:call-template name="genome_individual_box"/>
          </xsl:for-each>

        </div> <!-- end .objects -->
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <!-- box element for taxon, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_taxon_box">

    <xsl:comment>template: status/genome_taxon.xsl:genome_taxon_box</xsl:comment>

    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_taxon_16 span-7 last">Taxon</h3></div>
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
              <td class="name">Species Name:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='species_name']/value"/>
              </td>
            </tr>


            <tr>
              <td class="name">Species Latin Name:
              </td>
              <td class="value">
                <xsl:choose>
                  <xsl:when test="string(normalize-space(aspect[@name='species_latin_name']/value))">
                    <xsl:value-of select="normalize-space(aspect[@name='species_latin_name']/value)"/>
                  </xsl:when>
                  <xsl:otherwise>
                    --
                  </xsl:otherwise>
                </xsl:choose>
              </td>
            </tr>

            <tr>
              <td class="name">NCBI Taxon ID:
              </td>
              <td class="value">
                <xsl:choose>
                  <xsl:when test="string(normalize-space(aspect[@name='ncbi_taxon_id']/value))">
                    <xsl:value-of select="normalize-space(aspect[@name='ncbi_taxon_id']/value)"/>
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
