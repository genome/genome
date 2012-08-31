<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- full page display for a instrumentdata flowcell -->
  <xsl:template name="genome_instrumentdata_solexa" match="object[./types[./isa[@type='Genome::InstrumentData::Solexa']]]">
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Instrument Data Solexa:'" />
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

          <xsl:for-each select="aspect[@name='taxon']/object">
            <xsl:call-template name="genome_taxon_box"/>
          </xsl:for-each>

          <xsl:for-each select="aspect[@name='library']/object">
            <xsl:call-template name="genome_library_box"/>
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

            <tr>
              <td class="name">Subset Name:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='subset_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Project Name:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='project_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Target Region Set Name:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='target_region_set_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Index Sequence:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='index_sequence']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Read Length:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='read_length']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Forward Read Length:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='fwd_read_length']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Reverse Read Length:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='rev_read_length']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Clusters:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='clusters']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Forward Clusters:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='fwd_clusters']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Reverse Clusters:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='rev_clusters']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">Flow Cell ID:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='flow_cell']/object/aspect[@name='flow_cell_id']/value"/>
              </td>
            </tr>

          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>

</xsl:stylesheet>
