<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_sample" match="object[./types[./isa[@type='Genome::Sample']]]">
    <xsl:comment>template: /html/status/genome_sample.xsl  match: object[./types[./isa[@type='Genome::Sample']]]</xsl:comment>

    <script type='text/javascript' src='/res/js/app/status/genome_sample.js'></script>
    <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Sample:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_sample_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <!-- details for this sample -->
          <div class="span_8_box_masonry">
            <div class="box_header span-8 last rounded-top">
              <div class="box_title"><h3 class="nontyped span-7 last">Sample Attributes</h3></div>
              <div class="box_button">

              </div>
            </div>

            <div class="box_content rounded-bottom span-8 last">
              <table class="name-value">
                <tbody>

                  <tr>
                    <td class="name">Extraction Label:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='extraction_label']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='extraction_label']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Patient Name:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='individual_common_name']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='individual_common_name']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Extraction Type:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='extraction_type']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='extraction_type']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Tissue Type:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='tissue_type']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='tissue_type']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Organ Name:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='organ_name']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='organ_name']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>

                  <tr>
                    <td class="name">Extraction Type:
                    </td>
                    <td class="value">
                      <xsl:choose>
                        <xsl:when test="string(normalize-space(aspect[@name='extraction_type']/value))">
                          <xsl:value-of select="normalize-space(aspect[@name='extraction_type']/value)"/>
                        </xsl:when>
                        <xsl:otherwise>
                          --
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>
                  <xsl:if test="count(aspect[@name='models']/object/aspect[@name='region_of_interest_set_name']/value) > 0">
                    <tr>
                      <td class="name">Coverage Report:
                      </td>
                      <td class="value">

                        <!-- copied from object_button_link template in order to use non-id parameters -->
                        <a class="mini btn">
                          <xsl:attribute name="href">
                            <xsl:value-of select="$rest"/>
                            <xsl:text>/</xsl:text>
                            <xsl:value-of select="rest:typetourl('Genome::Model::Set')"/>
                            <xsl:text>/</xsl:text>
                            <xsl:value-of select="'coverage'"/>
                            <xsl:text>.</xsl:text>
                            <xsl:value-of select="'html'"/>
                            <xsl:text>?</xsl:text>
                            <xsl:for-each select="aspect[@name='models']/object[aspect[@name='region_of_interest_set_name']/value]">
                              <xsl:text>genome_model_id=</xsl:text>
                              <xsl:value-of select='@id'/>
                              <xsl:text>&amp;</xsl:text>
                            </xsl:for-each>
                          </xsl:attribute>
                          <span class="sm-icon sm-icon-extlink"><br/></span><xsl:text>Coverage Report</xsl:text>
                        </a>

                      </td>
                    </tr>
                  </xsl:if>
                </tbody>
              </table>
            </div>
          </div>

          <xsl:for-each select="aspect[@name='libraries']/object">
            <xsl:call-template name="genome_library_box"/>
          </xsl:for-each>

          <xsl:for-each select="aspect[@name='taxon']/object">
            <xsl:call-template name="genome_taxon_box"/>
          </xsl:for-each>

          <!-- TODO: figure out how to get UR to load the templates for objects contained in the 'source' node so
                     that the templates necessary to render the XSL below will be available. -->
          <!-- <xsl:for-each select="aspect[@name='source']/object"> -->
          <!--   <xsl:choose> -->
          <!--     <xsl:when test="@type = 'Genome::Individual'"> -->
          <!--       <xsl:call-template name="genome_individual_box"/> -->
          <!--     </xsl:when> -->
          <!--     <xsl:when test="@type = 'Genome::PopulationGroup'"> -->
          <!--       <xsl:call-template name="genome_populationgroup_box"/> -->
          <!--     </xsl:when> -->
          <!--   </xsl:choose> -->
          <!-- </xsl:for-each> -->


          <!--<xsl:for-each select="aspect[@name='projects']/object">-->
            <!--<xsl:call-template name="genome_project_box"/>-->
          <!--</xsl:for-each>-->


        </div> <!-- end .masonry -->
        <xsl:for-each select="aspect[@name='models']/object[./types[./isa[@type='Genome::Model']]]">
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


  <!-- box element for sample, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_sample_box">
    <xsl:comment>template: /html/status/genome_sample.xsl; name="genome_sample_box"</xsl:comment>
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_sample_16 span-7 last">Sample</h3></div>
        <div class="box_button">
        </div>
      </div>

      <div class="box_content rounded-bottom span-8 last">
        <table class="name-value">
          <tbody>
            <tr>
              <td class="name">ID:
              </td>
              <td class="value">
                <xsl:for-each select=".">
                  <xsl:call-template name="object_link_button">
                    <xsl:with-param name="linktext" select="@id" />
                    <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                  </xsl:call-template>
                </xsl:for-each>
              </td>
            </tr>

            <tr>
              <td class="name">
                Name:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Patient Common Name:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='individual_common_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Common Name:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='common_name']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Tissue Label:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='tissue_label']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Tissue Description:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='tissue_desc']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Extraction Type:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='extraction_type']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Extraction Label:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='extraction_label']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                Extraction Description:
              </td>
              <td class="value">
                <xsl:choose>
                  <xsl:when test="string(normalize-space(aspect[@name='extraction_desc']/value))">
                    <xsl:value-of select="aspect[@name='extraction_desc']/value"/>
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
