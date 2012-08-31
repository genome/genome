<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_druggablegene_drugnamereport_set" match="object[./types[./isa[@type='Genome::DruggableGene::DrugNameReport::Set']]]">

    <script type='text/javascript' src='/res/js/pkg/boxy/javascripts/jquery.boxy.js'></script>
    <link rel="stylesheet" href="/res/js/pkg/boxy/stylesheets/boxy.css" type="text/css" />
    <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:variable name='header_name'>
        <xsl:choose>
            <xsl:when test='count(aspect[@name="name"]/value)=1'>
                <xsl:value-of select='aspect[@name="members"]/object/aspect[@name="human_readable_name"]/value'/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select='./display_name'/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Drug'" />
      <xsl:with-param name="display_name" select='$header_name'/>
      <xsl:with-param name="icon" select="'genome_drugname_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
          <xsl:for-each select="aspect[@name='members']/object">
              <xsl:call-template name="genome_drugnamereport_box"/>
          </xsl:for-each>
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="genome_drugnamereport_box">

    <div class="span_12_box_masonry">
      <div class="box_header span-12 last rounded-top">
        <div class="box_title"><h3 class="genome_drugnamereport_16 span-7 last">
            <xsl:value-of select="normalize-space(aspect[@name='source_db_name']/value)"/>
        </h3></div>
      </div>

      <div class="box_content rounded-bottom span-12 last">
        <table class="name-value">
          <tbody>

            <tr>
              <td class="name">Source Name and Link:
              </td>
              <td class="value">
                  <a class='mini btn'>
                      <xsl:attribute name="href">
                          <xsl:value-of select="normalize-space(aspect[@name='original_data_source_url']/value)"/>
                      </xsl:attribute>
                      <xsl:value-of select="normalize-space(aspect[@name='name']/value)"/>
                  </a>
              </td>
            </tr>

            <tr>
              <td class="name">Source Database Name:
              </td>
              <td class="value">
                  <a class='mini btn'>
                      <xsl:attribute name="href">
                          <xsl:value-of select="normalize-space(aspect[@name='source_db_url']/value)"/>
                      </xsl:attribute>
                  <xsl:value-of select="normalize-space(aspect[@name='source_db_name']/value)"/>
                  </a>
              </td>
            </tr>

            <tr>
              <td class="name">Source Database Version:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='source_db_version']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">Alternate Names:
              </td>
              <td class="value">
                <ul>
                  <xsl:for-each select="aspect[@name='drug_alt_names']/object">
                    <li>
                      <xsl:value-of select="normalize-space(aspect[@name='alternate_name']/value)"/>
                      <xsl:text>  </xsl:text>
                      (<xsl:value-of select="normalize-space(aspect[@name='nomenclature']/value)"/>)
                    </li>
                  </xsl:for-each>
                </ul>
              </td>
            </tr>

            <tr>
              <td class="name">Categories:
              </td>
              <td class="value">
                <ul>
                  <xsl:for-each select="aspect[@name='drug_categories']/object">
                    <li>
                      <xsl:value-of select="normalize-space(aspect[@name='category_name']/value)"/>
                      <xsl:text> : </xsl:text>
                      <xsl:value-of select="normalize-space(aspect[@name='category_value']/value)"/>
                    </li>
                  </xsl:for-each>
                </ul>
              </td>
            </tr>

            <tr>
              <td class="name">Citation:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='citation']/object/aspect/value)"/>
              </td>
            </tr>

          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>

</xsl:stylesheet>
