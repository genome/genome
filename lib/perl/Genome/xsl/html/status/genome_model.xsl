<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <!-- full page template for models -->
  <xsl:template name="genome_model" match="object[./types[./isa[@type='Genome::Model']]]">
    <!-- These parameters are only valid when there is a last_complete_build; they are here for overriding by subclasses -->
    <xsl:param name="build_directory_url">
      <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='last_complete_build']/object/aspect[@name='data_directory']/value)" />
    </xsl:param>
    <xsl:param name="summary_report_url">
      <xsl:value-of select="$build_directory_url"/><xsl:text>/reports/Summary/report.html</xsl:text>
    </xsl:param>

    <xsl:comment>template: status/genome_model.xsl match: object[./types[./isa[@type='Genome::Model']]]</xsl:comment>

    <script type='text/javascript' src='/res/js/app/genome_model.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Model:'" />
      <xsl:with-param name="display_name" select="./aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_model_32'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

          <xsl:call-template name="genome_model_attributes_box"/>

          <xsl:for-each select="aspect[@name='processing_profile']/object">
            <xsl:call-template name="genome_processingprofile_box"/>
          </xsl:for-each>

          <xsl:for-each select="aspect[@name='last_complete_build_flagstat']/object[@type='UR::Value::HASH']/perldata/hashref">
            <xsl:call-template name="genome_model_flagstat_table"/>
          </xsl:for-each>


        </div> <!-- end .objects -->

        <xsl:if test="count(aspect[@name='inputs']) > 0 ">
          <xsl:for-each select="aspect[@name='inputs']">
            <xsl:call-template name="genome_model_input_table"/>
          </xsl:for-each>
        </xsl:if>

        <xsl:if test="count(aspect[@name='to_models'] | aspect[@name='from_models']) > 0">
          <xsl:call-template name="genome_model_link_table"/>
        </xsl:if>

        <xsl:if test="count(aspect[@name='builds']) > 0">
          <xsl:call-template name="genome_model_build_lister"/>
        </xsl:if>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <!-- initializes the dataTable plugin for model set views -->
  <xsl:template name="genome_model_set_table_init" match="object[./types[./isa[@type='Genome::Model']]]" mode="set_table_init">
    <xsl:comment>template: status/genome_model.xsl match: object[./types[./isa[@type='Genome::Model']]] mode: set_table_init</xsl:comment>
    <script type="text/javascript">
      <xsl:text disable-output-escaping="yes">
        <![CDATA[
                 $(document).ready(
                 window.setTable = $('#set').dataTable({
                 "sScrollX": "100%",
                 "sScrollInner": "110%",
                 "bJQueryUI": true,
                 "sPaginationType": "full_numbers",
                 "bStateSave": true,
                 "iDisplayLength": 25
                 })
                 );
        ]]>
      </xsl:text>
    </script>
  </xsl:template>

  <!-- describes the columns for model set views -->
  <xsl:template name="genome_model_set_header" match="object[./types[./isa[@type='Genome::Model']]]" mode="set_header">
    <xsl:comment>template: status/genome_model.xsl match: object[./types[./isa[@type='Genome::Model']]] mode: set_header</xsl:comment>
    <tr>
      <th>
        name
      </th>
      <th>
        ID
      </th>
      <th>
        creation date
      </th>
      <th>
        user name
      </th>
      <th>
        builds
      </th>
      <th>
        last complete build
      </th>
    </tr>
  </xsl:template>

  <!-- describes the row for model set views -->
  <xsl:template name="genome_model_set_row" match="object[./types[./isa[@type='Genome::Model']]]" mode="set_row">
    <xsl:comment>template: status/genome_model.xsl match: object[./types[./isa[@type='Genome::Model']]] mode: set_row</xsl:comment>

    <xsl:variable name="total_builds" select="count(aspect[@name='builds']/object)"/>
    <xsl:variable name="failed_builds" select="count(aspect[@name='builds']/object/aspect[@name='status'][value = 'Failed'])"/>
    <tr>
      <td>
        <xsl:value-of select="aspect[@name='name']"/>
      </td>
      <td>
        <xsl:call-template name="object_link_button">
          <xsl:with-param name="linktext" select="aspect[@name='genome_model_id']"/>
          <xsl:with-param name="icon" select="'sm-icon-extlink'" />
        </xsl:call-template>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='creation_date']"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='created_by']"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='run_as']"/>
      </td>
      <td>
        <xsl:value-of select="$total_builds"/>
        <xsl:if test="$failed_builds">
          <span style="color: red;"> (<xsl:value-of select="$failed_builds"/>)</span>
        </xsl:if>
      </td>

      <td>
        <xsl:choose>
          <xsl:when test="aspect[@name='last_complete_build']/object">
            <xsl:for-each select="aspect[@name='last_complete_build']/object">
              <xsl:call-template name="object_link_button">
                <xsl:with-param name="linktext" select="@id" />
                <xsl:with-param name="icon" select="'sm-icon-extlink'" />
              </xsl:call-template>
            </xsl:for-each>
          </xsl:when>
          <xsl:otherwise>
            --
          </xsl:otherwise>
        </xsl:choose>
      </td>

    </tr>
  </xsl:template>

  <xsl:template name="genome_model_build_lister">
    <xsl:comment>template: genome_model.xsl:genome_model_build_lister</xsl:comment>

    <div class="generic_lister">
      <div class="box_header span-24 last rounded-top">
        <div class="box_title"><h3 class="genome_model_build_16 span-24 last">Builds</h3></div>
      </div>
      <div class="box_content rounded-bottom span-24 last">
        <table class="lister">
          <thead>
            <tr>
              <th>build</th>
              <th>status</th>
              <th>scheduled</th>
              <th>completed</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <xsl:for-each select="aspect[@name='builds']/object">
                <xsl:call-template name="genome_model_builds_list_table_row" />
              </xsl:for-each>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>

  <xsl:template name="genome_model_attributes_box">
    <xsl:param name="build_directory_url">
      <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='last_complete_build']/object/aspect[@name='data_directory']/value)" />
    </xsl:param>

    <xsl:param name="summary_report_url">
      <xsl:value-of select="$build_directory_url"/><xsl:text>/reports/Summary/report.html</xsl:text>
    </xsl:param>
    <xsl:param name="metrics_chart_url">
      <xsl:text>/view/genome/model/build/set/chart.html?model_id=</xsl:text><xsl:value-of select="@id"/>
    </xsl:param>
    <xsl:param name="model_attributes_type_name">
      <xsl:value-of select="normalize-space(aspect[@name='processing_profile']/object/aspect[@name='type_name']/value)"/>
    </xsl:param>

    <xsl:comment>template: genome_model.xsl:genome_model_attributes_box</xsl:comment>
    <!-- details for this model -->
    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="nontyped span-7 last">Model Attributes</h3></div>
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
            <!-- notes button -->
            <xsl:variable name="note_count" select="count(aspect[@name='notes']/object)"/>
            <xsl:variable name="model_id" select="@id"/>
            <xsl:choose>
              <xsl:when test="$note_count &gt; 0">

                <tr>
                  <td class="name"></td>
                  <td class="value">
                    <a class="mini btn notes-popup">
                      <xsl:attribute name="title">Model <xsl:value-of select="$model_id"/> Notes</xsl:attribute>
                      <xsl:attribute name="id"><xsl:value-of select="model_id"/></xsl:attribute>
                      <span class="sm-icon sm-icon-newwin"><br/></span>notes (<xsl:value-of select="$note_count"/>)
                    </a>
                    <!-- div for instrument data -->
                    <div style="display: none;">
                      <xsl:attribute name="id">notes_subject_<xsl:value-of select="model_id"/></xsl:attribute>
                      <table class="lister" border="0" width="100%" cellspacing="0" cellpadding="0">
                        <colgroup>

                        </colgroup>
                        <thead>
                          <th>
                            header
                          </th>
                          <th>
                            date
                          </th>
                          <th>
                            editor id
                          </th>
                          <th>
                            <xsl:text> </xsl:text>
                          </th>

                        </thead>
                        <tbody>
                          <xsl:for-each select="aspect[@name='notes']/object">
                            <xsl:sort select="aspect[@name='entry_date_sort']/value" data-type="number" order="ascending"/>
                            <tr>
                              <td><strong><xsl:value-of select="aspect[@name='header_text']/value"/></strong></td>
                              <td><xsl:value-of select="aspect[@name='entry_date']/value"/></td>
                              <td><xsl:value-of select="aspect[@name='editor_id']/value"/></td>
                              <td class="buttons">
                                <xsl:for-each select="aspect[@name='subject']/object">
                                  <xsl:call-template name="object_link_button">
                                    <xsl:with-param name="linktext" select="'subject'"/>
                                    <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
                                  </xsl:call-template>
                                </xsl:for-each>
                              </td>

                            </tr>
                            <tr>
                              <td colspan="4" class="text"><xsl:value-of select="aspect[@name='body_text']/value"/></td>
                            </tr>

                          </xsl:for-each>
                        </tbody>
                      </table>

                    </div>
                  </td>
                </tr>
              </xsl:when>
            </xsl:choose>

            <tr>
              <td class="name">creation date:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='creation_date']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">created by:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='created_by']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">run as:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='run_as']/value"/>
              </td>
            </tr>

            <tr>
              <td class="name">status:
              </td>
              <td class="value"><xsl:value-of select="aspect[@name='status']/value"/>
              </td>
            </tr>

            <xsl:choose>
              <xsl:when test="aspect[@name='last_complete_build']/object">

                <tr>
                  <td class="name">last complete build:
                  </td>
                  <td class="value">
                    <xsl:for-each select="aspect[@name='last_complete_build']/object">
                      <xsl:call-template name="object_link_button">

                        <xsl:with-param name="linktext" select="@id" />
                        <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                      </xsl:call-template>
                    </xsl:for-each>
                  </td>
                </tr>

                <tr>
                  <td class="name"><br/>
                  </td>
                  <td class="value">

                    <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$build_directory_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>data directory</a>

                  </td>
                </tr>

                <tr>
                  <td class="name"><br/>
                  </td>
                  <td class="value">

                    <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$summary_report_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>summary report</a>

                  </td>
                </tr>
                <tr>
                  <td class="name"><br/>
                  </td>
                  <td class="value">
                    <xsl:for-each select="aspect[@name='last_complete_build_flagstat']/object[@type='UR::Value::HASH']/perldata/hashref">
                      <a class="mini btn" id="flagstat_button" href="#"><span class="sm-icon sm-icon-newwin"><br/></span>flagstat report</a>
                    </xsl:for-each>

                  </td>
                </tr>

              </xsl:when>
              <xsl:otherwise>
                <tr>
                  <td class="name">last complete build:
                  </td>
                  <td class="value">
                    --
                  </td>
                </tr>
              </xsl:otherwise>
            </xsl:choose>

            <xsl:if test="$model_attributes_type_name = 'benchmark'">
              <tr>
                <td class="name">metrics:
                </td>
                <td class="value"><a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$metrics_chart_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>metrics</a>
                </td>
              </tr>
            </xsl:if>

            <!--            <tr>
                 <td class="name">subject type:
                 </td>
                 <td class="value"><xsl:value-of select="normalize-space(aspect[@name='subject_class_name']/value)"/>
                 </td>
                 </tr>
            -->
            <tr>
              <td class="name">subject:
              </td>
              <td class="value">

                <xsl:call-template name="object_link_button">
                  <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                  <xsl:with-param name="type" select="normalize-space(/object/aspect[@name='subject']/object/aspect[@name='subclass_name']/value)"/>
                  <xsl:with-param name="id" select="normalize-space(/object/aspect[@name='subject']/object/aspect[@name='subject_id']/value)"/>
                  <xsl:with-param name="linktext">
                    <xsl:choose>
                      <xsl:when test="string(/object/aspect[@name='subject']/object/aspect[@name='name']/value)">
                        <xsl:value-of select="normalize-space(/object/aspect[@name='subject']/object/aspect[@name='name']/value)"/>
                      </xsl:when>
                      <xsl:otherwise>
                        <xsl:value-of select="normalize-space(/object/aspect[@name='subject']/object/aspect[@name='subject_id']/value)"/>
                      </xsl:otherwise>
                    </xsl:choose>
                  </xsl:with-param>
                </xsl:call-template>
                <!-- <xsl:if test="not(string(aspect[@name='subject_name']/value))"> -->
                <!--   <br/><span class="small">(subject_name unspecified)</span> -->
                <!-- </xsl:if> -->
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>

  <!-- box element for convergence model, intended for display in a jquery masonry layout -->
  <xsl:template name="genome_model_convergence_box" match="object[./types[./isa[@type='Genome::Model::Convergence']]]" mode="box">
    <xsl:param name="last_complete_build_directory_url">
      <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="normalize-space(aspect[@name='last_complete_build']/object/aspect[@name='data_directory']/value)" />
    </xsl:param>

    <xsl:param name="last_complete_build_summary_report_url">
      <xsl:value-of select="$last_complete_build_directory_url"/><xsl:text>/reports/Summary/report.html</xsl:text>
    </xsl:param>

    <xsl:comment>template: status/genome_model.xsl:genome_model_box; match: object[./types[./isa[@type='Genome::Model']]]; mode: box</xsl:comment>

    <div class="span_8_box_masonry">
      <div class="box_header span-8 last rounded-top">
        <div class="box_title"><h3 class="genome_model_16 span-7 last">Convergence Model</h3></div>
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
              <td class="name">ID:
              </td>
              <td class="value"><xsl:value-of select="@id"/>
              </td>
            </tr>

            <tr>
              <td class="name">Name:
              </td>
              <td class="value">
                <xsl:value-of select="normalize-space(aspect[@name='name']/value)"/>
              </td>
            </tr>

            <xsl:choose>
              <xsl:when test="aspect[@name='last_complete_build']/object">

                <tr>
                  <td class="name">Last Complete Build:
                  </td>
                  <td class="value">
                    <xsl:for-each select="aspect[@name='last_complete_build']/object">
                      <xsl:call-template name="object_link_button">

                        <xsl:with-param name="linktext" select="@id" />
                        <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                      </xsl:call-template>
                    </xsl:for-each>
                  </td>
                </tr>

                <tr>
                  <td class="name"><br/>
                  </td>
                  <td class="value">

                    <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$last_complete_build_directory_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>data directory</a>

                  </td>
                </tr>

                <tr>
                  <td class="name"><br/>
                  </td>
                  <td class="value">

                    <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$last_complete_build_summary_report_url'/></xsl:attribute><span class="sm-icon sm-icon-extlink"><br/></span>summary report</a>

                  </td>
                </tr>

              </xsl:when>
              <xsl:otherwise>
                <tr>
                  <td class="name">Last Complete Build:
                  </td>
                  <td class="value">
                    --
                  </td>
                </tr>
              </xsl:otherwise>
            </xsl:choose>
          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>



  <xsl:template name="genome_model_input_table">
    <xsl:comment>template: status/genome_model.xsl:genome_model_input_table</xsl:comment>
    <div class="generic_lister">
      <div class="box_header span-24 last rounded-top">
        <div class="box_title"><h3 class="genome_instrumentdata_16 span-24 last">Inputs</h3></div>
      </div>
      <div class="box_content rounded-bottom span-24 last">
        <table class="lister">
          <tbody>
            <xsl:for-each select="object[aspect[@name='name']/value!='instrument_data']">
              <xsl:sort select="aspect[@name='name']/value" data-type="text" order="ascending"/>
              <xsl:call-template name="genome_model_input_table_row"/>
            </xsl:for-each>

            <!-- It has been decided that instrument_data should be last -->
            <xsl:for-each select="object[aspect[@name='name']/value='instrument_data']">
              <xsl:sort select="aspect[@name='name']/value" data-type="text" order="ascending"/>
              <xsl:call-template name="genome_model_input_table_row"/>
            </xsl:for-each>

          </tbody>
        </table>
      </div> <!-- end box_content -->
    </div> <!-- end generic lister -->

  </xsl:template>

  <xsl:template name="genome_model_input_table_row">
    <xsl:comment>template: status/genome_model.xsl:genome_model_input_table_row</xsl:comment>
    <tr>

      <td><xsl:value-of select="normalize-space(aspect[@name='name']/value)"/></td>

      <td style="width: 100%">
        <xsl:call-template name="genome_model_input_value_link"/>
      </td>

    </tr>
  </xsl:template>

  <xsl:template name="genome_model_input_value_link">
    <xsl:call-template name="object_link_button">
      <xsl:with-param name="icon" select="'sm-icon-extlink'" />
      <xsl:with-param name="id">
        <xsl:value-of select="aspect[@name='value_id']/value"/>
      </xsl:with-param>

      <xsl:with-param name="type">
        <xsl:value-of select="aspect[@name='value_class_name']/value"/>
      </xsl:with-param>

      <xsl:with-param name="linktext">
        <xsl:choose>
          <xsl:when test="aspect[@name='value']/object/display_name">
            <xsl:value-of select="aspect[@name='value']/object/display_name"/>
          </xsl:when>
          <xsl:when test="display_name">
            <xsl:value-of select="display_name"/>
          </xsl:when>
          <xsl:otherwise>
            could not resolve display_name
          </xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>
  </xsl:template>

  <xsl:template name="genome_model_input_rows">
    <xsl:param name="report_missing"/>
    <xsl:param name="input"/>

    <xsl:comment>template: status/genome_model.xsl name: genome_model_input_rows</xsl:comment>

    <xsl:variable name="inputs" select="$input/object"/>
    <!-- select unique input names -->
    <xsl:for-each select="$input">
      <xsl:for-each select="object[not(preceding-sibling::*/aspect[@name='name']/value = aspect[@name='name']/value)]/aspect[@name='name']/value">
        <xsl:variable name="name" select="."/>
        <xsl:variable name="model_id" select="ancestor::object[types/isa[@type='Genome::Model']]/@id"/>
        <xsl:variable name="model_name" select="ancestor::object[types/isa[@type='Genome::Model']]/aspect[@name='name']/value"/>

        <tr>
          <xsl:choose>

            <xsl:when test="count($inputs[aspect[@name='name']/value = $name]) = 1">
              <!-- only one instrument data difference, just show the button -->
              <td><xsl:for-each select="$inputs[aspect[@name='name']/value = $name]"><xsl:call-template name="genome_model_input_value_link"/></xsl:for-each></td>
              <td><xsl:value-of select="$name"/></td>
            </xsl:when>

            <xsl:otherwise>
              <!-- more than one instrument data difference, create button + popup to show all -->
              <td>
                <a href="javascript:return(false)" class="mini btn popup-ajax-instrument-data">
                  <xsl:attribute name="title"><xsl:value-of select="$model_name"/> (<xsl:value-of select="$model_id"/>) instrument_data</xsl:attribute>
                  <xsl:attribute name="id"><xsl:value-of select="$model_id"/></xsl:attribute>
                  <span class="sm-icon sm-icon-newwin"><br/></span>instrument data (<xsl:value-of select="count($inputs[aspect[@name='name']/value = $name])"/>)
                </a>
                <!-- div for instrument data -->
                <div style="display: none;">
                  <xsl:attribute name="id">instrument_data_<xsl:value-of select="$model_id"/></xsl:attribute>
                  <table class="lister">
                    <xsl:for-each select="$inputs[aspect[@name='name']/value = $name]">
                      <xsl:call-template name="genome_model_input_table_row"/>
                    </xsl:for-each>
                  </table>
                </div>
              </td>
              <td><xsl:value-of select="$name"/></td>
            </xsl:otherwise>
          </xsl:choose>
        </tr>
      </xsl:for-each>
    </xsl:for-each>

    <!-- select unique input names (again)-->
    <xsl:if test="$report_missing">
      <xsl:variable name="missing_inputs" select="$report_missing/object"/>
      <xsl:for-each select="$missing_inputs[not(preceding-sibling::*/aspect[@name='name']/value = aspect[@name='name']/value)]/aspect[@name='name']/value">
        <tr>
          <td>
            <xsl:variable name="name" select="."/>
            <xsl:choose>
              <xsl:when test="count($missing_inputs[aspect[@name='name']/value = $name]) = 1">
                <span class="label"><xsl:value-of select="$name"/></span>:
                <span class="value"><em><xsl:text>missing</xsl:text></em></span>
              </xsl:when>
              <xsl:otherwise>
                <span class="label"><xsl:value-of select="$name"/></span>:
                <span class="value"><em><xsl:text>missing </xsl:text><xsl:value-of select="count($missing_inputs[aspect[@name='name']/value = $name])"/></em></span>
              </xsl:otherwise>
            </xsl:choose>
          </td>
        </tr>
      </xsl:for-each>
    </xsl:if>

  </xsl:template>


  <xsl:template name="genome_model_input_compact_list">
    <xsl:param name="input"/>
    <xsl:param name="report_missing"/>

    <xsl:variable name="inputs" select="$input/object"/>

    <!-- select unique input names -->
    <xsl:for-each select="$input">
      <xsl:variable name="build_id" select="ancestor::object[types/isa[@type='Genome::Model::Build']]/@id"/>
      <xsl:variable name="build_name" select="ancestor::object[types/isa[@type='Genome::Model::Build']]/aspect[@name='display_name']/value"/>

      <xsl:for-each select="object[not(preceding-sibling::*/aspect[@name='name']/value = aspect[@name='name']/value)]/aspect[@name='name']/value">
        <xsl:variable name="name" select="."/>
        <span class="input_diff">
          <xsl:choose>
            <!-- only one instrument data difference, just show the button -->
            <xsl:when test="count($inputs[aspect[@name='name']/value = $name]) = 1">
              <span class="label"><xsl:value-of select="$name"/></span>:
              <span class="value"><xsl:for-each select="$inputs[aspect[@name='name']/value = $name]"><xsl:call-template name="genome_model_input_value_link"/></xsl:for-each></span>;
            </xsl:when>

            <xsl:otherwise>
              <!-- more than one instrument data difference, create button + popup to show all -->
              <a href="javascript:return(false)" class="mini btn popup-ajax-input-diff">
                <xsl:attribute name="title">(<xsl:value-of select="$build_id"/>) input differences</xsl:attribute>
                <xsl:attribute name="id"><xsl:value-of select="$build_id"/></xsl:attribute>
                <span class="sm-icon sm-icon-newwin"><br/></span><xsl:value-of select="$name"/> (<xsl:value-of select="count($inputs[aspect[@name='name']/value = $name])"/>)
              </a>

              <!-- div for instrument data -->
              <div style="display: none;">
                <xsl:attribute name="id">input_diff_<xsl:value-of select="$build_id"/></xsl:attribute>
                <table class="lister">
                  <xsl:for-each select="$inputs[aspect[@name='name']/value = $name]">
                    <xsl:call-template name="genome_model_input_table_row"/>
                  </xsl:for-each>
                </table>
              </div>

            </xsl:otherwise>

          </xsl:choose>
        </span>
      </xsl:for-each>
    </xsl:for-each>

    <!-- select unique input names (again)-->
    <xsl:if test="$report_missing">
      <xsl:variable name="missing_inputs" select="$report_missing/object"/>
      <xsl:for-each select="$missing_inputs[not(preceding-sibling::*/aspect[@name='name']/value = aspect[@name='name']/value)]/aspect[@name='name']/value">
        <span class="input_diff">
          <xsl:variable name="name" select="."/>
          <xsl:choose>

            <xsl:when test="count($missing_inputs[aspect[@name='name']/value = $name]) = 1">
              <span class="label"><xsl:value-of select="$name"/></span>:
              <span class="value"><xsl:text>missing</xsl:text></span>
            </xsl:when>

            <xsl:otherwise>
              <span class="label"><xsl:value-of select="$name"/></span>:
              <span class="value"><xsl:text>missing </xsl:text><xsl:value-of select="count($missing_inputs[aspect[@name='name']/value = $name])"/></span>
            </xsl:otherwise>

            </xsl:choose>;
        </span>
      </xsl:for-each>
    </xsl:if>
  </xsl:template>

  <xsl:template name="genome_model_link_table">
    <xsl:comment>template: status/genome_model.xsl:genome_model_link_table</xsl:comment>
    <div class="generic_lister">
      <div class="box_header span-24 last rounded-top">
        <div class="box_title"><h3 class="genome_model_16 span-24 last">Model Links</h3></div>
      </div>
      <div class="box_content rounded-bottom span-24 last">
        <table class="lister">
          <thead>
            <tr>
              <th>direction</th>
              <th>model</th>
              <th><br/></th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="aspect[@name='to_models']/object | aspect[@name='to_builds']/object">
              <xsl:call-template name="genome_model_link_table_row">
                <xsl:with-param name="type">to</xsl:with-param>
              </xsl:call-template>
            </xsl:for-each>
            <xsl:for-each select="aspect[@name='from_models']/object | aspect[@name='from_builds']/object">
              <xsl:call-template name="genome_model_link_table_row">
                <xsl:with-param name="type">from</xsl:with-param>
              </xsl:call-template>
            </xsl:for-each>
          </tbody>
        </table>
      </div>
    </div>

  </xsl:template>

  <xsl:template name="genome_model_link_table_row">
    <xsl:param name="type" select="''" />

    <xsl:comment>template: status/genome_model.xsl:genome_model_link_table_table_row</xsl:comment>

    <tr>
      <td><xsl:value-of select="$type"/></td>
      <td>
        <xsl:value-of select="@type"/><xsl:text>: </xsl:text>
        <xsl:choose>
          <xsl:when test="aspect[@name='name']/value">
            <xsl:value-of select="aspect[@name='name']/value"/> (#<xsl:value-of select="@id"/>)
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="@id"/>
          </xsl:otherwise>
        </xsl:choose>
      </td>
      <td class="buttons">
        <xsl:call-template name="object_link_button_tiny">
          <xsl:with-param name="icon" select="'sm-icon-extlink'" />
        </xsl:call-template>
      </td>
    </tr>
  </xsl:template>

  <!-- genome model w/ build list template -->

  <xsl:template name="genome_model_builds_list_table">
    <xsl:variable name="is_default" select="aspect[@name='is_default']/value" />
    <xsl:variable name="build_requested" select="aspect[@name='build_requested']/value" />
    <xsl:variable name="build_needed" select="aspect[@name='build_needed']/value" />

    <xsl:comment>template: status/genome_model.xsl; name: genome_model_builds_list_table</xsl:comment>
    <div class="model_builds_list">
      <div class="box_header span-24 last rounded-top">
        <div class="box_title">
          <table cellpadding="0" cellspacing="0" border="0" class="name-value-row">
            <tr>
              <td class="model_name" colspan="6">
                <h3>
                  <xsl:attribute name="class">
                    <xsl:choose>
                      <xsl:when test="$is_default = 1">genome_model_default_16</xsl:when>
                      <xsl:otherwise>genome_model_16</xsl:otherwise>
                    </xsl:choose>
                    <xsl:text> </xsl:text>
                    <xsl:if test="$build_requested = 1">genome_model_build_requested</xsl:if>
                    <xsl:text> </xsl:text>
                    <xsl:if test="$build_needed = 1">genome_model_build_needed</xsl:if>
                  </xsl:attribute>

                  <xsl:value-of select="aspect[@name='name']/value"/>

                </h3>
              </td>
            </tr>
          </table>
        </div>
        <div class="box_button">
          <xsl:call-template name="object_link_button_tiny">
            <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
          </xsl:call-template>
        </div>
      </div> <!-- end box header -->
      <xsl:variable name="num_builds" select="count(aspect[@name='builds']/object)"/>

      <div class="box_header_details rounded-bottom-left">
        <div class="rounded" style="float: left; width: 29%;">
          <table class="name-value" cellpadding="0" cellspacing="0" border="0" style="margin-left: 10px;">
            <tr>
              <td class="name">
                model id:
              </td>
              <td class="value">
                <xsl:value-of select="@id"/>
              </td>
            </tr>
            <tr>
              <td class="name">
                subject name:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='subject']/object/display_name"/>
              </td>
            </tr>
            <tr>
              <td class="name">
                subject type:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='subject']/object/label_name"/>
              </td>
            </tr>
            <tr>
              <td class="name">
                run as:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='run_as']/value"/>
              </td>
            </tr>
            <tr>
              <td class="name">
                created by:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='created_by']/value"/>
              </td>
            </tr>
            <tr>

              <td class="name">
                created:
              </td>
              <td class="value">
                <xsl:value-of select="aspect[@name='creation_date']/value"/>
              </td>
            </tr>
            <xsl:if test="count(aspect[@name='_sample_subject']/object) > 0">
              <tr>
                <td class="name" colspan="2" style="text-align: left">
                  subject sample:
                </td>
              </tr>
              <tr>
                <td class="value" colspan="2">
                  <xsl:for-each select="aspect[@name='_sample_subject']/object">
                    <xsl:call-template name="object_link_button">
                      <xsl:with-param name="linktext" select="display_name" />
                      <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                    </xsl:call-template>
                  </xsl:for-each>
                </td>
              </tr>
            </xsl:if>

            <xsl:if test="normalize-space(aspect[@name='build_requested']/value) != '' and normalize-space(aspect[@name='build_requested']/value) != '0'">
              <tr>

                <td class="name">
                  build requested:
                </td>
                <td class="value">
                  yes
                </td>
              </tr>
            </xsl:if>
          </table>
        </div>

        <div style="float: right; width: 67%; margin-right: 27px;">
          <table class="name-value" cellpadding="0" cellspacing="0" border="0">
            <tbody>
              <colgroup>
                <col />
                <col width="100%"/>
              </colgroup>
              <tr>
                <td colspan="2" style="border-bottom: 1px solid #acaca3; padding-top: -1px;">
                  <span style="font-weight: bold;">
                    processing profile:
                  </span>
                  <xsl:call-template name="object_link_button">
                    <xsl:with-param name="linktext" select="aspect[@name='processing_profile']/object/@id" />
                    <xsl:with-param name="type" select="aspect[@name='processing_profile']/object/@type" />
                    <xsl:with-param name="id" select="aspect[@name='processing_profile']/object/@id" />
                    <xsl:with-param name="icon" select="'sm-icon-extlink'" />
                  </xsl:call-template>
                  <xsl:text> </xsl:text>
                  <xsl:value-of select="aspect[@name='processing_profile']/object/aspect[@name='name']/value"/>
                </td>
              </tr>
              <tr>
                <td colspan="2" style="border-bottom: 1px solid #acaca3; padding-top: -1px;">
                  <span style="font-weight: bold;">
                    inputs
                    <xsl:if test="normalize-space(aspect[@name='build_needed']/value) != ''">
                      <span style="color: #a40014;">(build needed)</span>
                    </xsl:if>
                  </span>
                </td>
              </tr>
              <xsl:choose>
                <xsl:when test="count(aspect[@name='inputs']) > 0 ">
                  <xsl:for-each select="aspect[@name='inputs']">
                    <xsl:call-template name="genome_model_input_rows">
                      <xsl:with-param name="input" select="."/>
                    </xsl:call-template>
                  </xsl:for-each>
                </xsl:when>
                <xsl:otherwise>
                  <tr><td><p>No model inputs found.</p></td></tr>
                </xsl:otherwise>
              </xsl:choose>
            </tbody>
          </table>

        </div>

      </div>

      <div class="span-23 prepend-1 last">
        <div style="margin-left: -5px; background-color: #e5e5e3; width: 5px; height:5px; float: left;">
          <div class="rounded-top-right" style="background-color: #f7f7f7; width: 5px; height: 5px;">
            <br/>
          </div>
        </div>
        <div class="box_content rounded-bottom span-23 last" style="margin-bottom: 20px;">

          <xsl:if test="$num_builds &gt; 0">
            <table cellpadding="0" cellspacing="0" border="0" class="lister" style="margin-bottom: 8px;">
              <col width="10%" />
              <col width="10%" />
              <col/>
              <col/>
              <col/>
              <thead>
                <tr>
                  <th>build</th>
                  <th class="center">status</th>
                  <th>scheduled</th>
                  <th>completed</th>
                  <th class="rounded-right"><br/></th>
                </tr>
              </thead>
              <tbody>
                <xsl:choose>
                  <xsl:when test="count(aspect[@name='builds']/object) > 0">
                    <xsl:for-each select="aspect[@name='builds']/object">
                      <xsl:call-template name="genome_model_builds_list_table_row" />
                    </xsl:for-each>
                  </xsl:when>
                  <xsl:when test="count(aspect[@name='last_succeeded_build']/object) > 0" >
                    <xsl:for-each select="aspect[@name='last_succeeded_build']/object">
                      <xsl:call-template name="genome_model_builds_list_table_row" />
                    </xsl:for-each>
                  </xsl:when>
                  <xsl:when test="count(aspect[@name='last_complete_build']/object) > 0" >
                    <xsl:for-each select="aspect[@name='last_complete_build']/object">
                      <xsl:call-template name="genome_model_builds_list_table_row" />
                    </xsl:for-each>
                  </xsl:when>
                  <xsl:otherwise>
                    <tr>
                      <td></td>
                      <td colspan="6">
                        <strong>No builds found for this model.</strong>
                      </td>
                    </tr>
                  </xsl:otherwise>
                </xsl:choose>

              </tbody>
            </table>
          </xsl:if>
        </div> <!-- end box_content -->
      </div>
    </div><!-- end model_builds_list -->

  </xsl:template>


  <xsl:template name="genome_model_builds_list_table_row">
    <xsl:variable name="build_directory_url">
      <xsl:choose>
        <xsl:when test="aspect[@name='data_directory']/value != ''">
          <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="aspect[@name='data_directory']/value" />
        </xsl:when>
        <xsl:otherwise>

        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="b_status" select="aspect[@name='status']/value"/>
    <xsl:variable name="lc_b_status" select="translate($b_status,'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz')"/>

    <tr>
      <td>
        <xsl:call-template name="object_link_button">

          <xsl:with-param name="linktext" select="@id" />
          <xsl:with-param name="icon" select="'sm-icon-extlink'" />
        </xsl:call-template>

      </td>
      <td>
        <xsl:attribute name="class"><xsl:text>status </xsl:text><xsl:value-of select="$lc_b_status"/></xsl:attribute><xsl:value-of select="$lc_b_status"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='date_scheduled']/value"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='date_completed']/value"/>
      </td>
      <td class="buttons">

        <!-- data directory button -->
        <xsl:if test="$build_directory_url != ''">
          <a class="mini btn"><xsl:attribute name="href"><xsl:value-of select='$build_directory_url'/></xsl:attribute>data directory</a>
        </xsl:if>

        <!-- notes button -->
        <xsl:variable name="note_count" select="count(aspect[@name='notes']/object)"/>
        <xsl:variable name="build_id" select="@id"/>

        <xsl:choose>
          <xsl:when test="$note_count &gt; 0">

            <a class="mini btn notes-popup">
              <xsl:attribute name="title">Build <xsl:value-of select="$build_id"/> Notes</xsl:attribute>
              <xsl:attribute name="id"><xsl:value-of select="$build_id"/></xsl:attribute>
              <span class="sm-icon sm-icon-newwin"><br/></span>notes (<xsl:value-of select="$note_count"/>)
            </a>
            <!-- div for instrument data -->
            <div style="display: none;">
              <xsl:attribute name="id">notes_subject_<xsl:value-of select="$build_id"/></xsl:attribute>
              <table class="lister" border="0" width="100%" cellspacing="0" cellpadding="0">
                <colgroup>

                </colgroup>
                <thead>
                  <th>
                    header
                  </th>
                  <th>
                    date
                  </th>
                  <th>
                    editor id
                  </th>
                  <th>
                    <xsl:text> </xsl:text>
                  </th>

                </thead>
                <tbody>
                  <xsl:for-each select="aspect[@name='notes']/object">
                    <xsl:sort select="aspect[@name='entry_date_sort']/value" data-type="number" order="ascending"/>
                    <tr>
                      <td><strong><xsl:value-of select="aspect[@name='header_text']/value"/></strong></td>
                      <td><xsl:value-of select="aspect[@name='entry_date']/value"/></td>
                      <td><xsl:value-of select="aspect[@name='editor_id']/value"/></td>
                      <td class="buttons">
                        <xsl:for-each select="aspect[@name='subject']/object">
                          <xsl:call-template name="object_link_button">
                            <xsl:with-param name="linktext" select="'subject'"/>
                            <xsl:with-param name="icon" select="'sm-icon-extlink'"/>
                          </xsl:call-template>
                        </xsl:for-each>
                      </td>

                    </tr>
                    <tr>
                      <td colspan="4" class="text"><xsl:value-of select="aspect[@name='body_text']/value"/></td>
                    </tr>

                  </xsl:for-each>
                </tbody>
              </table>

            </div>
          </xsl:when>
        </xsl:choose>
      </td>
    </tr>

    <!-- show input deltas row if found -->
    <xsl:if test="count(aspect[@name='build_input_differences_from_model']/object) + count(aspect[@name='delta_model_input_differences_from_model']/object) > 0">
      <tr>
        <td></td>
        <td style="background-color: #ffeff1; text-align: center; "><span style="color: #a40014; font-weight: bold;">inputs differ:</span></td>
        <td colspan="3" class="inputs_differ">
          <table class="name-value inputs_differ" cellpadding="0" cellspacing="0" border="0">
            <tr>
              <td class="value">

                <xsl:call-template name="genome_model_input_compact_list">
                  <xsl:with-param name="input" select="./aspect[@name='build_input_differences_from_model']"/>
                  <xsl:with-param name="report_missing" select="./aspect[@name='delta_model_input_differences_from_model']"/>
                </xsl:call-template>
              </td>
            </tr>
          </table>
        </td>
      </tr>

    </xsl:if>
  </xsl:template>

  <!-- DEPRECATED TEMPLATES -->
  <!-- MODEL w/ BUILDS LIST TEMPLATES  -->

  <xsl:template name="genome_model_build_table">
    <xsl:param name="want_builds" select="1" />

    <xsl:comment>template: status/genome_model.xsl:genome_model_build_table</xsl:comment>

    <!-- Called on a node containing one or more object nodes of type model -->
    <table id="model_list" class="list" width="100%" cellspacing="0" cellpadding="0" border="0" style="clear:both">
      <colgroup>
        <col width="40%" />
        <col />
        <col />
        <col />
      </colgroup>
      <tbody>
        <xsl:for-each select="object[./types[./isa[@type='Genome::Model']]]">
          <xsl:sort select="aspect[@name='name']/value" data-type="text" order="ascending"/>
          <xsl:variable name="is_default" select="aspect[@name='is_default']/value" />
          <tr class="model_row_header">
            <td class="model_name">
              <xsl:if test="$is_default = 1">
                <!-- if this is the default model, show a nice little star -->
                <img class="default_report_star" src="/res/old/report_resources/apipe_dashboard/images/icons/star_16.png" width="16" height="16" absmiddle="middle" alt="Default Model"/>
              </xsl:if>
              <xsl:call-template name="object_link"/>
            </td>
            <td>
              <strong>model id: </strong><xsl:value-of select="@id"/>
            </td>
            <td><strong>run as: </strong><xsl:value-of select="aspect[@name='run_as']/value"/></td>
            <td><strong>created by: </strong><xsl:value-of select="aspect[@name='created_by']/value"/></td>
            <td class="last"><strong>scheduled: </strong><xsl:value-of select="aspect[@name='creation_date']/value"/></td>
          </tr>
          <tr class="model_row_subheader">
            <td colspan="4">
              <xsl:call-template name="genome_model_input_table"/>
            </td>
          </tr>
          <xsl:if test="$want_builds = 1">
            <tr>
              <td colspan="4" class="subtable_cell">
                <xsl:call-template name="genome_model_build_table_section"/>
              </td>
            </tr>
          </xsl:if>
        </xsl:for-each>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="genome_model_build_table_section">
    <xsl:comment>template: status/genome_model.xsl:genome_ model_build_table_section</xsl:comment>
    <table width="100%" cellpadding="0" cellspacing="0" border="0" class="subtable">
      <colgroup>
        <col width="25%" />
        <col width="15%"/>
        <col width="15%"/>
        <col width="15%"/>
        <col width="15%"/>
        <col width="15%"/>
      </colgroup>
      <thead>
        <th class="subtable_label">BUILDS</th>
        <th>build id</th>
        <th>status</th>
        <th>date scheduled</th>
        <th>date completed</th>
      </thead>
      <tbody>
        <xsl:choose>
          <xsl:when test="count(aspect[@name='builds']/object) > 0">
            <xsl:for-each select="aspect[@name='builds']/object">
              <xsl:call-template name="genome_model_build_table_row" />
            </xsl:for-each>
          </xsl:when>
          <xsl:when test="count(aspect[@name='last_succeeded_build']/object) > 0" >
            <xsl:for-each select="aspect[@name='last_succeeded_build']/object">
              <xsl:call-template name="genome_model_build_table_row" />
            </xsl:for-each>
          </xsl:when>
          <xsl:when test="count(aspect[@name='last_complete_build']/object) > 0" >
            <xsl:for-each select="aspect[@name='last_complete_build']/object">
              <xsl:call-template name="genome_model_build_table_row" />
            </xsl:for-each>
          </xsl:when>
          <xsl:otherwise>
            <tr>
              <td></td>
              <td colspan="5">
                <strong>No builds found for this model.</strong>
              </td>
            </tr>
          </xsl:otherwise>
        </xsl:choose>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="genome_model_build_table_row">
    <xsl:comment>template: status/genome_model.xsl:genome_ model_build_table_row</xsl:comment>
    <tr onmouseover="this.className = 'hover'" onmouseout="this.className=''">
      <xsl:attribute name="onclick">
        <xsl:text>javascript:document.location.href='</xsl:text>
        <xsl:call-template name="object_link_href" />
        <xsl:text>'</xsl:text>
      </xsl:attribute>
      <td>

      </td>
      <td>
        <xsl:call-template name="object_link">
          <xsl:with-param name="linktext"><xsl:value-of select="@id"/></xsl:with-param>
        </xsl:call-template>
      </td>
      <td><xsl:attribute name="class"><xsl:text>status </xsl:text><xsl:value-of select="aspect[@name='status']/value"/></xsl:attribute>
      <xsl:value-of select="aspect[@name='status']/value"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='date_scheduled']/value"/>
      </td>
      <td>
        <xsl:value-of select="aspect[@name='date_completed']/value"/>
      </td>
    </tr>
  </xsl:template>


  <xsl:template name="genome_model_flagstat_table">
    <xsl:comment>template: genome_model.xsl:genome_model_flagstats</xsl:comment>
    <!-- creates flagstat table in a jQueryUI popup. See G/M/V/Resource/Html/js/app/genome_model.js -->

    <div id="flagstat_table">
      <table class="name-value" width="100%">
        <colgroup>
          <col width="60%"/>
          <col width="40%"/>
        </colgroup>
        <tbody>
          <xsl:for-each select="item">
            <xsl:variable name="value" select="."/>
            <xsl:variable name="name" select="@key"/>
            <tr>
              <td class="name" style="white-space: normal;"><xsl:value-of select="translate($name, '_', ' ')"/>
              </td>
              <td class="value">
                <xsl:choose>
                  <xsl:when test="contains($name, 'percent')">
                    <xsl:value-of select="$value"/>%
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="format-number($value, '#,##0')"/>
                  </xsl:otherwise>
                </xsl:choose>
              </td>
            </tr>
          </xsl:for-each>
        </tbody>
      </table>
    </div>


  </xsl:template>

</xsl:stylesheet>
