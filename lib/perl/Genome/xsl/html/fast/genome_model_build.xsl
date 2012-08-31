<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- full page display for a project -->
  <xsl:template name="genome_model_build" match="/object[./types[./isa[@type='Genome::Model::Build']]]">
    <xsl:comment>template: fast/genome_model_build.xsl match: object[./types[./isa[@type='Genome::Model::Build']]]</xsl:comment>

    <script type="text/javascript">
          window.page_data = function(){
            return {
                workflow: {
                "id": "<xsl:value-of select="/object/aspect[@name='_newest_workflow_instance']/object/@id"/>"
                },
                stages: {
                    "count": <xsl:value-of select="count(build/stages/stage)"/>
                }
            }
          }();
    </script>


    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Build'" />
      <xsl:with-param name="display_name" select="/object/display_name" />
      <xsl:with-param name="icon" select="'genome_model_build_32'" />
    </xsl:call-template>

    <script type="text/javascript" src="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.js"></script>
    <script type="text/javascript" src="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/js/jquery.dataTables.plugin.formatted-num.js"></script>
    <link rel="stylesheet" href="https://imp.gsc.wustl.edu/resources/report_resources/jquery/dataTables-1.5/media/css/gc_table.css" type="text/css" media="screen"></link>

    <link rel="stylesheet" href="/res/css/genome_model_build.css" type="text/css" />

    <script type='text/javascript' src='/res/js/app/fast/genome_model_build.js'></script>


    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

            <xsl:comment>left side attribute box</xsl:comment>
            <div class="span_8_box_masonry attrib-box">
                    <div class="box_header rounded-top span-8 last">
                        <div class="box_title"><h3 class="nontyped span-7 last">Attributes</h3></div>
                    </div>
                    <div class="box_content rounded-bottom span-8 last"><br/>
                        <dl>
                            <dt>status</dt>
                            <dd>
                                <xsl:value-of select="aspect[@name='the_master_event']/object/aspect[@name='event_status']/value"/><xsl:text>  </xsl:text>
                                <xsl:call-template name="notes_button"></xsl:call-template>
                            </dd>

                            <dt>run by</dt>
                            <dd><xsl:value-of select="aspect[@name='run_by']/value"/></dd>

                            <dt>model</dt>
                            <dd><xsl:value-of select="aspect[@name='model']/object/aspect[@name='name']"/>
                                &#160;<a>
                                <xsl:attribute name="href">
                                    <xsl:for-each select="/object/aspect[@name='model']/object">
                                        <xsl:call-template name="object_link_href"></xsl:call-template>
                                    </xsl:for-each> 
                                </xsl:attribute>
                                <xsl:value-of select="aspect[@name='model']/object/@id"/>
                                </a>
                            </dd>

                            <dt>subject</dt>
                            <dd><xsl:value-of select="aspect[@name='model']/object/aspect[@name='subject']/object/display_name"/>
                               &#160;<a>
                                <xsl:attribute name="href">
                                    <xsl:for-each select="aspect[@name='model']/object/aspect[@name='subject']/object">
                                        <xsl:call-template name="object_link_href"></xsl:call-template>
                                    </xsl:for-each> 
                                </xsl:attribute>
                                <xsl:value-of select="aspect[@name='model']/object/aspect[@name='subject']/object/@id"/>
                                </a>
                            </dd>

                            <xsl:variable name="data_directory" select="aspect[@name='data_directory']/value"/>
                            <dt>files</dt>
                            <dd>
                                <a class="mini btn">
                                    <xsl:attribute name="href">
                                        <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="$data_directory"/>
                                    </xsl:attribute>
                                    <xsl:text>data dir</xsl:text>
                                </a>
                                <a class="mini btn">
                                    <xsl:attribute name="href">
                                        <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="$data_directory"/>
                                        <xsl:text>/logs/</xsl:text>
                                        <xsl:value-of select="/object/aspect[@name='the_master_event']/object/@id"/>
                                        <xsl:text>.err</xsl:text>
                                    </xsl:attribute>
                                    <xsl:text>error log</xsl:text>
                                </a>
                                <a class="mini btn">
                                    <xsl:attribute name="href">
                                        <xsl:text>https://gscweb.gsc.wustl.edu/</xsl:text><xsl:value-of select="$data_directory"/>
                                        <xsl:text>/logs/</xsl:text>
                                        <xsl:value-of select="/object/aspect[@name='the_master_event']/object/@id"/>
                                        <xsl:text>.out</xsl:text>
                                    </xsl:attribute>
                                    <xsl:text>output log</xsl:text>
                                </a>
                            </dd>
                        </dl>
                    </div>
            </div>

            <xsl:comment>right side attribute box</xsl:comment>
            <div class="span_8_box_masonry attrib-box">
                    <div class="box_header rounded-top span-8 last">
                        <div class="box_title"><h3 class="nontyped span-7 last"></h3></div>
                    </div>
                    <div class="box_content rounded-bottom span-8 last"> 
                        <dl>
                            <dt>software</dt>
                            <dd><xsl:value-of select="aspect[@name='software_revision']/value"/></dd>

                            <dt>processing profile</dt>
                            <dd><xsl:value-of select="aspect[@name='model']/object/aspect[@name='processing_profile']/object/aspect[@name='name']/value"/>
                               &#160;<a>
                                <xsl:attribute name="href">
                                    <xsl:for-each select="aspect[@name='model']/object/aspect[@name='processing_profile']/object">
                                        <xsl:call-template name="object_link_href"></xsl:call-template>
                                    </xsl:for-each> 
                                </xsl:attribute>
                                <xsl:value-of select="aspect[@name='model']/object/aspect[@name='processing_profile']/object/@id"/>
                                </a>
                            </dd>

                            <dt>workflow</dt>
                            <dd>
                                <a>
                                <xsl:attribute name="href">
                                    <xsl:for-each select="/object/aspect[@name='_newest_workflow_instance']/object">
                                        <xsl:call-template name="object_link_href"></xsl:call-template>
                                    </xsl:for-each> 
                                </xsl:attribute>
                                <xsl:value-of select="/object/aspect[@name='_newest_workflow_instance']/object/@id"/>
                                </a>
                            </dd>

                            <dt>LSF job</dt>
                            <dd><xsl:value-of select="aspect[@name='the_master_event']/object/aspect[@name='lsf_job_id']/value"/></dd>
                        </dl>
                    </div>
            </div>

            <div style="clear:both">
            </div>

        </div> <!-- end .objects -->

            <xsl:comment>inputs box</xsl:comment>
            <xsl:if test="count(/object/aspect[@name='inputs']/object) > 0 ">
                <xsl:for-each select="/object/aspect[@name='inputs']">
                    <xsl:call-template name="genome_model_input_table"/>
                </xsl:for-each>
            </xsl:if>

        <xsl:call-template name="genome_model_build_tabs"/>

      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>


  <xsl:template name="notes_button">

        <xsl:variable name="note_count" select="count(/object/aspect[@name='notes']/object)"/>
        <xsl:variable name="build_id" select="build/@build-id"/>

        <xsl:choose>
          <xsl:when test="$note_count &gt; 0">
            <a class="notes-link mini btn notes-popup">
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
                  <xsl:for-each select="/object/aspect[@name='notes']/object">
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

  <xsl:template name="genome_model_build_tabs">

    <xsl:variable name="workflow" select="/object/aspect[@name='_newest_workflow_instance']/object"/>

        <div id="process_tabs" class="span-24 last">
          <ul>
            <li class="tab_header"><h3 class="genome_processingprofile_16">Build Process</h3></li>
            <li><a><xsl:attribute name="href">/viewajax/workflow/operation/instance/statuspopup.html?id=<xsl:value-of select="$workflow/@id"/></xsl:attribute>
                    <span class="spinner"><xsl:text> </xsl:text></span>workflow
            </a></li>

            <li><a href="#events">events</a></li>
          </ul>

          <div id="events">
            <div id="eventview">
              <xsl:choose>
                <xsl:when test="count(build/stages/stage) > 0">
                  <table border="0" cellpadding="0" cellspacing="0" class="stages" width="100%">
                    <tr>
                      <xsl:for-each select="build/stages/stage[count(command_classes/*) > 0]">
                        <td>
                          <table class="stage" border="0" cellpadding="0" cellspacing="0" width="100%">

                            <tr>
                              <th colspan="2">
                                <xsl:variable name="stage_name" select="@value"/>
                                <xsl:value-of select="translate($stage_name,'_', ' ')"/>
                              </th>
                            </tr>

                            <xsl:variable name="num_succeeded" select="count(descendant::*/event_status[text()='Succeeded'])"/>
                            <xsl:variable name="num_succeeded_label">
                              <xsl:choose>
                                <xsl:when test="$num_succeeded = 0">ghost</xsl:when>
                                <xsl:otherwise>label</xsl:otherwise>
                              </xsl:choose>
                            </xsl:variable>

                            <xsl:variable name="num_scheduled" select="count(descendant::*/event_status[text()='Scheduled'])"/>
                            <xsl:variable name="num_scheduled_label">
                              <xsl:choose>
                                <xsl:when test="$num_scheduled = 0">ghost</xsl:when>
                                <xsl:otherwise>label</xsl:otherwise>
                              </xsl:choose>
                            </xsl:variable>

                            <xsl:variable name="num_running" select="count(descendant::*/event_status[text()='Running'])"/>
                            <xsl:variable name="num_running_label">
                              <xsl:choose>
                                <xsl:when test="$num_running = 0">ghost</xsl:when>
                                <xsl:otherwise>label</xsl:otherwise>
                              </xsl:choose>
                            </xsl:variable>

                            <xsl:variable name="num_abandoned" select="count(descendant::*/event_status[text()='Abandoned'])"/>
                            <xsl:variable name="num_abandoned_label">
                              <xsl:choose>
                                <xsl:when test="$num_abandoned = 0">ghost</xsl:when>
                                <xsl:otherwise>label</xsl:otherwise>
                              </xsl:choose>
                            </xsl:variable>

                            <xsl:variable name="num_failed" select="count(descendant::*/event_status[text()='Crashed' or text()='Failed'])"/>
                            <xsl:variable name="num_failed_label">
                              <xsl:choose>
                                <xsl:when test="$num_failed = 0">ghost</xsl:when>
                                <xsl:otherwise>label</xsl:otherwise>
                              </xsl:choose>
                            </xsl:variable>


                            <tr><xsl:attribute name="class"><xsl:value-of select="$num_scheduled_label"/></xsl:attribute>
                            <td class="label">
                              Scheduled:
                            </td>
                            <td class="value">
                              <xsl:value-of select="$num_scheduled"/>
                            </td>

                            </tr>

                            <tr><xsl:attribute name="class"><xsl:value-of select="$num_running_label"/></xsl:attribute>
                            <td class="label">
                              Running:
                            </td>
                            <td class="value">
                              <xsl:value-of select="$num_running"/>
                            </td>
                            </tr>

                            <tr><xsl:attribute name="class"><xsl:value-of select="$num_succeeded_label"/></xsl:attribute>
                            <td class="label">
                              Succeeded:
                            </td>
                            <td class="value">
                              <xsl:value-of select="$num_succeeded"/>
                            </td>
                            </tr>

                            <tr><xsl:attribute name="class"><xsl:value-of select="$num_abandoned_label"/></xsl:attribute>
                            <td class="label">
                              Abandoned:
                            </td>
                            <td class="value">
                              <xsl:value-of select="$num_abandoned"/>
                            </td>
                            </tr>

                            <tr><xsl:attribute name="class"><xsl:value-of select="$num_failed_label"/></xsl:attribute>
                            <td class="label">
                              Crashed/Failed:
                            </td>
                            <td class="value">
                              <xsl:value-of select="$num_failed"/>
                            </td>
                            </tr>

                            <tr class="total">
                              <td class="label">
                                Total:
                              </td>
                              <td class="value">
                                <xsl:value-of select="count(descendant::*/event_status)"/>
                              </td>
                            </tr>
                          </table>
                        </td>
                      </xsl:for-each>
                    </tr>
                  </table>

                  <hr/>

                  <xsl:for-each select="//stage[count(command_classes/*) > 0 ]">
                    <h2>
                      <xsl:variable name="stage_name" select="@value"/>
                      <xsl:value-of select="translate($stage_name,'_', ' ')"/>
                    </h2>
                    <table class="lister" width="100%" cellspacing="0" cellpadding="0" border="0">
                      <xsl:attribute name="id"><xsl:value-of select="@value"/></xsl:attribute>
                      <colgroup>
                        <col width="40%"/>
                        <col/>
                        <col/>
                        <col/>
                        <col/>
                        <col/>
                      </colgroup>
                      <thead>
                        <tr>
                          <th>
                            <xsl:choose><xsl:when test="@value='alignment'">flow cell</xsl:when>
                            <xsl:otherwise>event</xsl:otherwise>
                            </xsl:choose>
                          </th>

                          <th>status</th>
                          <th>scheduled</th>
                          <th>completed</th>
                          <th class="last">elapsed</th>
                          <th><br/></th>
                        </tr>
                      </thead>
                      <tbody>
                        <xsl:for-each select="descendant::*/event">
                          <tr>
                            <xsl:variable name="command_class" select="@command_class"/>
                            <td>
                              <xsl:variable name="containing_command_class" select="../../@value"/>
                              <xsl:choose>
                                <!-- if command_class contains 'AlignReads' there should be instrument data associated -->
                                <xsl:when test="contains($command_class, 'AlignReads') or contains($containing_command_class, 'AlignReads')">
                                  <xsl:variable name="inst_data_id" select="instrument_data_id" />
                                  <xsl:variable name="inst_data_count" select="count(//instrument_data[@id=$inst_data_id])"/>
                                  <xsl:choose>
                                    <!-- if we have instrument data element(s), show flow cell and lane -->
                                    <xsl:when test="$inst_data_count > 0">
                                      <xsl:for-each select="//instrument_data[@id=$inst_data_id]" >
                                        <xsl:value-of select="flow_cell_id"/><xsl:text disable-output-escaping="yes"> Lane: </xsl:text><xsl:value-of select="lane"/>
                                      </xsl:for-each>
                                      <xsl:if test="filter_desc"><xsl:text disable-output-escaping="yes"> </xsl:text>(<xsl:value-of select="filter_desc"/>)</xsl:if>
                                    </xsl:when>
                                    <xsl:otherwise>
                                      <!-- no instrument data elements, show a warning message -->
                                      <span style="color: #933; font-weight: bold;">No instrument data found for this lane.</span>
                                    </xsl:otherwise>
                                  </xsl:choose>
                                </xsl:when>
                                <xsl:otherwise>
                                  <!-- event is not expected to have instrument data, so show command class -->
                                  <xsl:variable name="full_command_class" select="@command_class" />
                                  <xsl:value-of select="substring-after($full_command_class,'Genome::Model::Build::Command::')"/>
                                  <xsl:value-of select="substring-after($full_command_class,'Genome::Model::Event::Build::')"/>
                                </xsl:otherwise>
                              </xsl:choose>
                            </td>
                            <td>
                              <xsl:variable name="e_status" select="event_status"/>
                              <xsl:variable name="lc_e_status" select="translate($e_status,'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz')"/>
                              <xsl:attribute name="class"><xsl:text disable-output-escaping="yes">status </xsl:text><xsl:value-of select="$lc_e_status"/></xsl:attribute>
                              <xsl:value-of select="$lc_e_status"/>
                            </td>
                            <!-- <td><xsl:value-of select="event_status"/></td> -->
                            <td><xsl:value-of select="date_scheduled"/></td>
                            <td><xsl:value-of select="date_completed"/></td>
                            <td class="last"><xsl:value-of select="elapsed_time"/></td>
                            <td class="buttons">
                              <xsl:comment>
                                rest: <xsl:value-of select="$rest"/>
                              </xsl:comment>

                              <a class="mini btn event-popup-ajax">
                                <!-- <xsl:attribute name="onclick">
                                     <xsl:text>javascript:status_popup('</xsl:text><xsl:value-of select="@id"/><xsl:text>','</xsl:text><xsl:value-of select="$currentLink"/><xsl:text>','</xsl:text><xsl:value-of select="@id"/><xsl:text>');</xsl:text>
                                     </xsl:attribute> -->
                                <xsl:attribute name="href">
                                  <xsl:call-template name="object_link_href">
                                    <xsl:with-param name="type" select="$command_class"/>
                                    <xsl:with-param name="rest" select="'/viewajax'" />
                                    <xsl:with-param name="key" select="'id'"/>
                                    <xsl:with-param name="id" select="@id"/>
                                    <xsl:with-param name="perspective" select="'statuspopup'"/>
                                    <xsl:with-param name="toolkit" select="'html'"/>
                                  </xsl:call-template>
                                </xsl:attribute>
                                <xsl:attribute name="title">
                                  <xsl:value-of select="@id"/>
                                </xsl:attribute>
                              <span class="sm-icon sm-icon-newwin"><br/></span>event details</a>
                            </td>
                          </tr>
                        </xsl:for-each>
                      </tbody>
                    </table>
                  </xsl:for-each>

                </xsl:when>
                <xsl:otherwise>
                  <p style="padding-top: 25px; text-align: center"><strong>No events found.</strong></p>
                </xsl:otherwise>
              </xsl:choose>
            </div> <!-- end .eventview -->
          </div> <!-- end .events -->

        </div>

</xsl:template>




</xsl:stylesheet>
