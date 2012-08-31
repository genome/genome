<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_task" match="object[./types[./isa[@type='Genome::Task']]]">
    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Task '" />
      <xsl:with-param name="display_name" select="display_name" />
      <xsl:with-param name="icon" select="'genome_task_32'" />
    </xsl:call-template>

    <script src="/res/js/app/genome_task_updater.js"/>
    <span id="task-id" style="display:none"><xsl:value-of select="aspect[@name='id']/value"/></span>

    <div class="content rounded shadow">
      <div class="container">
        <div id="objects" class="span-24 last">

        <!-- details for this task-->
        <div class="span_16_box_masonry">
          <div class="box_header span-16 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-7 last">Task Attributes</h3></div>
            <div class="box_button">
            </div>
          </div>

          <div class="box_content rounded-bottom span-16 last">
            <table class="name-value">
              <tbody>
                <tr>
                  <td class="name">Task:
                  </td>
                  <td class="value"><xsl:value-of select="display_name"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Task Status:
                  </td>
                  <td class="value" id="task-status"><xsl:value-of select="aspect[@name='status']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Submitted By:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='user_name']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Submitted On:
                  </td>
                  <td class="value"><xsl:value-of select="aspect[@name='time_submitted']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Started On:
                  </td>
                  <td class="value" id="time-started"><xsl:value-of select="aspect[@name='time_started']/value"/>
                  </td>
                </tr>
                <tr>
                  <td class="name">Completed On:
                  </td>
                  <td class="value" id="time-finished"><xsl:value-of select="aspect[@name='time_finished']/value"/>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        <div class="span_24_box_masonry">
            <div class="box_header span-24 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-7 last">Output Messages</h3></div>
            <div class="box_button">
            </div>
        </div>

        <div class="box_content rounded-bottom span-24 last">
        <pre id="stdout-content">
            Output Messages Will Appear Here
        </pre>
        <span class="spinner" style="display:none">Updating... <img src="/res/img/spinner.gif"/></span>
        </div>
        </div>

        <div class="span_24_box_masonry">
            <div class="box_header span-24 last rounded-top">
            <div class="box_title"><h3 class="nontyped span-7 last">Error Messages</h3></div>
            <div class="box_button">
            </div>
        </div>

        <div class="box_content rounded-bottom span-24 last">
        <pre id="stderr-content">
            Error Messages, If Any, Will Appear Here
        </pre>
        <span class="spinner" style="display:none">Updating... <img src="/res/img/spinner.gif"/></span>
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


</xsl:stylesheet>
