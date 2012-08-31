<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_sys_service_set" match="object[./types[./isa[@type='Genome::Sys::Service::Set']]]">

    <script type='text/javascript' src='/res/js/pkg/boxy/javascripts/jquery.boxy.js'></script>
    <link rel="stylesheet" href="/res/js/pkg/boxy/stylesheets/boxy.css" type="text/css" />
    <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Services'" />
    </xsl:call-template>

    <div class="content rounded shadow">
      <div class="container">
          <xsl:for-each select="aspect[@name='members']/object">
              <xsl:call-template name="genome_service_box"/>
          </xsl:for-each>
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="genome_service_box">

    <div class="span_12_box_masonry">
      <div class="box_header span-12 last rounded-top">
          <div class="box_title"><h3 class="genome_service_16 span-7 last">
              <xsl:value-of select="normalize-space(aspect[@name='name']/value)"/>
          </h3></div>
      </div>

      <div class="box_content rounded-bottom span-12 last">
        <table class="name-value">
          <tbody>

            <tr>
              <td class="name">
                  <button class="restart">Restart</button>
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='restart_command']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">
                  <button class="stop">Stop</button>
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='stop_command']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">Host:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='host']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">Log Path:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='log_path']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">URL:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='url']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">Status:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='status']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">PID Name:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='pid_name']/value)"/>
              </td>
            </tr>

            <tr>
              <td class="name">PID Status:
              </td>
              <td class="value">
                  <xsl:value-of select="normalize-space(aspect[@name='pid_status']/value)"/>
              </td>
            </tr>

          </tbody>
        </table>
      </div>
    </div>
  </xsl:template>
</xsl:stylesheet>
