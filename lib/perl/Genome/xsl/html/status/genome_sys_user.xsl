<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_sys_user" match="object[./types[./isa[@type='Genome::Sys::User']]]">
    <xsl:comment>template: /html/status/genome_sys_user.xsl  match: object[./types[./isa[@type='Genome::Sys::User']]]</xsl:comment>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'User:'" />
      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
      <xsl:with-param name="icon" select="'genome_sys_user_32'" />
    </xsl:call-template>

    <script type="text/javascript" src="/res/js/app/genome_sys_user.js"></script>

    <link rel="stylesheet" href="/res/js/pkg/TableTools/media/css/TableTools.css" media="screen"/>
    <script type="text/javascript" src="/res/js/pkg/dataTables/media/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" src="/res/js/pkg/TableTools/media/js/TableTools.min.js"></script>

    <script type="text/javascript" charset="utf-8">
        window.username = "<xsl:value-of select="aspect[@name='username']/value"/>";
        window.email = "<xsl:value-of select="aspect[@name='email']/value"/>";
    </script>

    <div class="content rounded shadow">
      <div class="project_container container">

          <div id="user_data" class="box_content rounded-bottom last" style="width: 35%">
              <div class="padding10">

                <img alt="No photo available" style="width: 160px; height: 120px;">
                    <xsl:attribute name="src">http://gscweb.gsc.wustl.edu/p/<xsl:value-of select="aspect[@name='username']/value"/>.jpg
                    </xsl:attribute>
                </img>
                <br/><br/>

                Send mail to 
                <a>
                    <xsl:attribute name="href">
                    http://mailto:<xsl:value-of select="/object/aspect[@name='email']/value"/>
                    </xsl:attribute>
                <xsl:value-of select="/object/aspect[@name='email']/value"/>
                </a>
                <br/>

              </div><!-- end .padding10 -->
          </div> <!-- end user_data .box_content -->

       
        <br/> 
        <div> <h3>Projects:</h3> </div>

        <table id="projects" class="lister datatable">
            <thead><tr id="table-header"><th></th><th>Name</th><th># Items</th></tr></thead>
            <tbody>
                <xsl:for-each select="/object/aspect[@name='projects']/object">
                    <tr>
                        <td><input class="selectionCheckbox" type="checkbox"/></td>
                        <td>
                            <a>
                                <xsl:attribute name="href">
                                    /view/genome/project/status.html?id=<xsl:value-of select="@id"/>
                                </xsl:attribute>
                            <xsl:value-of select="aspect[@name='name']/value"/>
                            </a>
                        </td>
                        <td>
                            <xsl:value-of select="aspect[@name='parts_count']/value"/>
                        </td>
                    </tr>
                </xsl:for-each>
            </tbody>
        </table>
<br/>
<a href="/view/genome/project/set/help.html">How do I use projects?</a>

      </div> <!-- end container -->
    </div> <!-- end content -->

  </xsl:template>

</xsl:stylesheet>



