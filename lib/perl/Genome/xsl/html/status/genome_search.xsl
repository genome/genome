<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">


  <xsl:template name="genome_search" match="object[./types[./isa[@type='Genome::Search']]]">

    <xsl:comment> ***************************************************** </xsl:comment>
    <xsl:comment> Interpreter: <xsl:value-of select="/object/aspect[@name='interpreter']/value"/> </xsl:comment>
    <xsl:comment> Genome: <xsl:value-of select="/object/aspect[@name='genome_path']/value"/> </xsl:comment>
    <xsl:comment> UR: <xsl:value-of select="/object/aspect[@name='ur_path']/value"/> </xsl:comment>
    <xsl:comment> ***************************************************** </xsl:comment>

    <script type="text/javascript" src="/res/js/app/genome_search.js"></script>

    <xsl:call-template name="control_bar_app"/>

    <div class="content rounded shadow" style="padding-top: 0;">
      <xsl:call-template name="app_header">
        <xsl:with-param name="app_name" select="'Analysis Search'"/>
        <xsl:with-param name="icon" select="'app_analysis_search_32'"/>
      </xsl:call-template>

      <div class="container">
        <div class="span-12">
          <div class="main_search">
            <form id="searchForm" method="get" action="/view/genome/search/query/status.html">
            <br/>

              <table cellpadding="0" cellspacing="0" border="0" class="search_elements">
                <tr>
                  <td>
                    <input class="query_box rounded" type="text" id="searchBox" name="query"/>
                  </td>
                  <td>
                    <input id="searchButton" type="submit" class="button" value="Search"/>
                  </td>
                </tr>
              </table>



            </form>
          </div>
<br/>
        </div> <!-- end .span-12 -->
        <div class="span-12 last">
          <br/>
        </div>
        <hr class="space"/>

        <div class="main_search_hints clearfix">
          <div class="box_header span-8 last rounded-top">
            <div class="box_title"><h3 class="nontyped last">Analysis Tools</h3></div>
          </div>

          <div class="box_content rounded-bottom span-24 last">
            <div style="width: 50%; float: left;">
              <div class="padding10">
                <a id="project_link">Your projects</a> - Use projects to organize your analysis related data. 
                (<a href="/view/genome/project/set/help.html">read more</a>)<br/>

                <a href="/view/genome/nomenclature/set/status.html">Nomenclatures</a> - Nomenclatures describe different formats for clinical data<br/>

                <a href="/view/genome/subject/set/create.html">Import subject data</a> - Upload key/value pair clinical data about samples<br/>

<br/>
<br/>
Email <a href="mailto:apipe@genome.wustl.edu">apipe@genome.wustl.edu</a> if you have any questions or suggestions.
              </div><!-- end .padding10 -->
            </div>
          </div> <!-- end .box_content -->
        </div><!-- end .main_search_hints -->

      </div> <!-- end .container  -->
    </div> <!-- end .content  -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>
