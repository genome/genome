<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_search_query" match="/solr-results">
    <xsl:call-template name="control_bar_app"/>

    <div class="content rounded shadow" style="padding-top: 0;">
      <xsl:call-template name="app_header">
        <xsl:with-param name="app_name" select="'Analysis Search'"/>
        <xsl:with-param name="icon" select="'app_analysis_search_32'"/>
      </xsl:call-template>

<!--    <script type="text/javascript" src="/res/js/pkg/jquery.js"></script> -->
<!-- <script src="/res/js/pkg/jquery-1.5.1.min.js" type="text/javascript"></script>
<script src="/res/js/pkg/jquery-ui-1.8.11.custom.min.js" type="text/javascript"></script>
-->

<script type="text/javascript" src="/res/js/app/genome_search_query.js"></script>
<script type="text/javascript" src="/res/js/app/genome_projectbox.js"></script>


<style>
    .projectContainer {
        text-align: left;
    }

    .clickable {
        cursor: pointer;
    }

    .ui-accordion-content {
        background: #CCCCCC;
        width: 400px;
        height: 200px;
    }

    .ui-accordion-header {
        background: #FAAFBA;
        width: 400px;
    }
</style>

        <!-- background-color: #CBE0FF; -->


      <div class="container">
        <div class="span-5">
          <br/>
        </div>
        <div class="span-19 last">
          <div id="search_box" class="main_search" style="margin-left: 15px; margin-bottom: 15px;">
            <form action="">
              <table cellpadding="0" cellspacing="0" border="0" class="search_elements">
                <tr>
                  <td>
                    <input class="query_box rounded" type="text" name="query"><xsl:attribute name="value"><xsl:value-of select="//@query"/></xsl:attribute></input>
                  </td>
                  <td>
                    <input type="submit" class="button" value="Search"/>
                  </td>
                </tr>
              </table>
            </form>
            <p class="small">Didn't find what you're looking for? Please <a href="mailto:apipe@genome.wustl.edu">let us know</a>.</p>
            <div id="errorAccordion" class="errorAccordion">
                <h3><a href="#">Error Messages</a></h3>
                <div id="errorMsg">Something went wrong.</div>
            </div>

          </div>
        </div>

        <div class="span-5">
          <xsl:choose>
            <xsl:when test="/solr-results/@num-found &gt; 0">
              <div class="sidebar_search rounded-right">
                <h4>Show:</h4>
                  <xsl:variable name="selected_facet" select="//@facet-name"/>
                <ul>

                  <li>
                    <xsl:if test="not(//@facet-name)">
                      <xsl:attribute name="class">current</xsl:attribute>
                    </xsl:if>
                    <div class="category">
                      <a>
                        <xsl:attribute name="href">/view/genome/search/query/status.html?query=<xsl:value-of select="/solr-results/@query-no-types"/></xsl:attribute>
                        All
                      </a>
                        (<xsl:value-of select="//@facet-total"/>)
                    </div>
                  </li>

                  <xsl:for-each select="facets/field">

                    <xsl:variable name="current_facet" select="@name"/>
                    <li>
                      <xsl:if test="($selected_facet = $current_facet)">
                        <xsl:attribute name="class">current</xsl:attribute>
                      </xsl:if>

                      <div><xsl:attribute name="class">icon16 <xsl:value-of select="@icon-prefix"/>_16</xsl:attribute><br/></div>
                      <div class="category">
                        <a>
                          <xsl:attribute name="href">/view/genome/search/query/status.html?query=<xsl:value-of select="/solr-results/@query-no-types"/>&amp;fq=type:"<xsl:value-of select="@name"/>"
                          </xsl:attribute>
                          <xsl:value-of select="@label"/>
                        </a>
                      (<xsl:value-of select="@count"/>)</div>
                    </li>

                  </xsl:for-each>
                </ul>
              </div>
            </xsl:when>
            <xsl:otherwise>
              <br/>
            </xsl:otherwise>
          </xsl:choose>

<br/>

            <div id="myProjectBox" class="sidebar_project rounded-right">
                <h4 style="float: left" id="myProjectsCount"></h4>
                <br/>

                <ul id="projectBox"> </ul>

                <div style="margin: 2px 0px 0px 10px; color: red" id="loadingStatus" />

                <a class="clickable" id="createProjectLink">Create a new project</a> <br/>
                <div id="createProjectDiv">
                    <form id="createProject">
                    New project name:<br/>
                    <input id="createProjectName" name="name" type="text"/> <input id="submitCreateProject" type="submit" value="create"/>
                    </form>
                </div>
            </div>

        </div>


        <div class="span-19 last">
          <div style="margin: 0 15px;">

            <!--            <xsl:apply-templates/> -->
            <xsl:for-each select="/solr-results/doc">
                <xsl:call-template name="universal_search_result"/>
            </xsl:for-each>


            <div class="pager">
              <div class="nav">

                <xsl:choose>
                  <xsl:when test="string(page-info/@previous-page)">
                    <a class="mini btn">
                      <xsl:attribute name="href">
                            <xsl:text>/view/genome/search/query/status.html?query=</xsl:text>
                            <xsl:value-of select="@query"/><xsl:if test="//@facet-name">&amp;fq=type:<xsl:value-of select="//@facet-name"/></xsl:if>
                            <xsl:text>&amp;page=</xsl:text><xsl:value-of select="page-info/@previous-page" />
                      </xsl:attribute>
                      <span class="sm-icon sm-icon-triangle-1-w"><br/></span>
                    </a>
                  </xsl:when>
                  <xsl:otherwise>
                    <a href="#" class="grey mini btn"><span class="sm-icon sm-icon-triangle-1-w"><br/></span></a>
                  </xsl:otherwise>
                </xsl:choose>
              </div>

              <div class="position">Page <xsl:value-of select="page-info/@current-page" /> of <xsl:value-of select="page-info/@last-page" /></div>

              <div class="nav">
                <xsl:choose>
                  <xsl:when test="string(page-info/@next-page)">
                    <a class="mini btn">
                        <xsl:attribute name="href">
                            <xsl:text>/view/genome/search/query/status.html?query=</xsl:text>
                            <xsl:value-of select="@query"/><xsl:if test="//@facet-name">&amp;fq=type:"<xsl:value-of select="//@facet-name"/>"</xsl:if>
                            <xsl:text>&amp;page=</xsl:text><xsl:value-of select="page-info/@next-page" />
                        </xsl:attribute>
                        <span class="sm-icon sm-icon-triangle-1-e"><br/></span>
                    </a>
                  </xsl:when>
                  <xsl:otherwise>
                    <a href="#" class="grey mini btn"><span class="sm-icon sm-icon-triangle-1-e"><br/></span></a>
                  </xsl:otherwise>
                </xsl:choose>
              </div>
            </div> <!-- end pager -->

          </div>
        </div> <!-- end span-20 -->
      </div> <!-- end container  -->
    </div> <!-- end content  -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

<!--
<xsl:template match="doc">
<p>hi</p>
</xsl:template>
-->

</xsl:stylesheet>
