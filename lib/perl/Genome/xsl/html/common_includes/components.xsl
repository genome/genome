<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:rest="urn:rest"
                xmlns:str="http://xsltsl.org/string">

  <xsl:strip-space elements="*"/>

  <xsl:template name="html_head_page">
    <xsl:comment>template: /html/common_includes/components.xsl name="html_head_page"</xsl:comment>

    <xsl:variable name="title">
      <xsl:choose>
        <!-- check to see if we're provided default UR generated title vars -->
        <xsl:when test="/object/display_name">
          <xsl:value-of select="/object/display_name" />
          <!-- <xsl:value-of select="/object/label_name" /><xsl:text> </xsl:text> -->
          <xsl:text> &#8211; </xsl:text>
          <xsl:call-template name="str:capitalise">
            <xsl:with-param name="text">
              <xsl:value-of select="$currentPerspective"/>
            </xsl:with-param>
          </xsl:call-template>
        </xsl:when>

        <!-- otherwise use XSL vars (see UR::Object::View::Default::Xsl) -->
        <xsl:otherwise>
          <xsl:value-of select="$displayName" />
          <!-- add object ID in parans if it's not already in the display name -->
          <xsl:if test="$objectId">
            <xsl:variable name="displayNameContainsId">
              <xsl:call-template name="str:string-match">
                <xsl:with-param name="text" select="$objectId"/>
                <xsl:with-param name="pattern" select="$displayName"/>
              </xsl:call-template>
            </xsl:variable>
            <xsl:if test="not($displayNameContainsId)">
              <xsl:text> (</xsl:text>
              <xsl:value-of select="$objectId"/>
              <xsl:text>)</xsl:text>
            </xsl:if>
          </xsl:if>
          <!-- <xsl:value-of select="$labelName" /> -->
          <xsl:text> &#8211; </xsl:text>
          <xsl:call-template name="str:capitalise">
            <xsl:with-param name="text" select="$currentPerspective"/>
          </xsl:call-template>
          <xsl:text>: </xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <title><xsl:value-of select="normalize-space($title)"/></title>
    <link rel="shortcut icon" href="/res/img/gc_favicon.png" type="image/png" />

    <!-- blueprint CSS -->
    <link rel="stylesheet" href="/res/css/blueprint/screen.css" type="text/css" media="screen, projection"/>
    <link rel="stylesheet" href="/res/css/blueprint/print.css" type="text/css" media="print"/>
    <!--[if lt IE 8]><link rel="stylesheet" href="/css/blueprint/ie.css" type="text/css" media="screen, projection"><![endif]-->

    <!-- report styles  -->
    <link rel="stylesheet" href="/res/css/master.css" type="text/css" media="screen, projection"/>
    <link rel="stylesheet" href="/res/css/buttons.css" type="text/css" media="screen, projection"/>
    <link rel="stylesheet" href="/res/css/icons.css" type="text/css" media="screen, projection"/>
    <link rel="stylesheet" href="/res/css/forms.css" type="text/css" media="screen, projection"/>

    <!-- GENOME_DEV_MODE is removed, if this is to be renabled, check for window.console instead -->
    <!--    <xsl:if test="$GENOME_DEV_MODE != 1"> -->
    <!-- prevent console() statements from being called on non-dev servers -->
    <!--      <script type="text/javascript" src="/res/js/app/console.js"></script> -->
    <!--    </xsl:if> -->

    <!-- jquery and jquery UI -->
    <link type="text/css" href="/res/js/pkg/jquery-ui-1.8.1.custom/css/gsc-theme/jquery-ui-1.8.1.custom.css" rel="stylesheet" />
    <link href="/res/css/jquery-ui-overrides.css" type="text/css" rel="stylesheet" media="screen, projection"/>
    <script type="text/javascript" src="/res/js/pkg/jquery.js"></script>
    <script type="text/javascript" src="/res/js/pkg/jquery-ui.js"></script>

    <!-- jquery.easydate for natural time displays -->
    <script type="text/javascript" src="/res/js/pkg/jquery.easydate.js"></script>

    <!-- jquery.masonry to arrange the object info boxes-->
    <script type="text/javascript" src="/res/js/pkg/jquery.masonry.min.js"></script>

    <!-- jquery.dataTables for spiffy feature-laden lists -->
    <script type="text/javascript" src="/res/js/pkg/dataTables/media/js/jquery.dataTables.min.js"></script>
    <link rel="stylesheet" href="/res/css/dataTables.css" type="text/css" media="screen, projection"/>

    <script type="text/javascript">
      <![CDATA[
               (function($) {
               var cache = [];
               // Arguments are image paths relative to the current page.
               $.preLoadImages = function() {
               var args_len = arguments.length;
               for (var i = args_len; i--;) {
               var cacheImage = document.createElement('img');
               cacheImage.src = arguments[i];
               cache.push(cacheImage);
               }
               }
$(document).data('updatedOn', new Date(]]><xsl:copy-of select="$currentTime"/><![CDATA[));
})(jQuery)

]]>
    </script>

    <!-- fire up spiffy UI scripts-->
    <script type="text/javascript" src="/res/js/app/ui-init.js"></script>

  </xsl:template>

  <!-- page header for apps -->
  <xsl:template name="app_header">
    <xsl:param name="app_name"/>
    <xsl:param name="icon"/>
    <span id="authUser" style="visibility:hidden">-</span>

    <xsl:comment>template: /html/common_includes/components.xsl:app_header</xsl:comment>

    <div class="header rounded-top gradient-grey">
      <div class="container">
        <div><xsl:attribute name="class">title span-24 last <xsl:copy-of select="$icon"/></xsl:attribute>
        <h1><xsl:value-of select="$app_name"/></h1>
        </div>
      </div>
    </div>
  </xsl:template>

  <!-- page header for views -->
  <xsl:template name="view_header">
    <xsl:param name="label_name"/>
    <xsl:param name="display_name"/>
    <xsl:param name="perspective" select="$currentPerspective"/>
    <xsl:param name="icon"/>

    <xsl:comment>template: /html/common_includes/components.xsl:view_header</xsl:comment>

    <div class="header rounded-bottom gradient-grey shadow">
      <div class="container">
        <div><xsl:attribute name="class">title span-24 last <xsl:copy-of select="$icon"/></xsl:attribute>
        <h1>
          <xsl:value-of select="$label_name"/><xsl:text> </xsl:text>
          <xsl:value-of select="$display_name"/>
        </h1>
        </div>
      </div>
    </div>

  </xsl:template>

  <!-- page header for sets -->
  <xsl:template name="set_header">
    <xsl:param name="display_name"/>

    <xsl:comment>template: /html/common_includes/components.xsl:view_header</xsl:comment>

    <div class="header rounded-bottom gradient-grey shadow">
      <div class="container">
        <div><xsl:attribute name="class">title span-24 last</xsl:attribute>
        <h1 class="no_icon" style="margin-left: 0;"><xsl:value-of select="$display_name"/></h1>
        </div>
      </div>
    </div>

  </xsl:template>

  <!-- app control bar  -->
  <xsl:template name="control_bar_app">
    <xsl:comment>template: /html/common_includes/components.xsl name:control_bar_app</xsl:comment>

    <div class="control_bar app rounded-bottom shadow">
      <div class="control_bar_menu" id="bar_menu">
        <xsl:call-template name="control_bar_cache_info"/>
        <xsl:call-template name="control_bar_menu"/>
      </div>

      <div class="control_bar_base" id="bar_base">&#160;</div>

    </div>
  </xsl:template>

  <!-- view control bar  -->
  <xsl:template name="control_bar_view">
    <xsl:comment>template: /html/common_includes/components.xsl name:control_bar_view</xsl:comment>
    <div class="control_bar view shadow">
      <div class="control_bar_menu" id="bar_menu">
        <xsl:call-template name="control_bar_cache_info"/>
        <xsl:call-template name="control_bar_menu"/>
      </div>

      <div class="control_bar_base" id="bar_base">&#160;</div>

    </div>
  </xsl:template>

  <!-- application menu for control bars -->
  <xsl:template name="control_bar_menu">
    <xsl:comment>template: /html/common_includes/components.xsl name:control_bar_menu</xsl:comment>

    <ul class="app_menu" style="width: 450px">
      <li style="width: 50px; float: right">
        <a id="authUserLink"><div id="authUser" style="color: white">-</div></a>
      </li>
    <xsl:choose>
        <xsl:when test="$objectClassName='Genome::Search'">
        </xsl:when>
        <xsl:when test="$objectClassName='Genome::Search::Query'">
            <li style="width: 30px; float: right">
                <a id="home_link" href="/view/genome/search/status.html">Main</a>
            </li>
        </xsl:when>
        <xsl:otherwise><xsl:comment>items in control bar when not the search main page</xsl:comment>
            <li style="width: 210px; float: right">
            <form id="search" action="/view/genome/search/query/status.html">
                <input name="query"/>
                <a id="search_submit" style="padding: 0px 5px 0px 5px" class="app btn shadow" onClick="javascript:$('#search').submit();">
                Search</a>
            </form>
            </li>

            <li style="width: 30px; float: right">
                <a id="home_link" href="/view/genome/search/status.html">Main</a>
            </li>

            <li style="width: 70px; float: right">
              <a id="perspective_switcher" title="Available Perspectives">
                <xsl:attribute name="alt"><xsl:value-of select="$objectId"/></xsl:attribute>
                <xsl:attribute name="href"><xsl:text>/viewajax/u-r/object/type/available-views.html?class_name=</xsl:text><xsl:value-of select="$objectClassName"/></xsl:attribute>
              Perspectives
              </a>
            </li>
        </xsl:otherwise>
    </xsl:choose>
    </ul>
  </xsl:template>

  <!-- cache info block for control bars -->
  <xsl:template name="control_bar_cache_info">
    <xsl:comment>template: /html/common_includes/components.xsl name:control_bar_cache_info</xsl:comment>

    <div class="cache_info">
      <div class="titlebar">
      </div>
      <div class="cache_time">
        <p>
          View generated<br/>
          <strong><span id="updatedOn"></span></strong>
        </p>
      </div>
      <div class="cache_refresh">
        <a class="btn_cache_refresh rounded" id="refreshCache" title="Refresh Cache"></a>
      </div>
    </div>
  </xsl:template>

  <!-- basic footer -->
  <xsl:template name="footer">
    <xsl:param name="footer_text"/>

    <xsl:comment>template: /html/common_includes/components.xsl:footer</xsl:comment>

    <div class="footer rounded shadow span-24 last gradient-grey">
      <div class="container">
        <p class="small"><xsl:copy-of select="$footer_text"/></p>
      </div>
    </div>
    <br/>
    <br/>
  </xsl:template>

</xsl:stylesheet>
