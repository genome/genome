<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
xmlns:rest="urn:rest">

  <xsl:template name="ur_object" match="object">
    <xsl:comment>template: /html/search-result/ur_object.xsl; match="object" </xsl:comment>

    <div class="result">
    <table width="100%" cellpadding="0" cellspacing="0" border="0" class="result"><tbody><tr>
      <td>
        <div class="icon">
          <xsl:call-template name="object_link">
            <xsl:with-param name="linktext">
              <img width="32" height="32" src="/res/old/report_resources/apipe_dashboard/images/icons/eye_16.png" />
            </xsl:with-param>
          </xsl:call-template>
        </div>
      </td><td>
        <div class="description">
        <h2 class="name">
          <span class="label">
            <xsl:value-of select="label_name"/>:
          </span>
          <span class="title">
            <xsl:call-template name="object_link">
              <xsl:with-param name="linktext" select="display_name"/>
            </xsl:call-template>
          </span>
        </h2>
      </div>
      </td></tr></tbody></table>
    </div>
  </xsl:template>

</xsl:stylesheet>
