<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
xmlns:rest="urn:rest">

  <xsl:template name="workflow_operation_instance" match="object[./types[./isa[@type='Workflow::Operation::Instance']]]">

    <xsl:if test="count(ancestor::aspect) = 0">
      <script type='text/javascript' src='/res/js/boxy/javascripts/jquery.boxy.js'></script>
      <link rel="stylesheet" href="/res/js/boxy/stylesheets/boxy.css" type="text/css" />
      <script type="text/javascript">
        <![CDATA[
          function status_popup(title,typeurl,id) {
            Boxy.load("/viewajax/" + typeurl + "/statuspopup.html?id=" + id, {cache: true, title: title, afterShow: function() { this.center; this.resize(400,400); this.getContent().css('overflow','auto'); } });
            return false;
          }
        ]]>
      </script>

      <table class="lister" border="0" width="100%" cellspacing="0" cellpadding="0">
        <colgroup>
          <col/>
          <col width="40%"/>
          <col/>
          <col/>
          <col/>
          <col/>
        </colgroup>
        <thead>
          <tr>
            <th>idx</th>
            <th>operation</th>
            <th>status</th>
            <th>started</th>
            <th>ended</th>
            <th class="last">elapsed</th>
          </tr>
        </thead>
        <tbody>
          <xsl:call-template name="inner_woi"/>
        </tbody>
      </table>

      <img><xsl:attribute name="src"><xsl:text>/view/workflow/cache/graph.png?id=</xsl:text><xsl:value-of select="aspect[@name='cache_workflow_id']/value"/></xsl:attribute></img>
    </xsl:if>
    <xsl:if test="count(ancestor::aspect) > 0">
      <xsl:call-template name="inner_woi"/>
    </xsl:if>

  </xsl:template>

  <xsl:template name="inner_woi">

    <xsl:if test="count(aspect[@name='operation_type']/object/types/isa[@type='Workflow::OperationType::Command']) + count(aspect[@name='operation_type']/object/types/isa[@type='Workflow::OperationType::Event']) > 0">

      <xsl:variable name="currentLink">
        <xsl:value-of select="rest:typetourl(aspect[@name='current']/object[1]/@type)" />
      </xsl:variable>

      <tr onmouseover="this.className = 'hover'" onmouseout="this.className=''">
        <xsl:attribute name="onclick">
          <xsl:text>javascript:status_popup('</xsl:text><xsl:value-of select="aspect[@name='name']/value"/><xsl:text>','</xsl:text><xsl:value-of select="$currentLink"/><xsl:text>','</xsl:text><xsl:value-of select="aspect[@name='current']/object[1]/@id"/><xsl:text>');</xsl:text>
        </xsl:attribute>
        <xsl:for-each select="aspect[@name='parallel_index']">
          <td><xsl:value-of select="value"/></td>
        </xsl:for-each>
        <xsl:for-each select="aspect[@name='name']">
          <td><xsl:value-of select="value"/></td>
        </xsl:for-each>
        <xsl:for-each select="aspect[@name='status']">
          <td><xsl:attribute name="class"><xsl:text>status </xsl:text><xsl:value-of select="value"/></xsl:attribute><xsl:value-of select="value"/></td>
        </xsl:for-each>
        <xsl:for-each select="aspect[@name='start_time']">
          <td><xsl:value-of select="value"/></td>
        </xsl:for-each>
        <xsl:for-each select="aspect[@name='end_time']">
          <td><xsl:value-of select="value"/></td>
        </xsl:for-each>
        <xsl:for-each select="aspect[@name='elapsed_time']">
          <td class="last"><xsl:value-of select="value"/></td>
        </xsl:for-each>
      </tr>
    </xsl:if>

    <xsl:if test="count(aspect[@name='related_instances']) > 0">
      <xsl:for-each select="aspect[@name='related_instances']">
        <xsl:apply-templates/>
      </xsl:for-each>
    </xsl:if>


  </xsl:template>

</xsl:stylesheet>

