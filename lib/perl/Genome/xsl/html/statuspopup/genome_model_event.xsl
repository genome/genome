<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
xmlns:rest="urn:rest">

  <xsl:template name="genome_model_event" match="genome-model-event">
     <table class="boxy_info" cellpadding="0" cellspacing="0" border="0" width="300"><tbody>
       <xsl:for-each select="child::*">
         <xsl:choose>
           <!-- handle output and error log file nodes -->
           <xsl:when test="contains(local-name(),'_file')">
             <xsl:call-template name="genome_model_event_row">
               <xsl:with-param name="label">
                 <xsl:value-of select="local-name()"/>
               </xsl:with-param>
               <xsl:with-param name="value">
	             <a><xsl:attribute name="href"><xsl:text>http://gscweb</xsl:text><xsl:value-of select="current()"/></xsl:attribute>
	               <xsl:call-template name="substring-after-last">
	                 <xsl:with-param name="input" select="current()"/>
	                 <xsl:with-param name="substr" select="'/'"/>
	               </xsl:call-template>
	             </a>
               </xsl:with-param>>
             </xsl:call-template>
           </xsl:when>

           <!-- handle alignment_directory node(s) -->
           <xsl:when test="contains(local-name(),'alignment_directory')">
             <xsl:choose>
               <xsl:when test="starts-with(current(), '/')">
                 <!-- starts with a /, so most likely is a directory string -->
                 <xsl:call-template name="genome_model_event_row">
                   <xsl:with-param name="label">
                     <xsl:value-of select="local-name()"/>
                   </xsl:with-param>
                   <xsl:with-param name="value">
                     <a><xsl:attribute name="href"><xsl:text>http://gscweb</xsl:text><xsl:value-of select="current()"/></xsl:attribute>
                       <xsl:call-template name="substring-after-last">
                         <xsl:with-param name="input" select="current()"/>
                         <xsl:with-param name="substr" select="'/'"/>
                       </xsl:call-template>
                     </a>
                   </xsl:with-param>
                 </xsl:call-template>
               </xsl:when>
               <xsl:otherwise>
                 <!-- doesn't start with a /, so probably is a message of some kind -->
                 <xsl:call-template name="genome_model_event_row">
                   <xsl:with-param name="label">
                     <xsl:value-of select="local-name()"/>
                   </xsl:with-param>
                   <xsl:with-param name="value">
                     <xsl:value-of select="normalize-space(current())"/>
                   </xsl:with-param>
                 </xsl:call-template>
               </xsl:otherwise>
             </xsl:choose>
           </xsl:when>

           <xsl:otherwise>
             <xsl:call-template name="genome_model_event_row">
               <xsl:with-param name="label">
                 <xsl:value-of select="local-name()"/>
               </xsl:with-param>
               <xsl:with-param name="value">
                 <xsl:value-of select="current()"/>
               </xsl:with-param>
             </xsl:call-template>
           </xsl:otherwise>
         </xsl:choose>
       </xsl:for-each>
       <!-- see if we have any instrument data and append that to the event object if we do -->
       <xsl:if test="instrument_data_id!=''">
         <xsl:variable name="inst_data_id" select="instrument_data_id" />
         <xsl:for-each select="//instrument_data[@id=$inst_data_id]/*" >
           <xsl:choose>
             <xsl:when test="local-name() != 'gerald_directory' and local-name() != 'lane'">
               <xsl:call-template name="genome_model_event_row">
                 <xsl:with-param name="label">
                   <xsl:value-of select="local-name()"/>
                 </xsl:with-param>
                 <xsl:with-param name="value">
                    <xsl:value-of select="current()"/>
                 </xsl:with-param>
               </xsl:call-template>
             </xsl:when>
           </xsl:choose>
         </xsl:for-each>
       </xsl:if>
    </tbody></table>
  </xsl:template>
  
  <xsl:template name="genome_model_event_row">
    <xsl:param name="label" />
    <xsl:param name="value" />
    <tr>
      <td class="label">
        <xsl:copy-of select="$label"/>
      </td>
      <td class="value">
        <xsl:copy-of select="$value"/>
      </td>
    </tr>
  </xsl:template>
</xsl:stylesheet>
