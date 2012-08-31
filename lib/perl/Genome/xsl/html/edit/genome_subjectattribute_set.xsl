<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome_subjectattribute_set" match="object[./types[./isa[@type='Genome::SubjectAttribute::Set']]]">

    <script type="text/javascript" src="/res/js/pkg/dataTables/media/js/jquery.dataTables.js"></script>

    <script type="text/javascript" language="javascript"></script>

    <xsl:call-template name="control_bar_view"/>

    <xsl:call-template name="view_header">
      <xsl:with-param name="label_name" select="'Sample Attributes '" />
      <xsl:with-param name="display_name" select="'test_sample0'" />
      <xsl:with-param name="icon" select="'genome_sample_32'" />
    </xsl:call-template>

<!--      <xsl:with-param name="display_name" select="aspect[@name='name']/value" />
-->

    <div class="content rounded shadow">
      <div class="container">

            <!-- details for this sample -->

                <fieldset>
                    <legend>Nomenclature: myNomenclature</legend>

                    <xsl:for-each select="/object/aspect[@name='members']/object[1]/aspect[@name='all_nomenclature_fields']/object">
                        <xsl:variable name='field_name' select="aspect[@name='name']/value"/>
                        <xsl:variable name='field_type' select="aspect[@name='type']/value"/>
                        <xsl:variable name='field_value' select="/object/aspect[@name='members']/object/aspect[value = $field_name]/../aspect[@name='attribute_value']/value"/>

                            <p>
                            <label> <xsl:value-of select="$field_name"/>:</label><br/>

                            <xsl:choose>
                                <xsl:when test="$field_type = 'enumerated'">
                                    <select>
                                        <xsl:for-each select="aspect[@name='enumerated_values']/object">
                                            <option> 
                                                <xsl:if test="$field_value = aspect[@name='value']/value"> 
                                                    <xsl:attribute name="selected">true</xsl:attribute> 
                                                </xsl:if>
                                                <xsl:value-of select="aspect[@name='value']/value"/> 
                                            </option>
                                        </xsl:for-each>
                                    </select>
                                </xsl:when> 
                                <xsl:otherwise>
                                    <input class="text"> <xsl:attribute name="value"><xsl:value-of select="$field_value"/> </xsl:attribute> </input>
                                </xsl:otherwise>
                            </xsl:choose>

                            </p>

                    </xsl:for-each>

                    <p><button>Save</button></p>
                </fieldset>
      </div> <!-- end container -->
    </div> <!-- end content -->

    <xsl:call-template name="footer">
      <xsl:with-param name="footer_text">
        <br/>
      </xsl:with-param>
    </xsl:call-template>

  </xsl:template>

</xsl:stylesheet>



