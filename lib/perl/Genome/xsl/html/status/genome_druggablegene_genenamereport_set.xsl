<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template name="gene_go_results_table">
  <xsl:param name="name"/>
  <xsl:param name="definite_go_results"/>

  <div class="box_header span-24 last rounded-top">
    <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last"><xsl:value-of select="$name"/></h3></div>
  </div>

  <div class="box_content rounded-bottom span-24 last">
    <table class='dataTable' id='interactions'>
        <thead>
          <tr>
            <th>Gene</th>
            <th>GO Category</th>
            <th>Search Term</th>
          </tr>
        </thead>
        <tbody>
          <xsl:for-each select="$definite_go_results/item">
            <tr>
              <th><xsl:value-of select="gene_group_name"/></th>
              <th><xsl:value-of select="category_name"/></th>
              <th><xsl:value-of select="search_terms"/></th>
            </tr>
          </xsl:for-each>
        </tbody>
    </table>
  </div>

</xsl:template>

<xsl:template name="gene_go_summary_table">
  <xsl:param name="name"/>
  <xsl:param name="go_results_summary"/>

  <div class="box_header span-24 last rounded-top">
    <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last"><xsl:value-of select="$name"/></h3></div>
  </div>

  <div class="box_content rounded-bottom span-24 last">
    <table class='dataTable' id='interactions'>
        <thead>
          <tr>
            <th>GO Category</th>
            <th>Gene Count</th>
            <th>Genes</th>
          </tr>
        </thead>
        <tbody>
          <xsl:for-each select="$go_results_summary/item">
            <tr>
              <th><xsl:value-of select="go_category_name"/></th>
              <th><xsl:value-of select="gene_count"/></th>
              <th><xsl:value-of select="gene_names"/></th>
            </tr>
          </xsl:for-each>
        </tbody>
    </table>
  </div>

</xsl:template>

<xsl:template name="drug_gene_interactions_table">
  <xsl:param name="name"/>
  <xsl:param name="interactions"/>

  <div class="box_header span-24 last rounded-top">
    <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last"><xsl:value-of select="$name"/></h3></div>
  </div>

  <div class="box_content rounded-bottom span-24 last">
    <table class='dataTable' id='interactions'>
      <thead>
        <tr>
          <th>Source</th>
          <th>Drug</th>
          <th>Interaction Type</th>
          <th>Gene</th>
          <th>Search Term</th>
          <xsl:if test="$interactions/item/number_of_matches">
            <th>Matches</th>
          </xsl:if>
        </tr>
      </thead>
      <tbody>
        <xsl:for-each select="$interactions/item">
          <tr>
            <th><a class='mini btn'>
                <xsl:attribute name="href">
                  <xsl:value-of select='source_url'/>
                </xsl:attribute>
                <xsl:value-of select='source' />
            </a></th>
            <th>
              <xsl:call-template name='object_link_button'>
                <xsl:with-param name='type' select="'Genome::DruggableGene::DrugNameReport::Set'"/>
                <xsl:with-param name="key" select="'name'"/>
                <xsl:with-param name="id" select="drug"/>
                <xsl:with-param name="linktext">
                  <xsl:value-of select="substring(human_readable_drug_name, 1, 40)"/>
                  <xsl:if test="string-length(human_readable_drug_name) &gt; 40">
                    <xsl:text>...</xsl:text>
                  </xsl:if>
                </xsl:with-param>
              </xsl:call-template>
            </th>
            <th>
              <xsl:call-template name='object_link_button'>
                <xsl:with-param name='type' select="'Genome::DruggableGene::DrugGeneInteractionReport::Set'"/>
                <xsl:with-param name="keys" select='.'/>
                <xsl:with-param name="linktext">
                  <xsl:choose>
                    <xsl:when test="interaction_type = 'na'">
                      <xsl:text>n/a</xsl:text>
                    </xsl:when>
                    <xsl:otherwise>
                      <xsl:value-of select="interaction_type"/>
                    </xsl:otherwise>
                  </xsl:choose>
                </xsl:with-param>
              </xsl:call-template>
            </th>
            <th>
              <xsl:choose>
                <xsl:when test="group">
                  <xsl:call-template name='object_link_button'>
                    <xsl:with-param name='type' select="'Genome::DruggableGene::GeneNameGroup'"/>
                    <xsl:with-param name="key" select="'name'"/>
                    <xsl:with-param name="id" select="group"/>
                    <xsl:with-param name="linktext" select="group"/>
                  </xsl:call-template>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:call-template name='object_link_button'>
                    <xsl:with-param name='type' select="'Genome::DruggableGene::GeneNameReport::Set'"/>
                    <xsl:with-param name="key" select="'name'"/>
                    <xsl:with-param name="id" select="gene"/>
                    <xsl:with-param name="linktext" select="group"/>
                  </xsl:call-template>
                </xsl:otherwise>
              </xsl:choose>
            </th>
            <th><xsl:value-of select='search_terms' /></th>
            <xsl:if test="number_of_matches">
              <th><xsl:value-of select="number_of_matches"/></th>
            </xsl:if>
          </tr>
        </xsl:for-each>
      </tbody>
    </table>
  </div>
</xsl:template>

<xsl:template name="go_results" match="go_results">
  <script type='text/javascript' src='/res/js/pkg/boxy/javascripts/jquery.boxy.js'></script>
  <link rel="stylesheet" href="/res/js/pkg/boxy/stylesheets/boxy.css" type="text/css" />
  <link rel="stylesheet" href="/res/css/genome_druggablegene_genenamereport_set.css" type="text/css" />
  <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>
  <script class="jsbin" src="http://datatables.net/download/build/jquery.dataTables.nightly.js"></script>

  <xsl:call-template name="control_bar_view"/>

  <xsl:call-template name="view_header">
    <xsl:with-param name="label_name" select="'The Drug-Gene Interactions Database Potentially Druggable Gene Families'" />
    <xsl:with-param name="icon" select="'genome_genename_32'" />
  </xsl:call-template>

  <div class="content rounded shadow">
    <div class="container">
      <div class="span_24_box_masonry">
        <xsl:call-template name="gene_go_results_table">
          <xsl:with-param name="name" select="'Definite GO Results'" />
          <xsl:with-param name="definite_go_results" select="definite_go_results" />
        </xsl:call-template>
      </div>
      <div class="span_24_box_masonry">
        <xsl:call-template name="gene_go_results_table">
          <xsl:with-param name="name" select="'Filtered Definite GO Results'" />
          <xsl:with-param name="definite_go_results" select="filtered_definite_go_results" />
        </xsl:call-template>
      </div>

      <div class="span_24_box_masonry">
        <xsl:call-template name="gene_go_summary_table">
          <xsl:with-param name="name" select="'GO Results Summary'" />
          <xsl:with-param name="go_results_summary" select="go_results_summary" />
        </xsl:call-template>
      </div>
    </div>
  </div>

  <script type="text/javascript">
    $(document).ready(function(){
    $('.dataTable').dataTable({
    "sScrollX": "100%",
    "sScrollInner": "110%",
    "bJQueryUI": true,
    "sPaginationType": "full_numbers",
    "bStateSave": true,
    "iDisplayLength": 25,
    'oLanguage': { 'sSearch': 'Filter results:' },
    });
    }
    );
  </script>

  <xsl:call-template name="footer">
    <xsl:with-param name="footer_text">
      <FORM><INPUT TYPE="button" VALUE="Back" onClick="history.go(-1);return true;" /></FORM>
    </xsl:with-param>
  </xsl:call-template>

</xsl:template>

<xsl:template name="drug_gene_interaction" match="drug_gene_interaction">
  <script type='text/javascript' src='/res/js/pkg/boxy/javascripts/jquery.boxy.js'></script>
  <link rel="stylesheet" href="/res/js/pkg/boxy/stylesheets/boxy.css" type="text/css" />
  <link rel="stylesheet" href="/res/css/genome_druggablegene_genenamereport_set.css" type="text/css" />
  <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>
  <script class="jsbin" src="http://datatables.net/download/build/jquery.dataTables.nightly.js"></script>

  <xsl:call-template name="control_bar_view"/>

  <xsl:call-template name="view_header">
    <xsl:with-param name="label_name" select="'The Drug-Gene Interactions Database'" />
    <xsl:with-param name="icon" select="'genome_genename_32'" />
  </xsl:call-template>

  <div class="content rounded shadow">
    <div class="container">

      <div class="span_24_box_masonry">
        <xsl:call-template name="drug_gene_interactions_table">
          <xsl:with-param name="interactions" select="interactions" />
          <xsl:with-param name="name" select="'Interactions'" />
        </xsl:call-template>
        <xsl:call-template name="drug_gene_interactions_table">
          <xsl:with-param name="interactions" select="ambiguous_interactions" />
          <xsl:with-param name="name" select="'Ambiguously Matched Interactions'" />
        </xsl:call-template>
        <xsl:call-template name="drug_gene_interactions_table">
          <xsl:with-param name="interactions" select="filtered_interactions" />
          <xsl:with-param name="name" select="'Filtered Interactions'" />
        </xsl:call-template>
        <xsl:call-template name="drug_gene_interactions_table">
          <xsl:with-param name="interactions" select="filtered_ambiguous_interactions" />
          <xsl:with-param name="name" select="'Filtered Ambiguously Matched Interactions'" />
        </xsl:call-template>
      </div>

      <div id='objects' class='span-24 last'>

        <div class="span_12_box_masonry">
          <div class="box_header span-12 last rounded-top">
            <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last">Search Terms Without Matches</h3></div>
          </div>
          <div class="box_content rounded-bottom span-12 last">
            <br />
            <ul>
              <xsl:for-each select="search_terms_without_groups/item">
                <li><xsl:value-of select='.' /></li>
              </xsl:for-each>
            </ul>
          </div>
        </div>

        <div class="span_12_box_masonry">
          <div class="box_header span-12 last rounded-top">
            <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last">Genes Without Interactions</h3></div>
          </div>

          <div class="box_content rounded-bottom span-12 last">
            <br />
            <ul>
              <xsl:for-each select="missing_interactions/item">
                <li>
                  <xsl:call-template name='object_link_button'>
                    <xsl:with-param name='type' select="'Genome::DruggableGene::GeneNameGroup'"/>
                    <xsl:with-param name="key" select="'name'"/>
                    <xsl:with-param name="id" select="group"/>
                    <xsl:with-param name="linktext" select="group"/>
                  </xsl:call-template>
                  <xsl:text> for search term </xsl:text>
                  <xsl:value-of select='search_terms' />
                </li>
              </xsl:for-each>
            </ul>
          </div>
        </div>

        <div class="span_12_box_masonry">
          <div class="box_header span-12 last rounded-top">
            <div class="box_title"><h3 class="genome_genenamereport_16 span-7 last">Ambiguously Matched Genes Without Interactions</h3></div>
          </div>

          <div class="box_content rounded-bottom span-12 last">
            <br />
            <ul>
              <xsl:for-each select="missing_ambiguous_interactions/item">
                <li>
                  <xsl:call-template name='object_link_button'>
                    <xsl:with-param name='type' select="'Genome::DruggableGene::GeneNameGroup'"/>
                    <xsl:with-param name="key" select="'name'"/>
                    <xsl:with-param name="id" select="group"/>
                    <xsl:with-param name="linktext" select="group"/>
                  </xsl:call-template>
                  <xsl:text> for search term </xsl:text>
                  <xsl:value-of select='search_terms' />
                  <xsl:text> with quantity of matches </xsl:text>
                  <xsl:value-of select='number_of_matches' />
                </li>
              </xsl:for-each>
            </ul>
          </div>
        </div>

      </div>

    </div> <!-- end container -->
  </div> <!-- end content -->


  <script type="text/javascript">
    $(document).ready(function(){
    $('.dataTable').dataTable({
    "sScrollX": "100%",
    "sScrollInner": "110%",
    "bJQueryUI": true,
    "sPaginationType": "full_numbers",
    "bStateSave": true,
    "iDisplayLength": 25,
    'oLanguage': { 'sSearch': 'Filter results:' },
    });
    }
    );
  </script>

  <xsl:call-template name="footer">
    <xsl:with-param name="footer_text">
      <FORM><INPUT TYPE="button" VALUE="Back" onClick="history.go(-1);return true;" /></FORM>
    </xsl:with-param>
  </xsl:call-template>

</xsl:template>

<xsl:template name="genome_druggablegene_genenamereport_set" match="object[./types[./isa[@type='Genome::DruggableGene::GeneNameReport::Set']]]">

  <script type='text/javascript' src='/res/js/pkg/boxy/javascripts/jquery.boxy.js'></script>
  <link rel="stylesheet" href="/res/js/pkg/boxy/stylesheets/boxy.css" type="text/css" />
  <script type='text/javascript' src='/res/js/app/genome_model_build_list.js'></script>

  <xsl:call-template name="control_bar_view"/>

  <xsl:variable name='header_name'>
    <xsl:choose>
      <xsl:when test='count(aspect[@name="name"]/value)=1'>
        <xsl:value-of select='aspect[@name="name"]/value'/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select='./display_name'/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <xsl:call-template name="view_header">
    <xsl:with-param name="label_name" select="'genename:'" />
    <xsl:with-param name="display_name" select='$header_name'/>
    <xsl:with-param name="icon" select="'genome_genename_32'" />
  </xsl:call-template>

  <div class="content rounded shadow">
    <div class="container">
      <xsl:for-each select="aspect[@name='members']/object">
        <xsl:call-template name="genome_genenamereport_box"/>
      </xsl:for-each>
    </div> <!-- end container -->
  </div> <!-- end content -->

  <xsl:call-template name="footer">
    <xsl:with-param name="footer_text">
      <br/>
    </xsl:with-param>
  </xsl:call-template>

</xsl:template>

</xsl:stylesheet>
