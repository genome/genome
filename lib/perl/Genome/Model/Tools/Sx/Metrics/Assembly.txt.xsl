<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
  <xsl:variable name="contigs_length" select="//aspect[@name='contigs_length']/value"/>
  <xsl:variable name="reads_processed" select="//aspect[@name='reads_processed']/value"/>
  <xsl:variable name="reads_processed_length_q20" select="//aspect[@name='reads_processed_length_q20']/value"/>
  <xsl:variable name='reads_assembled' select="//aspect[@name='reads_assembled']/value"/>
  <xsl:variable name='reads_assembled_unique' select="//aspect[@name='reads_assembled_unique']/value"/>
  <xsl:variable name='reads_assembled_duplicate' select="//aspect[@name='reads_assembled_duplicate']/value"/>
  <xsl:variable name='reads_not_assembled' select="$reads_processed - $reads_assembled_unique"/>
  <xsl:variable name='tier_one' select="//aspect[@name='tier_one']/value"/>
  <xsl:variable name='tier_two' select="$tier_one * 2"/>
  <xsl:template match="/">*** SIMPLE READ STATS ***
Total input reads: <xsl:value-of select="$reads_processed"/>
Total input bases: <xsl:value-of select="//aspect[@name='reads_processed_length']/value"/> bp
Total Q20 bases: <xsl:choose><xsl:when test='$reads_processed_length_q20 != 0'><xsl:value-of select="$reads_processed_length_q20"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose> bp
Average Q20 bases per read: <xsl:choose><xsl:when test='$reads_processed_length_q20 != 0'><xsl:value-of select="format-number(($reads_processed_length_q20 div $reads_processed), '#')"/></xsl:when><xsl:otherwise>NaN</xsl:otherwise></xsl:choose> bp
Average read length: <xsl:value-of select="//aspect[@name='reads_processed_average_length']/value"/> bp
Placed reads: <xsl:choose><xsl:when test='$reads_assembled_unique != 0'><xsl:value-of select="$reads_assembled_unique"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>     % of total input reads: <xsl:value-of select="//aspect[@name='reads_assembled_success_percent']/value"/>
  (reads in scaffolds: <xsl:choose><xsl:when test='$reads_assembled != 0'><xsl:value-of select="$reads_assembled"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>)
  (unique reads: <xsl:choose><xsl:when test='$reads_assembled_unique != 0'><xsl:value-of select="$reads_assembled_unique"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>)
  (duplicate reads: <xsl:choose><xsl:when test='$reads_assembled_unique != 0'><xsl:value-of select="$reads_assembled_duplicate"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>)
Unplaced reads: <xsl:choose><xsl:when test='$reads_assembled_unique != 0'><xsl:value-of select="$reads_not_assembled"/></xsl:when><xsl:otherwise>NaN</xsl:otherwise></xsl:choose>     % of total input reads: <xsl:value-of select="//aspect[@name='reads_not_assembled_percent']/value"/>
Chaff rate: <xsl:choose><xsl:when test='$reads_assembled_unique != 0'> <xsl:value-of select="format-number(($reads_not_assembled div $reads_processed), '#.##%')"/></xsl:when><xsl:otherwise>NaN</xsl:otherwise></xsl:choose>
Q20 base redundancy: <xsl:choose><xsl:when test="$reads_processed_length_q20 != '0'"><xsl:value-of select="format-number(($reads_processed_length_q20 div $contigs_length), '#.#')"/></xsl:when><xsl:otherwise>NaN</xsl:otherwise></xsl:choose>X


*** Contiguity: Contig ***
Total Contig number: <xsl:value-of select="//aspect[@name='contigs_count']/value"/>
Total Contig bases: <xsl:value-of select="$contigs_length"/> bp
Average Contig length: <xsl:value-of select="//aspect[@name='contigs_average_length']/value"/> bp
Maximum Contig length: <xsl:value-of select="//aspect[@name='contigs_maximum_length']/value"/> bp
N50 Contig length: <xsl:value-of select="//aspect[@name='contigs_n50_length']/value"/> bp
N50 contig number: <xsl:value-of select="//aspect[@name='contigs_n50_count']/value"/>

Major Contig (>= <xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/> bp) number: <xsl:value-of select="//aspect[@name='contigs_major_count']/value"/>
Major Contig bases: <xsl:value-of select="//aspect[@name='contigs_major_length']/value"/> bp
Major_Contig avg contig length: <xsl:value-of select="//aspect[@name='contigs_major_average_length']/value"/> bp
Major_Contig N50 contig length: <xsl:value-of select="//aspect[@name='contigs_major_n50_length']/value"/> bp
Major_Contig N50 contig number: <xsl:value-of select="//aspect[@name='contigs_major_n50_count']/value"/>

% of assembly in major contig: <xsl:value-of select="//aspect[@name='contigs_major_percent']/value"/>
Placed reads in major contig: <xsl:choose><xsl:when test="//aspect[@name='contigs_major_read_count']/value &gt; 0 or //aspect[@name='contigs_minor_read_count']/value &gt; 0"><xsl:value-of select="//aspect[@name='contigs_major_read_count']/value"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>, <xsl:value-of select="//aspect[@name='contigs_major_read_percent']/value"/>%

Minor Contig (up to <xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/> bp) number: <xsl:value-of select="//aspect[@name='contigs_minor_count']/value"/>
Minor Contig bases: <xsl:value-of select="//aspect[@name='contigs_minor_length']/value"/> bp
Minor Contig avg contig length: <xsl:value-of select="//aspect[@name='contigs_minor_average_length']/value"/>
Minor Contig N50 contig length: <xsl:value-of select="//aspect[@name='contigs_minor_n50_length']/value"/>
Minor Contig N50 contig number: <xsl:value-of select="//aspect[@name='contigs_minor_n50_count']/value"/>

Placed reads in Minor contig: <xsl:choose><xsl:when test="//aspect[@name='contigs_major_read_count']/value &gt; 0 or //aspect[@name='contigs_minor_read_count']/value &gt; 0"><xsl:value-of select="//aspect[@name='contigs_minor_read_count']/value"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>, <xsl:value-of select="//aspect[@name='contigs_minor_read_percent']/value"/>%

Top tier (up to <xsl:value-of select="$tier_one"/> bp): 
  Contig number: <xsl:value-of select="//aspect[@name='contigs_t1_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='contigs_t1_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='contigs_t1_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='contigs_t1_length']/value"/> bp
  Top tier N50 contig length: <xsl:value-of select="//aspect[@name='contigs_t1_n50_length']/value"/> bp
  Top tier N50 contig number: <xsl:value-of select="//aspect[@name='contigs_t1_n50_count']/value"/>
Middle tier (<xsl:value-of select="$tier_one"/> bp -- <xsl:value-of select="$tier_two"/> bp): 
  Contig number: <xsl:value-of select="//aspect[@name='contigs_t2_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='contigs_t2_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='contigs_t2_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='contigs_t2_length']/value"/> bp
  Middle tier N50 contig length: <xsl:value-of select="//aspect[@name='contigs_t2_n50_length']/value"/> bp
  Middle tier N50 contig number: <xsl:value-of select="//aspect[@name='contigs_t2_n50_count']/value"/>
Bottom tier (<xsl:value-of select="$tier_two"/> bp -- end): 
  Contig number: <xsl:value-of select="//aspect[@name='contigs_t3_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='contigs_t3_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='contigs_t3_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='contigs_t3_length']/value"/> bp
  Bottom tier N50 contig length: <xsl:value-of select="//aspect[@name='contigs_t3_n50_length']/value"/> bp
  Bottom tier N50 contig number: <xsl:value-of select="//aspect[@name='contigs_t3_n50_count']/value"/>


*** Contiguity: Supercontig ***
Total Supercontig number: <xsl:value-of select="//aspect[@name='supercontigs_count']/value"/>
Total Supercontig bases: <xsl:value-of select="//aspect[@name='supercontigs_length']/value"/> bp
Average Supercontig length: <xsl:value-of select="//aspect[@name='supercontigs_average_length']/value"/> bp
Maximum Supercontig length: <xsl:value-of select="//aspect[@name='supercontigs_maximum_length']/value"/> bp
N50 Supercontig length: <xsl:value-of select="//aspect[@name='supercontigs_n50_length']/value"/> bp
N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_n50_count']/value"/>

Major Supercontig (>= <xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/> bp) number: <xsl:value-of select="//aspect[@name='supercontigs_major_count']/value"/>
Major Supercontig bases: <xsl:value-of select="//aspect[@name='supercontigs_major_length']/value"/> bp
Major_Supercontig avg contig length: <xsl:value-of select="//aspect[@name='supercontigs_major_average_length']/value"/> bp
Major_Supercontig N50 contig length: <xsl:value-of select="//aspect[@name='supercontigs_major_n50_length']/value"/> bp
Major_Supercontig N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_major_n50_count']/value"/>

% of assembly in major supercontig: <xsl:value-of select="//aspect[@name='supercontigs_major_percent']/value"/>
Placed reads in major supercontig: <xsl:choose><xsl:when test="//aspect[@name='supercontigs_major_read_count']/value &gt; 0 or //aspect[@name='supercontigs_minor_read_count']/value &gt; 0"><xsl:value-of select="//aspect[@name='supercontigs_major_read_count']/value"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>, <xsl:value-of select="//aspect[@name='supercontigs_major_read_percent']/value"/>%

Minor Supercontig (up to <xsl:value-of select="//aspect[@name='major_contig_threshold']/value"/> bp) number: <xsl:value-of select="//aspect[@name='supercontigs_minor_count']/value"/>
Minor Supercontig bases: <xsl:value-of select="//aspect[@name='supercontigs_minor_length']/value"/> bp
Minor Supercontig avg contig length: <xsl:value-of select="//aspect[@name='supercontigs_minor_average_length']/value"/>
Minor Supercontig N50 contig length: <xsl:value-of select="//aspect[@name='supercontigs_minor_n50_length']/value"/>
Minor Supercontig N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_minor_n50_count']/value"/>

Placed reads in minor supercontig: <xsl:choose><xsl:when test="//aspect[@name='supercontigs_major_read_count']/value &gt; 0 or //aspect[@name='supercontigs_minor_read_count']/value &gt; 0"><xsl:value-of select="//aspect[@name='supercontigs_minor_read_count']/value"/></xsl:when><xsl:otherwise>NA</xsl:otherwise></xsl:choose>, <xsl:value-of select="//aspect[@name='supercontigs_minor_read_percent']/value"/>%

Scaffolds > 1M: <xsl:value-of select="//aspect[@name='scaffolds_1M']/value"/>
Scaffold 250K--1M: <xsl:value-of select="//aspect[@name='scaffolds_250K_1M']/value"/>
Scaffold 100K--250K: <xsl:value-of select="//aspect[@name='scaffolds_100K_250K']/value"/>
Scaffold 10--100K: <xsl:value-of select="//aspect[@name='scaffolds_10K_100K']/value"/>
Scaffold 5--10K: <xsl:value-of select="//aspect[@name='scaffolds_5K_10K']/value"/>
Scaffold 2--5K: <xsl:value-of select="//aspect[@name='scaffolds_2K_5K']/value"/>
Scaffold 0--2K: <xsl:value-of select="//aspect[@name='scaffolds_0K_2K']/value"/>

Top tier (up to <xsl:value-of select="$tier_one"/> bp): 
  Contig number: <xsl:value-of select="//aspect[@name='supercontigs_t1_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='supercontigs_t1_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='supercontigs_t1_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='supercontigs_t1_length']/value"/> bp
  Top tier N50 contig length: <xsl:value-of select="//aspect[@name='supercontigs_t1_n50_length']/value"/> bp
  Top tier N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_t1_n50_count']/value"/>
Middle tier (<xsl:value-of select="$tier_one"/> bp -- <xsl:value-of select="$tier_two"/> bp): 
  Contig number: <xsl:value-of select="//aspect[@name='supercontigs_t2_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='supercontigs_t2_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='supercontigs_t2_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='supercontigs_t2_length']/value"/> bp
  Middle tier N50 contig length: <xsl:value-of select="//aspect[@name='supercontigs_t2_n50_length']/value"/> bp
  Middle tier N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_t2_n50_count']/value"/>
Bottom tier (<xsl:value-of select="$tier_two"/> bp -- end): 
  Contig number: <xsl:value-of select="//aspect[@name='supercontigs_t3_count']/value"/>
  Average length: <xsl:value-of select="//aspect[@name='supercontigs_t3_average_length']/value"/> bp
  Longest length: <xsl:value-of select="//aspect[@name='supercontigs_t3_maximum_length']/value"/> bp
  Contig bases in this tier: <xsl:value-of select="//aspect[@name='supercontigs_t3_length']/value"/> bp
  Bottom tier N50 contig length: <xsl:value-of select="//aspect[@name='supercontigs_t3_n50_length']/value"/> bp
  Bottom tier N50 contig number: <xsl:value-of select="//aspect[@name='supercontigs_t3_n50_count']/value"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:text>&#10;</xsl:text>
  <xsl:variable name="content_gc" select="//aspect[@name='content_gc']/value"/>
  <xsl:if test="$content_gc != '0'">
    <xsl:variable name="content_at" select="//aspect[@name='content_at']/value"/>
    <xsl:variable name="content_nx" select="//aspect[@name='content_nx']/value"/>
*** Genome Contents ***
Total GC count: <xsl:value-of select="$content_gc"/>, (<xsl:value-of select="format-number(($content_gc div $contigs_length), '#.0%')"/>)
Total AT count: <xsl:value-of select="$content_at"/>, (<xsl:value-of select="format-number(($content_at div $contigs_length), '#.0%')"/>)
Total NX count: <xsl:value-of select="$content_nx"/>, (<xsl:value-of select="format-number(($content_nx div $contigs_length), '#.#%')"/>)
Total: <xsl:value-of select="$contigs_length"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:text>&#10;</xsl:text>
  </xsl:if>
  <xsl:variable name="coverage_1x" select="//aspect[@name='coverage_1x']/value"/>
  <xsl:if test="$coverage_1x != 'NA'">
  <xsl:variable name="coverage_2x" select="//aspect[@name='coverage_2x']/value"/>
  <xsl:variable name="coverage_3x" select="//aspect[@name='coverage_3x']/value"/>
  <xsl:variable name="coverage_4x" select="//aspect[@name='coverage_4x']/value"/>
  <xsl:variable name="coverage_5x" select="//aspect[@name='coverage_5x']/value"/>
  <xsl:variable name="coverage_0x" select="//aspect[@name='coverage_0x']/value"/>
*** Core Gene Survey Result ***
Percentage of core genes present in this assembly: <xsl:value-of select="//aspect[@name='core_gene_present_percent']/value"/> %
Number of Core Groups present in this assembly: <xsl:value-of select="//aspect[@name='core_gene_group_present_count']/value"/>
Coregene test: <xsl:value-of select="//aspect[@name='core_gene_survey_result']/value"/>


*** Read Depth Info ***
Total consensus bases: <xsl:value-of select="$contigs_length"/>
Depth >= 5: <xsl:value-of select="$coverage_5x"/><xsl:text>&#9;</xsl:text><xsl:value-of select="format-number(($coverage_5x div $contigs_length), '#.###############')"/>
Depth >= 4: <xsl:value-of select="$coverage_4x"/><xsl:text>&#9;</xsl:text><xsl:value-of select="format-number(($coverage_4x div $contigs_length), '#.###############')"/>
Depth >= 3: <xsl:value-of select="$coverage_3x"/><xsl:text>&#9;</xsl:text><xsl:value-of select="format-number(($coverage_3x div $contigs_length), '#.###############')"/>
Depth >= 2: <xsl:value-of select="$coverage_2x"/><xsl:text>&#9;</xsl:text><xsl:value-of select="format-number(($coverage_2x div $contigs_length), '#.###############')"/>
Depth >= 1: <xsl:value-of select="$coverage_1x"/><xsl:text>&#9;</xsl:text><xsl:value-of select="format-number(($coverage_1x div $contigs_length), '#.###############')"/>
  <xsl:if test="$coverage_0x != 'NA'">
Uncovered:  <xsl:value-of select="$coverage_0x"/>
  </xsl:if>
  <xsl:text>&#10;</xsl:text>
  <xsl:text>&#10;</xsl:text>
  </xsl:if>
  <xsl:variable name="contigs_length_5k" select="//aspect[@name='contigs_length_5k']/value"/>
  <xsl:if test="$content_gc != '0'">
*** 5 Kb and Greater Contigs Info ***
Total lengths of all contigs: <xsl:value-of select="$contigs_length"/>
Total lengths of contigs 5 Kb and greater: <xsl:value-of select="$contigs_length_5k"/>
Percentage of genome: <xsl:value-of select="format-number(($contigs_length_5k div $contigs_length), '#.#%')"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:text>&#10;</xsl:text>
  </xsl:if>
  </xsl:template>
</xsl:stylesheet>
