#eclark: Ignore table list should be reviewed.  I've seen classes for at few of them.

use strict;
use warnings;

package Genome::DataSource::GMSchemaOracle;

use Genome;

use Cwd;

class Genome::DataSource::GMSchemaOracle {
    is => ['UR::DataSource::Oracle'],
    type_name => 'genome datasource gmschema',
};

sub server {
    "dwrac";
}

sub login {
    "mguser";}

sub auth {
    "mguser_prd";}

sub owner {
    "MG";}


sub table_and_column_names_are_upper_case { 1; }

sub _sync_database {
    my $self = shift;

    $DB::single = 1;

    my $dbh = $self->get_default_handle;
    unless ($dbh->do("alter session set NLS_DATE_FORMAT = 'YYYY-MM-DD HH24:MI:SS'")
            and
            $dbh->do("alter session set NLS_TIMESTAMP_FORMAT = 'YYYY-MM-DD HH24:MI:SSXFF'"))
    {
        Carp::croak("Can't set date format: $DBI::errstr");
    }
    $self->SUPER::_sync_database(@_);
}

# A list of the old GM schema tables that Genome::Model should ignore
my @OLD_GM_TABLES = qw(
ALL_ALLELE_TYPE
ANALYSIS_METHOD
CHROMOSOME
COLLABORATOR
COLLABORATOR_SAMPLE
CONSERVATION_SCORE
EGI_TYPE
EXTERNAL_GENE_ID
GENE
GENE_EXPRESSION
GENE_GENE_EXPRESSION
GENE_GENOMIC_FEATURE
GENE_GENOTYPE
GENOMIC_FEATURE
GENOTYPE
GENOTYPE_VARIATION
GE_DETECTION
GE_TECH_TYPE
GF_FEATURE_TYPE
GO_XREF
GROUP_INFO
GROUP_TYPES
HISTOLOGY
IPRO_GENE_TRANSCRIPT_XREF
IPRO_RESULTS
MAF
META_GROUP
META_GROUP_STATUS
PF_FEATURE_TYPE
PP_DETECTION_SOFT
PP_MAPPING_REFERENCE
PP_TECH_TYPE
PROCESS_PROFILE
PROCESS_PROFILE_SOURCE
PROTEIN
PROTEIN_FEATURE
PROTEIN_VARIATION
PROTEIN_VARIATION_SCORE
PVS_SCORE_TYPE
PVS_SOFTWARE
READ_GROUP
READ_GROUP_GENOTYPE
READ_GROUP_GENOTYPE_OLD
READ_GROUP_INFO
READ_INFO
REPEAT_INFO
RGGI_INFO_TYPE
RGG_INFO
RGG_INFO_OLD
RI_REPEAT_TYPE
SAMPLE
SAMPLE_GENE
SAMPLE_GENE_EXPRESSION
SAMPLE_GENOTYPE
SAMPLE_GROUP_INFO
SAMPLE_HISTOLOGY
SAMPLE_SUBTYPE
SAMPLE_TISSUE
SAMPLE_TYPE
SEQUENCE_DIFF
SEQUENCE_DIFF_EVAL
SEQUENCE_DIFF_PART
SGIAM_SCORE_TYPE
SGI_ANALYSIS_METHOD
SG_INFO_TYPE
STRAND
SUBMITTER
SUBMITTER_METHOD
TERM
TISSUE
TRANSCRIPT
TRANSCRIPT_SOURCE
TRANSCRIPT_STATUS
TRANSCRIPT_SUB_STRUCTURE
TRANSCRIPT_VARIATION
TSS_STRUCTURE_TYPE
TV_OLD
TV_TYPE
VARIANT_REVIEW_DETAIL_OLD
VARIANT_REVIEW_LIST_FILTER
VARIANT_REVIEW_LIST_MEMBER
VARIATION
VARIATION_FREQUENCY
VARIATION_GROUP
VARIATION_INSTANCE
VARIATION_ORIG
VARIATION_ORIG_VARIATION_TYPE
VARIATION_SCORE
VS_SCORE_TYPE
VS_SOFTWARE
);

sub _ignore_table {
    my($self,$table_name) = @_;

    return 1 if $self->SUPER::_ignore_table($table_name);

    return scalar(grep { $_ eq $table_name } @OLD_GM_TABLES);
}

sub _resolve_class_name_for_table_name_fixups {
    my($self,@words) = @_;

    if ($words[0] eq 'Genome') {
        # Everything is already under Genome
        shift @words;
        if ($words[0] eq 'Model') {
            splice(@words, 1, 0, '::');  # Make it Model::Blah instead of ModelBlah
        }
    }
    return @words;
}


1;

