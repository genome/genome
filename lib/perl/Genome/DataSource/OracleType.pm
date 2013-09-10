package Genome::DataSource::OracleType;

use Genome;
use List::MoreUtils qw(any);

# Parent class for datasources connecting to our Oracle DB

class Genome::DataSource::OracleType {
    is => [
        'UR::DataSource::RDBMSRetriableOperations',
        'UR::DataSource::Oracle',
    ],
};

my @retriable_operations = (
    qr(ORA-25408), # can not safely replay call
    qr(ORA-03135), # connection lost contact
    qr(ORA-03113), # end-of-file on communication channel
);
sub should_retry_operation_after_error {
    my($self, $sql, $dbi_errstr) = @_;
    return any { $dbi_errstr =~ /$_/ } @retriable_operations;
}

sub table_and_column_names_are_upper_case { 1; }

sub _init_created_dbh {
    my ($self, $dbh) = @_;
    return unless defined $dbh;

    $self->SUPER::_init_created_dbh($dbh);

    $dbh->do('alter session set "_hash_join_enabled"=TRUE');

    # stores program name as "MODULE" and user name as "ACTION"
    $self->set_userenv(
        'dbh' => $dbh,
        'module' => substr(Cwd::abs_path($0), -48, 48)
    ); # our oracle module variable is 48 characters

    return $dbh;
}


sub _get_sequence_name_for_table_and_column {
    my ($self, $table_name, $column_name) = @_;
    if ($table_name =~ /PROCESSING_PROFILE/) {
        return 'PROCESSING_PROFILE_SEQ';
    }
    elsif($table_name =~ /GENOME_MODEL_BUILD/) {
        return 'GENOME_MODEL_EVENT_SEQ';
    }
    elsif($table_name =~ /SOFTWARE_RESULT/) {
        return 'GENOME_MODEL_EVENT_SEQ';
    }
    elsif($table_name =~ /MISC_NOTE/) {
        return 'GENOME_MODEL_SEQ';
    }
    elsif ($table_name =~ /INSTRUMENT_DATA/) {
        return 'SEQ_SEQ';
    }
    elsif ($column_name eq 'ID') {
        return $table_name . '_SEQ';
    }
    elsif($table_name =~ /GSC./) {
        return 'IMPORTED_INSTRUMENT_DATA_SEQ';
    }
    else {
        $self->SUPER::_get_sequence_name_for_table_and_column($table_name, $column_name);
    }
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
GENOME_MODEL_1
GENOME_UPDATE_HISTORY
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

1;

