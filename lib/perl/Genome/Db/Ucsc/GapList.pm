package Genome::Db::Ucsc::GapList;

use strict;
use warnings;

use Genome;

use Exception::Class (
    'DatabaseResolutionFailure' => {
        fields => [
            'reference_name'
        ],
    },
);

our @_HEADINGS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "ix",
    "n",
    "size",
    "type",
    "bridge",
    "bin",
);

our %_DB_TABLES_NAMES = (
    'GRCh37-lite-build37' => ['hg19', ['gap']],
    'g1k-human-build37' => ['hg19', ['gap']],
    'NCBI-human-build36' => ['hg18', [
        'chr1_gap',
        'chr2_gap',
        'chr2_gap',
        'chr3_gap',
        'chr4_gap',
        'chr5_gap',
        'chr6_gap',
        'chr7_gap',
        'chr8_gap',
        'chr9_gap',
        'chr10_gap',
        'chr11_gap',
        'chr12_gap',
        'chr13_gap',
        'chr14_gap',
        'chr15_gap',
        'chr16_gap',
        'chr17_gap',
        'chr18_gap',
        'chr19_gap',
        'chr20_gap',
        'chr21_gap',
        'chr22_gap',
        'chrX_gap',
        'chrY_gap',
    ]],
);


class Genome::Db::Ucsc::GapList {
    is => 'Command::V2',
    has_input => [
        gap_filename => {
            is => 'Text',
            is_output => 1,
        },
        reference_name => {
            is => 'Text',
        },
    ],
    has_optional_param => [
        db_user => {
            is => 'Text',
            default_value => 'genomep',
        },
        db_password => {
            is => 'Text',
            default_value => 'password',
        },
        db_host => {
            is => 'Text',
            default_value => 'genome-mysql.cse.ucsc.edu',
        },
    ],
    doc => "Fetches gap data for a reference name",
};

sub execute {
    my $self = shift;

    my $temp_output_file = eval { $self->_fetch_data_from_ucsc(); };
    if (my $e = Exception::Class->caught('DatabaseResolutionFailure')) {
        $self->status_message(sprintf(
                "Could not find gap database for reference name '%s'",
                $e->reference_name));
        return;
    }

    $self->_add_header_to_output_file();
    $self->_sort_into_output_file($temp_output_file);
    return 1;
}

sub _fetch_data_from_ucsc {
    my $self = shift;

    my ($database_name, $table_names) = $self->_resolve_database_and_table_names();
    my $query = $self->_resolve_query($table_names);

    my $temp_output_file = Genome::Sys->create_temp_file_path();
    my $command = sprintf(
        'mysql --user=%s --password=%s --host=%s -N -A -D %s -e \'%s\' > %s',
        $self->db_user, $self->db_password, $self->db_host,
        $database_name, $query, $temp_output_file);

    unless(Genome::Sys->shellcmd(cmd => $command)) {
        $self->error_message($!);
        die $self->error_message;
    }

    return $temp_output_file
}

sub _resolve_database_and_table_names {
    my $self = shift;

    my ($database_name, $table_names) = @{
        $_DB_TABLES_NAMES{$self->reference_name}};
    unless ($database_name && $table_names) {
        Exception::Class->throw('DatabaseResolutionFailure',
            reference_name => $self->reference_name);
    }

    return ($database_name, $table_names);
}

sub _resolve_query {
    my ($self, $table_names) = @_;

    my @select_statements;
    for my $tn (@{$table_names}) {
        push(@select_statements,
            sprintf("select %s from %s", join(", ", @_HEADINGS), $tn));
    }
    return join("; ", @select_statements);
}

sub _add_header_to_output_file {
    my $self = shift;

    my $header = sprintf("#%s\n", join("\t", @_HEADINGS));

    open my $fh, '>', $self->gap_filename;
    print $fh $header;
    close $fh;
}

sub _sort_into_output_file {
    my ($self, $unsorted_filename) = @_;

    unless (-s $unsorted_filename) {
        $self->error_message("Empty file passed to sort.");
        die $self->error_message;
    }

    my $command = sprintf('sort -V %s >> %s',
        $unsorted_filename, $self->gap_filename);

    Genome::Sys->shellcmd(cmd => $command);

    return $self->gap_filename;
}

1;
