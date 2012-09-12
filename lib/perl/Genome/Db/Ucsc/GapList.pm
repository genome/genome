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
    'NCBI-human-build36' => ['hg18', 'gap'],
    'GRCh37-lite-build37' => ['hg19', 'gap'],
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
    return $self->gap_filename;
}

sub _fetch_data_from_ucsc {
    my $self = shift;

    my ($database_name, $table_name) = $self->_resolve_database_and_table_names();
    my $query = $self->_resolve_query($table_name);

    my $temp_output_file = Genome::Sys->create_temp_file_path();
    my $command = sprintf(
        "mysql --user=%s --password=%s --host=%s -N -A -D %s -e '%s' > %s",
        $self->db_user, $self->db_password, $self->db_host,
        $database_name, $query, $temp_output_file);

    Genome::Sys->shellcmd(cmd => $command);

    return $temp_output_file
}

sub _resolve_database_and_table_names {
    my $self = shift;

    my ($database_name, $table_name) = @{
        $_DB_TABLES_NAMES{$self->reference_name}};
    unless ($database_name && $table_name) {
        Exception::Class->throw('DatabaseResolutionFailure',
            reference_name => $self->reference_name);
    }

    return ($database_name, $table_name);
}

sub _resolve_query {
    my ($self, $table_name) = @_;

    return sprintf("select %s from %s", join(", ", @_HEADINGS), $table_name);
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

    my $command = sprintf("sort -V > %s",
        $unsorted_filename, $self->gap_filename);

    Genome::Sys->shellcmd(cmd => $command);

    return $self->gap_filename;
}

1;
