package Genome::Model::Tools::Breakdancer::TigraToBreakdancer;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::Tools::Breakdancer::TigraToBreakdancer {
    is => 'Command',
    has => [
        original_breakdancer_file => {
            is => 'FilePath',
            is_input => 1,
            doc => 'Original breakdancer file, compared against tigra output file',
        },
        tigra_output_file => {
            is => 'FilePath',
            is_input => 1,
            doc => 'File produced by tigra validation, SVs in this file have passed the filter',
        },
        pass_filter_file => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Breakdancer predictions that are in the tigra output file are placed in this file',
        },
        fail_filter_file => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Breakdancer predictions that are NOT in the tigra output file are placed here',
        },
    ],
    doc => 'Produces a pass and fail file in breakdancer format using an original breakdancer file and a tigra output file',
};

sub help_detail {
    return <<EOS
Given an original breakdancer file and a tigra validation output file, this tool will put breakdancer predictions that have a matching line in the tigra file into the pass_filter_file and those that don't have a match into the fail_filter_file. These output files are in breakdancer format.
EOS
}

sub help_brief {
    return "Creates filtered breakdancer files using tigra validation output";
}

sub help_synopsis {
    return "Creates filtered breakdancer files using tigra validation output";
}

sub execute {
    my $self = shift;
    unless (-e $self->original_breakdancer_file) {
        confess 'Found no original breakdancer file at ' . $self->original_breakdancer_file;
    }
    unless (-e $self->tigra_output_file) {
        confess 'Found no tigra output file at ' . $self->tigra_output_file;
    }

    my $bd_fh = IO::File->new($self->original_breakdancer_file, 'r');
    confess 'Could not get file handle for ' . $self->original_breakdancer_file unless $bd_fh;
    my $tigra_fh = IO::File->new($self->tigra_output_file, 'r');
    confess 'Could not get file handle for ' . $self->tigra_output_file unless $tigra_fh;

    if (-e $self->pass_filter_file) {
        $self->warning_message('Removing existing pass filter file at ' . $self->pass_filter_file);
        unlink $self->pass_filter_file;
    }
    if (-e $self->fail_filter_file) {
        $self->warning_message('Removing existing fail filter file at ' . $self->fail_filter_file);
        unlink $self->fail_filter_file;
    }

    my $pass_fh = IO::File->new($self->pass_filter_file, 'w');
    confess 'Could not get file handle for ' . $self->pass_filter_file unless $pass_fh;
    my $fail_fh = IO::File->new($self->fail_filter_file, 'w');
    confess 'Could not get file handle for ' . $self->fail_filter_file unless $fail_fh;

    #make a lookup table
    my %tigra_match;
    while (my $tigra_line = $tigra_fh->getline) {
        next if $tigra_line =~ /^#/;
        my @tigra_columns = split /\s+/, $tigra_line;

        my ($tigra_pos1) = $tigra_columns[1] =~ /\d+\((\d+)\)/;
        my ($tigra_pos2) = $tigra_columns[3] =~ /\d+\((\d+)\)/;
        my ($tigra_size) = $tigra_columns[5] =~ /\d+\((\d+)\)/;
        my ($tigra_type) = $tigra_columns[6] =~ /\S+\((\S+)\)/;

        my $match_key = join '-', $tigra_columns[0], $tigra_pos1, $tigra_columns[2], $tigra_pos2, $tigra_type, abs($tigra_size);

        if (exists $tigra_match{$match_key}) {
            $self->warning_message("Line: $tigra_line seems duplicated and likely caused by assembly with diff orientation");
            #die;
        }
        else {
            $tigra_match{$match_key} = 1;
        }
    }

    while (my $bd_line = $bd_fh->getline) {
        next if $bd_line =~ /^#/;
        my @bd_columns = split /\s+/, $bd_line;

        my $match_key = join '-', $bd_columns[0], $bd_columns[1], $bd_columns[3], $bd_columns[4], $bd_columns[6], abs($bd_columns[7]);

        if (exists $tigra_match{$match_key}) {
            $pass_fh->print($bd_line);
            delete $tigra_match{$match_key};
        }
        else {
            $fail_fh->print($bd_line);
        }
    }

    if (%tigra_match) {
        my $msg = join "\n", keys %tigra_match;
        $self->warning_message("Some Tigra output can not match to breakdancer output:\n$msg");
        #die;
    }

    $bd_fh->close;
    $tigra_fh->close;
    $pass_fh->close;
    $fail_fh->close;
    return 1;
}

1;

