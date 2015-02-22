package Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam { 
    is => 'Command::V2',
    has_input => {
        bam_path => {
            is => 'Text',
            doc => 'The path of the dirty bam to clean.',
        }
    },
    has_output => {
        output_bam_path => {
            is => 'Text',
            calculate_from => [qw/ bam_path /],
            calculate => q{
                $bam_path =~ s/(\.bam)$/.clean$1/;
                return $bam_path;
            },
            doc => 'The path of the clean bam.',
        },
    },
    has_transient => {
        _sanitize_metrics => { is => 'Hash', },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('Sanitize bam...');

    my $sanitize_ok = $self->_sanitize_bam;
    return if not $sanitize_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    my $cleanup_ok = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->remove_paths_and_auxiliary_files($self->bam_path);
    return if not $cleanup_ok;

    $self->debug_message('Sanitize bam...done');
    return 1;
}

sub _sanitize_bam {
    my $self = shift;

    my $bam_path = $self->bam_path;
    $self->debug_message("Dirty bam path: $bam_path");

    my $output_bam_path = $self->output_bam_path;
    $self->debug_message("Clean bam path: $output_bam_path");

    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    my $sanitize_err = $tmp_dir.'/sanitize.err';
    my $cmd = "/usr/bin/seq-grind sanitize --input $bam_path --output $output_bam_path 2> $sanitize_err";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $output_bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run samtools sanitize!');
        return;
    }

    my $load_metrics = $self->_load_sanitize_metrics($sanitize_err);
    return if not $load_metrics;

    return 1;
}

sub _load_sanitize_metrics {
    my ($self, $sanitize_err) = @_;

    my $fh = eval{ Genome::Sys->open_file_for_reading($sanitize_err); };
    if ( not $fh ) {
        $self->error_message('Failed to open sanitize error file to retrieve metrics.');
        return;
    }

    my %metrics;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my ($name, $value) = $line =~ /^seq-grind sanitize sequences (\w+): (\d+)$/;
        next if not $name;
        $metrics{$name} = $value;
    }

    $fh->close;
    $self->_sanitize_metrics(\%metrics);

    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $dirty_flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $dirty_flagstat;

    my $clean_flagstat = $helpers->load_or_run_flagstat($self->output_bam_path);
    return if not $clean_flagstat;

    $self->debug_message('Clean bam read count: '.$clean_flagstat->{total_reads});
    $self->debug_message('Dirty bam read count: '.$dirty_flagstat->{total_reads});
    my $sanitize_metrics = $self->_sanitize_metrics;
    $self->debug_message('Sanitize filtered:   '.$sanitize_metrics->{filtered});

    if ( $sanitize_metrics->{input} != $dirty_flagstat->{total_reads} ) {
        $self->error_message('Reads input into sanitize does not match reads in dirty bam: %s vs. %s', $sanitize_metrics->{input}, $dirty_flagstat->{total_reads});
        return;
    }

    if ( $self->_sanitize_metrics->{output} != $clean_flagstat->{total_reads} ) {
        $self->error_message('Reads output by sanitize does not match reads in clean bam: %s vs. %s', $sanitize_metrics->{output}, $clean_flagstat->{total_reads});
        return;
    }

    my $filtered = $dirty_flagstat->{total_reads} - $clean_flagstat->{total_reads};
    if ( $self->_sanitize_metrics->{filtered} != $filtered ) {
        $self->error_message('Reads filtered by sanitize does not match the difference of dirty - clean reads: %s vs. %s', $sanitize_metrics->{filtered}, $filtered);
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

1;

