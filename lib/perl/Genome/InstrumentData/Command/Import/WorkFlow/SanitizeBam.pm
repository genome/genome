package Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam { 
    is => 'Command::V2',
    has_input => [
        dirty_bam_path => {
            is => 'Text',
            doc => 'The path of the dirty bam to clean.',
        }
    ],
    has_output => [ 
        clean_bam_path => {
            is => 'Text',
            calculate_from => [qw/ dirty_bam_path /],
            calculate => q{
                $dirty_bam_path =~ s/(\.bam)$/.clean$1/;
                $dirty_bam_path; 
                return $dirty_bam_path;
            },
            doc => 'The path of the clean bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Sanitize bam...');

    my $sort_ok = $self->_sort_bam;
    return if not $sort_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    my $cleanup_ok = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->remove_paths_and_auxiliary_files($self->dirty_bam_path);
    return if not $cleanup_ok;

    $self->debug_message('Sanitize bam...done');
    return 1;
}

sub _sort_bam {
    my $self = shift;

    my $dirty_bam_path = $self->dirty_bam_path;
    $self->debug_message("Dirty bam path: $dirty_bam_path");

    my $clean_bam_path = $self->clean_bam_path;
    $self->debug_message("Clean bam path: $clean_bam_path");

    my $cmd = "/usr/bin/seq-grind sanitize --input $dirty_bam_path --output $clean_bam_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $clean_bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run samtools sort!');
        return;
    }

    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $dirty_flagstat = $helpers->load_or_run_flagstat($self->dirty_bam_path);
    return if not $dirty_flagstat;

    my $clean_flagstat = $helpers->load_or_run_flagstat($self->clean_bam_path);
    return if not $clean_flagstat;

    $self->debug_message('Clean bam read count: '.$clean_flagstat->{total_reads});
    $self->debug_message('Dirty bam read count: '.$dirty_flagstat->{total_reads});

    if ( $clean_flagstat->{total_reads} != $dirty_flagstat->{total_reads} ) {
        $self->error_message('Clean and dirty bam read counts do not match!');
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

1;

