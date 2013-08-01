package Genome::InstrumentData::AlignedBamResult;

use Genome;

use warnings;
use strict;

class Genome::InstrumentData::AlignedBamResult {
    is => 'Genome::SoftwareResult::Stageable',
    is_abstract => 1,
    attributes_have => [
        is_output => { is => 'Boolean', is_optional => 1, },
    ],
    has_input => [
        reference_build => { # PROVIDES fasta VIA full_consensus_path('fa')
            is => 'Genome::Model::Build::ImportedReferenceSequence',
        },
    ],
    has_constant => [
        # from inputs
        reference_fasta => { 
            calculate_from => [qw/ reference_build /],
            calculate => q| return $reference_build->full_consensus_path('fa'); |, 
        },
        # from output
        bam_path => { # this is called bam_file in merged
            is_output => 1,
            calculate => q| return $self->output_dir.'/'.$self->id.'.bam'; |, 
        },
        bam_file => { # alias
          is => 'Text',
          via => '__self__',
          to => 'bam_path',
        },
        # flagstat
        bam_flagstat_path => {
            calculate_from => [qw/ bam_path /],
            calculate => q| return $bam_path.'.flagstat'; |,
        },
        bam_flagstat_file => { # alias
          is => 'Text',
          via => '__self__',
          to => 'bam_flagstat_path',
        },
        # md5
        bam_md5_path => {
            calculate_from => [qw/ bam_path /],
            calculate => q| return $bam_path.'.md5'; |,
        },
    ],
};

sub run_flagstat_on_output_bam_path {
    my $self = shift;
    $self->status_message('Run flagstat on output bam file...');

    my $bam_path = $self->bam_path;
    if ( not $bam_path or not -s $bam_path ) {
        $self->error_message('Bam file not set or does not exist!');
        return;
    }

    my $flagstat_path = $self->bam_flagstat_path;
    $self->status_message("Flagstat file: $flagstat_path");
    my $cmd = "samtools flagstat $bam_path > $flagstat_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    $self->status_message('Flagstat output:');
    $self->status_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );
    if ( not $flagstat->{total_reads} > 0 ) {
        $self->error_message('Flagstat determined that there are no reads in bam! '.$bam_path);
        return;
    }

    $self->status_message('Run flagstat on output bam file...done');
    return $flagstat;
}

sub run_md5sum_on_output_bam_path {
    my $self = shift;
    $self->status_message('Run MD5SUM on output bam file...');

    my $md5sum = $self->load_md5sum;
    if ( $md5sum ) {
        $self->status_message("MD5SUM: $md5sum");
        return $md5sum;
    }

    my $bam_path = $self->bam_path;
    if ( not $bam_path or not -s $bam_path ) {
        $self->error_message('Bam file does not exist!');
        return;
    }
    $self->status_message('Bam path: '.$bam_path);

    my $md5_path = $self->bam_md5_path;
    $self->status_message('MD5SUM path: '.$md5_path);


    my $cmd = "md5sum $bam_path > $md5_path";
    my $rv = eval{ 
        Genome::Sys->shellcmd(
            cmd                        => $cmd, 
            input_files                => [$bam_path],
            output_files               => [$md5_path],
            skip_if_output_is_present  => 0,
        ); 
    };
    if ( not $rv or not -s $md5_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run md5sum!');
        return;
    }

    $md5sum = $self->load_md5sum;
    return if not $md5sum;
    $self->status_message("MD5SUM: $md5sum");

    $self->status_message('Run MD5SUM on output bam file...done');
    return $md5sum;
}

sub load_md5sum {
    my $self = shift;

    my $md5_path = $self->bam_md5_path;
    return if not -s $md5_path;

    my $md5_fh = eval{ Genome::Sys->open_file_for_reading($md5_path); };
    if ( not $md5_fh ) {
        $self->error_message('Failed to open md5 path! '.$md5_path);
        return;
    }
    my ($md5sum) = split(/\s+/, $md5_fh->getline);
    $md5_fh->close;

    return $md5sum;
}

1;

