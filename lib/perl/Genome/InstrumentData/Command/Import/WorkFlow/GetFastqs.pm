package Genome::InstrumentData::Command::Import::WorkFlow::GetFastqs;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::GetFastqs { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for fastqs.',
        },
        source_fastq_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'Fastqs to transfer.',
        },
    ],
    has_output => [
        fastq_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'Retrieved fastq paths.',
        }, 
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Get fastqs...');

    my @fastq_paths;
    for my $source_fastq_path ( $self->source_fastq_paths ) {
        my $fastq_path = $self->_transfer_fastq($source_fastq_path);
        return if not $fastq_path;
        push @fastq_paths, $fastq_path;
    }
    $self->fastq_paths(\@fastq_paths);

    my $read_counts_ok = $self->_verify_read_counts;
    return if not $read_counts_ok;

    $self->status_message('Get fastqs...done');
    return 1;
}

sub _transfer_fastq {
    my ($self, $fastq) = @_;

    my $dir = $self->working_directory;
    $self->status_message("Fastq: $fastq");
    my $fastq_base_name = File::Basename::basename($fastq);
    $fastq_base_name =~ s/\.gz$//;
    my $destination_file = $dir.'/'.$fastq_base_name;
    $self->status_message("Destination file: $destination_file");

    my $cmd;
    if ( $fastq =~ /^http/ ) {
        my $log = $dir.'/'.$fastq_base_name.'.log';
        $cmd = "wget -o $log -O - $fastq | zcat - | tee $destination_file";
    }
    elsif ( $fastq =~ /\.gz$/ ) { # zcat
        $cmd .= "zcat $fastq | tee $destination_file";
    }
    else {
        $cmd .= "tee $destination_file < $fastq";
    }

    my $line_count_file = $destination_file.'.count';
    $cmd .= " | wc -l > $line_count_file";
    $self->status_message("Line count file: $line_count_file");

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $destination_file ) {
        $self->error_message('Failed to transfer source file to tmp directory!');
        return;
    }

    return $destination_file;
}

sub _verify_read_counts {
    my $self = shift;
    $self->status_message('Verify read counts...');

    my @fastq_paths = $self->fastq_paths;
    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $fastq_paths_and_read_counts = $helpers->load_read_count_for_fastq_paths(@fastq_paths);
    return if not $fastq_paths_and_read_counts;

    $self->status_message('Source file read counts: '.Dumper($fastq_paths_and_read_counts));
    my @read_counts = List::MoreUtils::uniq( values %$fastq_paths_and_read_counts );
    if ( @read_counts > 1 ) {
        $self->error_message('Read counts are not the same for source fastqs!');
        return;
    }

    $self->status_message('Verify read counts...done');
    return 1;
}


1;

