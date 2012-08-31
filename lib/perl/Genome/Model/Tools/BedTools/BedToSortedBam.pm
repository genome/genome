package Genome::Model::Tools::BedTools::BedToSortedBam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BedTools::BedToSortedBam {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input BED file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output BAM file',
            is_optional => 1,
        },
        genome => {
            is => 'Text',
            doc => 'The genome reference.  See BEDTools User Manual for documentation of format.',
            is_optional => 1,
            default_value => Genome::Model::Tools::BedTools::BedToBam->default_genome,
        },
        bed12 => {
            is => 'Boolean',
            doc => 'The BED file is in BED12 format.  The BAM CIGAR string will reflect BED blocks.',
            default_value => Genome::Model::Tools::BedTools::BedToBam->default_bed12,
            is_optional => 1,
        },
        mapq => {
            is => 'Integer',
            doc => 'Set the mappinq quality for the BAM records. (INT) Default: 255',
            is_optional => 1,
        },
        samtools_version => {
            is => 'Text',
            is_optional => 1,
            default_value => Genome::Model::Tools::Sam->default_samtools_version,
        },
    ],
};

sub help_brief {
    'Converts BED/GFF/VCF features to sorted BAM format file.',
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools bed-to-sorted-bam --input-file=example.bed --output-file=example.bam
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;

    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->input_file,qw/\.bed/);
    unless ($basename && $suffix) {
        die('Failed to parse BED file path:  '. $self->input_file);
    }
    my $unsorted_bam = Genome::Sys->create_temp_file_path($basename .'_unsorted.bam');
    unless ($self->output_file) {
        $self->output_file($dirname .'/'. $basename .'.bam');
    }
    my %params = (
        input_file => $self->input_file,
        output_file => $unsorted_bam,
        use_version => $self->use_version,
        genome => $self->genome,
    );
    if (defined($self->bed12)) {
        $params{bed12} = $self->bed12;
    }
    if (defined($self->mapq)) {
        $params{mapq} = $self->mapq;
    }
    my $BedToBam = Genome::Model::Tools::BedTools::BedToBam->execute(%params);
    unless ($BedToBam) {
        die('Failed to run BedToBam command with params:  '. Data::Dumper::Dumper(%params));
    }
    my $SortBam = Genome::Model::Tools::Sam::SortBam->execute(
        file_name => $unsorted_bam,
        output_file => $self->output_file,
        use_version => $self->samtools_version,
    );
    unless ($SortBam) {
        die('Failed to sort BAM file:  '. $unsorted_bam);
    }
    my $IndexBam = Genome::Model::Tools::Sam::IndexBam->execute(
        bam_file => $self->output_file,
        use_version => $self->samtools_version,
    );
    unless ($IndexBam) {
        die('Failed to index BAM file:  '. $self->output_file);
    }
    return 1;
}

1;
