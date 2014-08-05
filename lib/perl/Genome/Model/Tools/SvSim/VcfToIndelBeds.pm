package Genome::Model::Tools::SvSim::VcfToIndelBeds;

use Genome;
use strict;
use warnings;

my $JOINX_PATH = "/usr/bin/joinx1.9";

class Genome::Model::Tools::SvSim::VcfToIndelBeds {
    is => "Command::V2",
    has_input => [
        input_vcf => {
            is => "Text",
            doc => "Input vcf to process",
        },

        fasta_path => {
            is => "Text",
            doc => "The fasta file for the input vcf",
        },

        insertion_output_bed => {
            is => "Text",
            doc => "Output bed file of insertions",
        },

        deletion_output_bed => {
            is => "Text",
            doc => "Output bed file of deletions",
        },
    ],
    has_transient_optional => [
        _insertion_fh => {},
        _deletion_fh => {},
    ]
};

sub _vcf_to_raw {
    my $self = shift;
    my $temp_path = Genome::Sys->create_temp_file_path;
    my @cmd = (
        $JOINX_PATH, "vcf2raw",
        "-f", $self->fasta_path,
        "-o", $temp_path,
        $self->input_vcf
        );

    my $cmd = join(" ", @cmd);
    Genome::Sys->shellcmd(cmd => $cmd);

    return $temp_path;
}

sub _output_insertion {
    my ($self, $chr, $pos, $ref, $alt) = @_;
    # Convert to space-based coordinates
    $pos -= 1;

    my $svlen = length($alt);
    $self->_insertion_fh->printf("%s\t%d\t%d\t%d\n", $chr, $pos, $pos, $svlen);
}


sub _output_deletion {
    my ($self, $chr, $pos, $ref, $alt) = @_;
    # Convert to space-based coordinates
    $pos -= 1;
    my $svlen = length($ref);
    my $stop = $pos + $svlen;

    $self->_deletion_fh->printf("%s\t%d\t%d\t%d\n", $chr, $pos, $stop, $svlen);
}

sub execute {
    my $self = shift;
    my $raw_path = $self->_vcf_to_raw($self->input_vcf);
    my $fh = Genome::Sys->open_file_for_reading($raw_path);

    my $ins_path = $self->insertion_output_bed;
    my $del_path = $self->deletion_output_bed;
    $self->_insertion_fh(Genome::Sys->open_file_for_writing($ins_path));
    $self->_deletion_fh(Genome::Sys->open_file_for_writing($del_path));

    while (my $line = $fh->getline) {
        chomp $line;
        my ($chr, $pos, $call) = split(/\t/, $line);
        my ($ref, $alt) = split("/", $call);
        if ($ref eq '*') {
            $self->_output_insertion($chr, $pos, $ref, $alt);
        }
        elsif ($alt eq '*') {
            $self->_output_deletion($chr, $pos, $ref, $alt);
        }
    }
    $self->_insertion_fh->close();
    $self->_deletion_fh->close();

    return 1;
}

1;
