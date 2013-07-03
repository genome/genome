#extract reads from a set of bams in breakdancer predicted regions
package Genome::Model::Tools::DetectVariants2::Filter::TigraValidationFull;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::TigraValidationFull {
    is  => 'Genome::Model::Tools::DetectVariants2::Filter::TigraValidation',
    has_optional => [
        specify_chr => {
            is => 'String',
            is_input => 1,
            default => 'all',
        },
    ],
};

#this should go in some sort of BAM class

sub _chromosome_list_from_fasta {
    my $self = shift;
    my $fasta = shift;
    my @chr_list;
    my $fh = Genome::Sys->open_file_for_reading($fasta . ".fai");
    while(my $line = $fh->getline) {
        my ($chr) = split /\t/, $line;
        push @chr_list, $chr;
    }
    return @chr_list;
}

sub _full_chromosome_list {
     my ($self) = @_;
     return $self->_chromosome_list_from_fasta($self->reference_sequence_input);
}

sub _get_split_object {
    my ($self) = @_;
    return Genome::Model::Tools::Breakdancer::SplitFiles->create(
            input_file           => $self->_breakdancer_input,
            output_directory     => $self->_temp_staging_directory,
            output_file_template => 'svs.hq.tigra.CHR',
            create_other         => 0,
        );
}

1;
