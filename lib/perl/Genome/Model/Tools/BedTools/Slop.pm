package Genome::Model::Tools::BedTools::Slop;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BedTools::Slop {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The BED file to add slop',
        },
        genome_file => {
            is => 'Text',
            doc => 'The genome file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output BED file with slop',
            is_output => 1,
        },
        both => {
            is => 'Integer',
            doc => 'Increase the BED/GFF/VCF entry by base pairs in each direction. - (Integer) or (Float, e.g. 0.1) if used with -pct.',
            is_optional => 1,
        },
        left => {
            is => 'Integer',
            doc => 'The number of base pairs to subtract from the start coordinate. - (Integer) or (Float, e.g. 0.1) if used with -pct.',
            is_optional => 1,
        },
        right => {
            is => 'Integer',
            doc => 'The number of base pairs to add to the end coordinate. - (Integer) or (Float, e.g. 0.1) if used with -pct.',
            is_optional => 1,
        },
        strand => {
            is => 'Boolean',
            doc => 'Define -l and -r based on strand.  E.g. if used, -l 500 for a negative-stranded feature, it will add 500 bp downstream.  Default = false.',
            default_value => 0,
            is_optional => 1,
        },
        percent => {
            is => 'Boolean',
            doc => 'Define -l and -r as a fraction of the feature\'s length.  E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp "upstream".  Default = false.',
            default_value => 0,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    "Add requested base pairs of \"slop\" to each feature",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools slop --input-file example.bed --output-file sloppy.bed --genome-file all_sequences.genome --base-pair=250
EOS
}

sub help_detail {                           
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;

    my @input_files = ($self->input_file, $self->genome_file);
    my $cmd = $self->bedtools_path .'/bin/slopBed -i '. $self->input_file .' -g '. $self->genome_file;
    if ($self->both) {
        if ($self->left || $self->right) {
            die('Do not define left or right with the both option!');
        }
        $cmd .= ' -b '. $self->both;
    } else {
        if ($self->left) {
            $cmd .= ' -l '. $self->left;
        }
        if ($self->right) {
            $cmd .= ' -r '. $self->right;
        }
    }
    if ($self->strand) {
        $cmd .= ' -s';
    }
    if ($self->percent) {
        $cmd .= ' -pct';
    }
    my @output_files = ($self->output_file);
    $cmd .= ' > '. $self->output_file;

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );

    return 1;
}

1;
