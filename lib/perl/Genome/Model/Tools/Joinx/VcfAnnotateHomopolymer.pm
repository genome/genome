package Genome::Model::Tools::Joinx::VcfAnnotateHomopolymer;

use strict;
use warnings;

use Genome;

our $MINIMUM_JOINX_VERSION = 1.9; 

class Genome::Model::Tools::Joinx::VcfAnnotateHomopolymer {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file => {
            is  => 'Text',
            doc => 'Vcf File to filter',
            shell_args_position => 1,
        },
        homopolymer_file => {
            is  => 'Text',
            doc => 'Bed File containing homopolymer',
            shell_args_position => 2,
        },
    ],
    has_optional_input => [
        output_file => {
            is  => 'Text',
            doc => 'The output file (defaults to stdout)',
            is_output => 1,
        },
        max_length => {
            is  => 'Integer',
            doc => 'maximum indel length to annotate as in homopolymer',
            default_value => 2,
        },
        info_field => {
            is  => 'Text',
            doc => 'Info field id for homopolymer',
            default => 'HOMP_FILTER',
        },
        use_bgzip => {
            is  => 'Boolean',
            doc => 'zcats the input files into stdin, and bgzips the output',
            default => 0,
        },
    ],
};

sub help_brief {
    "Annotate information from one homopolymer bed file"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    This tool wraps the joinx vcf-annotate-homopolymers command which places homopolymer information from one bed file into vcf INFO field.
EOS
}


sub execute {
    my $self = shift;

    my $max_length       = $self->max_length;
    my $info_field       = $self->info_field;
    my $output_file      = $self->output_file;
    my $input_file       = $self->input_file;
    my $use_bgzip        = $self->use_bgzip;
    my $homopolymer_file = $self->homopolymer_file;

    if ($self->use_version < $MINIMUM_JOINX_VERSION) {
        die $self->error_message("This module requires joinx version $MINIMUM_JOINX_VERSION or higher to function correctly.");
    }

    if ($use_bgzip && ! $output_file){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth."); 
    }

    my $output = "-";
    $output = $output_file if $output_file;

    unless(-s $input_file) {
        die $self->error_message("$input_file does not exist");
    }

    unless (-s $homopolymer_file) {
        die $self->error_message("Homopolymer file: $homopolymer_file does not exist");
    }

    if ($use_bgzip){
        $input_file = "<(zcat $input_file)";
    }

    my $cmd = $self->joinx_path . " vcf-annotate-homopolymers --vcf-file $input_file --bed-file $homopolymer_file --max-length $max_length --info-field-name $info_field";
    
    if ($output_file && not defined $use_bgzip) {
        $cmd .= " --output-file $output" if $output_file;
    } 
    elsif (defined $use_bgzip && $output_file) {
        $cmd .= " | bgzip -c > $output";
        $cmd = "bash -c \"$cmd\"";
    }
        
    my %params = (cmd => $cmd);
    $params{output_files} = [$output] if $output ne "-";
    $params{skip_if_output_is_present} = 0;

    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
