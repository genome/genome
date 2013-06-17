package Genome::Model::Tools::Joinx::VcfNormalizeIndels;

use strict;
use warnings;

use Genome;

our $MINIMUM_JOINX_VERSION = 1.7;

class Genome::Model::Tools::Joinx::VcfNormalizeIndels {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Vcf file to normalize',
            shell_args_position => 1,
        },
        reference => {
            is => 'Text',
            doc => 'Reference sequence fasta file',
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The output file (defaults to stdout)',
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'zcats the input file into stdin, and bgzips the output',
            default => 0,
        },
    ],
};

sub help_brief {
    'Normalizes indels in a vcf file.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt joinx vcf-normalize-indels a.vcf --reference b.fasta --output-file normalized.vcf
EOS
}

sub execute {
    my ($self) = @_;
    unless ($self->input_file) {
        die $self->error_message("Your must provide an input_file");
    }

    $self->_validate_inputs;
    my %params = $self->_generate_joinx_command;
    Genome::Sys->shellcmd(%params);
    return 1;
}

sub _validate_inputs {
    my $self = shift;
    
    if($self->use_version < $MINIMUM_JOINX_VERSION) {
        die $self->error_message("This module requires joinx version $MINIMUM_JOINX_VERSION or higher to function correctly.");
    }

    my $input_file = $self->input_file;
    unless(-s $input_file) {
        die $self->error_message("$input_file does not exist");
    }

    my $reference = $self->reference;
    unless(-s $reference) {
        die $self->error_message("$reference does not exist");
    }

    if($self->use_bgzip && not defined($self->output_file)){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth."); 
    }

    return 1;
}

sub _generate_joinx_command {
    my $self = shift;
    my $output = defined $self->output_file ? $self->output_file : '-';
    my $input = $self->use_bgzip ? join('', '<(zcat ', $self->input_file, ')') 
                                 : $self->input_file; 
    my $reference = $self->reference;
    my $cmd = $self->joinx_path . " vcf-normalize-indels" . " --input-file $input" . " --fasta $reference";
    if ($self->output_file && not($self->use_bgzip)){
        $cmd .= " --output-file $output";
    } elsif( $self->use_bgzip && $self->output_file){
        $cmd .= " | bgzip -c > $output";
    }
    $self->status_message($cmd);

    my %params = (
        cmd => $cmd,
    );
    $params{output_files} = [$output] if $output ne '-';
    return %params;
}

1;
