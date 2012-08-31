package Genome::Model::Tools::DetectVariants::Somatic;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants::Somatic {
    is => 'Genome::Model::Tools::DetectVariants',
    has => [
        control_aligned_reads_input => {
            is => 'Text',
            doc => 'Location of the control aligned reads file to which the input aligned reads file should be compared',
            shell_args_position => '2',
            is_input => 1,
            is_output => 1,
        },
    ]
};

sub help_brief {
    "A selection of somatic variant detectors.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants somatic ...
EOS
}

sub help_detail {
    return <<EOS 
Tools to run somatic variant detectors with a common API and output their results in a standard format.
EOS
}

sub _verify_inputs {
    my $self = shift;
    
    my $control_aligned_reads_file = $self->control_aligned_reads_input;
    unless(Genome::Sys->check_for_path_existence($control_aligned_reads_file)) {
        $self->error_message("control aligned reads input $control_aligned_reads_file was not found.");
        return;
    }
    
    return $self->SUPER::_verify_inputs;
}

1;
