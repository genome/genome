package Genome::Model::DeNovoAssembly::Command::Assemble;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Command::Assemble {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1},
    ],    
    has_constant => [
        lsf_resource => {
            is => 'Text',
            default_value => "-R 'select[type==LINUX64 && mem>30000] rusage[mem=30000] span[hosts=1]' -M 30000000",
        },
    ],
};

sub execute {
    my $self = shift;


    my $build = $self->build;
    my $processing_profile = $build->processing_profile;

    $self->debug_message('Assemble '.$build->__display_name__);

    my $assembler_class = $processing_profile->assembler_class;
    $self->debug_message('Assembler class: '. $assembler_class);

    my %assembler_params = $build->assembler_params;
    $self->debug_message('Assembler params: '.Data::Dumper::Dumper(\%assembler_params));


    my $before_assemble = $build->before_assemble;
    if ( not $before_assemble ) {
        $self->error_message('Failed to run before assemble for '.$build->__display_name__);
        return;
    }

    my $assemble = $assembler_class->create(%assembler_params);
    unless ($assemble) {
        $self->error_message("Failed to create de-novo-assemble");
        return;
    }
    $self->debug_message("Created assembler for '$assembler_class'.\n");

    eval {
        unless ($assemble->execute) {
            $self->error_message("Failed to execute de-novo-assemble execute");
            return;
        }
        $self->debug_message('Assemble...OK');
    };
    if ($@) {
        $self->error_message($@);
        die $@;
    }

    my $after_assemble = $build->after_assemble;
    if ( not $after_assemble ) {
        $self->error_message('Failed to run after assemble for '.$build->__display_name__);
        return;
    }

    return 1;
}

1;
