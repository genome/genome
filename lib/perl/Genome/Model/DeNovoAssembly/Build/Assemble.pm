package Genome::Model::DeNovoAssembly::Build::Assemble;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Build::Assemble {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1,
        },
        sx_results => {
            is => 'Genome::InstrumentData::SxResult',
            is_many => 1,
        },
    ],
};


sub execute {
    my $self = shift;
    $self->debug_message('Assemble...');

    my $build = $self->build;
    $self->debug_message('Build: '.$build->__display_name__);

    my $assembler_class = $build->processing_profile->assembler_class;
    $self->debug_message('Assembler class: '. $assembler_class);

    my %assembler_params = $build->assembler_params;
    $self->debug_message('Assembler params: '.Data::Dumper::Dumper(\%assembler_params));

    my @sx_results = $self->sx_results;
    for my $sx_result ( @sx_results ) {
        $self->debug_message('SX result: '.$sx_result->__display_name__);
    }

    my $before_assemble = $build->before_assemble(@sx_results);
    if ( not $before_assemble ) {
        $self->error_message('Failed to run before assemble for '.$build->__display_name__);
        return;
    }

    my $assembler = $assembler_class->create(%assembler_params);
    unless ($assembler) {
        $self->error_message("Failed to create de-novo-assemble");
        return;
    }
    $self->debug_message("Created assembler for '$assembler_class'.\n");

    my $assemble_ok = eval { $assembler->execute; };
    if ( not $assemble_ok ) {
        $self->error_message($@) if $@;
        $self->error_message("Failed to execute de-novo-assemble execute");
        return;
    }

    my $after_assemble = $build->after_assemble(@sx_results);
    if ( not $after_assemble ) {
        $self->error_message('Failed to run after assemble for '.$build->__display_name__);
        return;
    }

    $self->debug_message('Assemble...OK');
    return 1;
}

1;
