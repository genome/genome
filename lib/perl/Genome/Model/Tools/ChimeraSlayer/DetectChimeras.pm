package Genome::Model::Tools::ChimeraSlayer::DetectChimeras;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::Model::Tools::ChimeraSlayer::DetectChimeras {
    is => 'Command::V2',
    has => [
        sequences => {
            is => 'Text',
            doc => 'sequence to run',
        },
        nastier_params => {
            is => 'Text',
            doc => 'String of params to pass to nastier',
            is_optional => 1,
        },
        chimera_slayer_params => {
            is => 'Text',
            doc => 'String of params to pass to chimera-slayer',
            is_optional => 1,
        },
        chimeras => {
            is => 'Text',
            doc => 'Link the chimeras output file (<$SEQUENCES>.out.CPS.CPC) to this file name.',
            is_mutable => 1,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    return 'Run nastier then with the output run chimera slayer .. for details on each, try gmt nastier --h and gmt chimera-slayer --h',
}

sub help_detail {
    return <<"EOS"
This command takes a sequences file and runs nastier then uses the nastier output file to run chimera slayer
EOS
}

sub execute {
    my $self = shift;
    
    unless ( -s $self->sequences ) {
        $self->error_message("Failed to find sequences file or file is zero size: ".$self->sequences );
        return;
    }

    # Build nastier/chimeraSlayer params
    my %nastier_params;
    if ( not %nastier_params = $self->build_nastier_params ) {
        $self->error_message("Failed to build nastier params");
        return;
    }
    my %chimera_slayer_params;
    if ( not %chimera_slayer_params = $self->build_chimera_slayer_params ) {
        $self->error_message("Failed to build chimera slayer params");
        return;
    }

    # Validate/create params/class
    my $nastier = Genome::Model::Tools::Nastier->create( %nastier_params );
    if ( not $nastier ) {
        $self->error_message("Failed to create Nastier class using param: ".$self->nastier_params);
        return;
    }
    my $chimera_slayer = Genome::Model::Tools::ChimeraSlayer->create( %chimera_slayer_params );
    if ( not $chimera_slayer ) {
        $self->error_message("Failed to create ChimeraSlayer class using param: ".$self->chimera_slayer_params);
        return;
    }

    # Nastier
    $self->debug_message("Run nastier...");
    if ( not $nastier->execute ) {
        $self->error_message("Failed to execute Nastier using params: ".$self->nastier_params);
        return;
    }
    $self->debug_message("Run nastier...Done");

    # Check nastier output...
    if ( not -s $nastier_params{output_file} ) {
        $self->debug_message('Nastier ran successfully, but no alignments found. Cannot run chimera slayer. Exitting.');
        return 1;
    }

    # Chimera slayer
    $self->debug_message("Run chimera slayer...");
    if ( not $chimera_slayer->execute ) {
        $self->error_message("Failed to execute ChimeraSlayer using param: ".$self->chimera_slayer_params);
        return;
    }
    $self->debug_message("Run chimera slayer...Done");

    # Check output
    my $chimeras = $self->sequences.'.out.CPS.CPC';
    if ( not -s $chimeras ) {
        $self->error_message("Failed to find chimera slayer output file or file is empty: ".$self->chimeras);
        return;
    }

    # Link chimera file to requested output
    if ( $self->chimeras ) {
        Genome::Sys->create_symlink($chimeras, $self->chimeras);
    }

    return 1;
}

sub build_nastier_params {
    my $self = shift;

    my %params;
    if ( $self->nastier_params ) {
        %params = Genome::Utility::Text::param_string_to_hash( $self->nastier_params );
        if ( $params{query_FASTA} ) {
            $self->debug_message("query_FASTA for nastier will be automatically set by this program");
            delete $params{query_FASTA};
        }
        if ( $params{output_file} ) {
            $self->debug_message("output_file for nastier will be automatically set by this program");
            delete $params{output_file};
        }
    }

    $params{query_FASTA} = $self->sequences;
    $params{output_file} = $self->sequences.'.out';
    #print Data::Dumper::Dumper \%params;

    return %params;
}

sub build_chimera_slayer_params {
    my $self = shift;

    my %params;
    if ( $self->chimera_slayer_params ) {
        %params = Genome::Utility::Text::param_string_to_hash( $self->chimera_slayer_params );
        if ( $params{exec_dir} ) {
            $self->debug_message("exec_dir for ChimeraSlayer will be automatically set by this program");
            delete $params{exec_dir};
        }
        if ( $params{query_NAST} ) {
            $self->debug_message("query_NAST for ChimeraSlayer will be automatically set by this program");
            delete $params{query_NAST};
        }
    }

    $params{query_NAST} = $self->sequences.'.out';
    $params{exec_dir} = File::Basename::dirname( $self->sequences );
    #print Data::Dumper::Dumper \%params;
    
    return %params;
}

1;

