package Genome::Model::DeNovoAssembly::Command::PostAssemble;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Command::PostAssemble {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1},
    ],
    has_constant => [

        lsf_resource => { is => 'Text',
            default_value => "-R 'select[type==LINUX64 && mem>30000] rusage[mem=30000] span[hosts=1]' -M 30000000",
        },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $processing_profile = $self->build->processing_profile;
    unless ( $processing_profile->post_assemble ) { #shouldn't happen
        $self->error_message(
            "Can't run de-novo-assembly post-assemble without declearing it in processing profile");
        return;
    }

    for my $post_assemble_part ( $processing_profile->post_assemble_parts ) {
        $self->_execute_tool ( $post_assemble_part );
    }

    return 1;
}

sub _execute_tool {
    my ( $self, $post_assemble_part ) = @_;

    $self->debug_message('Post assemble: '.$post_assemble_part);

    #get tool name, convert to class name
    my ( $tool_name ) = $post_assemble_part =~ /^(\S+)/;
    $tool_name =~ s/-/ /g;
    my $class_name = Genome::Utility::Text::string_to_camel_case( $tool_name );
    my $class = 'Genome::Model::Tools::'
        . ucfirst $self->build->processing_profile->assembler_base_name
            . '::' . $class_name;

    $self->debug_message('Class: '.$class);

    #get params string, convert to hash, append assembly directory
    my ( $params_string ) = $post_assemble_part =~ /^\S+\s+(.*)/;
    my %params;
    if ( $params_string ) {
        %params = Genome::Utility::Text::param_string_to_hash( $params_string );
    }

    for my $param_name ( keys %params ) {
        my $value = $params{$param_name};
        if ( $value =~ /:/ ) {
            my ( $class, $method ) = split ( ':', $value );
            # class should just be build class for now .. enforced at pp creation
            if ( not $class eq 'build' ) {
                $self->debug_message("Skipping running tool, $class_name.  Can not derive tool param value from class, $class.  Expected class to be genome::model::Build but got $class");
                next;
            }
            my $derived_value;
            if ( not $derived_value = $self->$class->$method ) {
                $self->debug_message("Skipping running tool, $class_name.  Failed to get param_value from class, $class, using method, $method");
                next;
            }
            $derived_value = lc $derived_value if not ref $derived_value; # lc value if values are not HASH OR ARRAY REFS
            $params{$param_name} = $derived_value
        }
    }

    $params{assembly_directory} = $self->build->data_directory;

    $self->debug_message('Params: '.Data::Dumper::Dumper(\%params));

    my $tool = $class->create(%params);
    unless ( $tool ) {
        $self->error_message("Failed to create post assemble process: $post_assemble_part");
        return;
    }

    $tool->dump_status_messages(1);

    unless ( $tool->execute ) {
        $self->error_message("Failed to execute post assemble process: $post_assemble_part");
        return;
    }

    return 1;
}

1;
