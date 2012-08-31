package Genome::Model::Tools::AmpliconAssembly::Assemble;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::AmpliconAssembly::Assemble{
    is => 'Genome::Model::Tools::AmpliconAssembly',
    has_optional => [ 
    assembler => {
        is => 'Text',
        is_optional => 1,
        default_value => __PACKAGE__->default_assembler,
        doc => 'The assembler to use. Currently supported assemblers: '.join(', ', valid_assemblers()).' (default: '.__PACKAGE__->default_assembler.').',
    },
    assembler_params => {
        is => 'Text',
        doc => 'String of parameters for the assembler',
    },
    ],
};
#< Helps >#
sub help_detail {
    return <<EOS;
Assembles the amplicons in an amplicon assembly with the given assembler and assembler parameters (assembler_params).  Currently supported assemblers are: .  The assembler parameters should be a string that would normally be used for a command line program.  Encapsulate the string in qoutes, and list attributes with a minus (-) sign followed by a space, and then the values.  See synopsis for examples.  See above for supported assemblers.
EOS
}

sub help_synopsis {
    return join(
        "\n", 
        map { $_[0]->command_name.' '.$_ }
        (
            "--directory . --assembler phredphrap --assembler-params '-'",

        )
    );
}

#< Assemblers #>
sub valid_assemblers {
    return (qw/ phredphrap /);
}

sub default_assembler {
    return (valid_assemblers)[0];
}

sub valid_assemblers_as_string {
    return join(', ', valid_assemblers());
}

#< Command >#
sub sub_command_sort_position { 30; }

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( grep { $self->assembler eq $_ } valid_assemblers() ) {
        $self->error_message(
            sprintf(
                'Invalid assembler (%s)', 
                $self->assembler,
                $self->valid_assemblers_as_string,
            )
        );
        $self->delete;
        return;
    }
    
    unless ( $self->_get_hashified_assembler_params ) {
        # this may not be necessary for all assemblers, but handle that later
        $self->delete;
        return;
    }

    return $self;
}

sub execute {
    my $self = shift;

    my $amplicons = $self->get_amplicons
        or return;

    my $method = '_assemble_amplicon_with_'.$self->assembler;
    
    for my $amplicon ( @$amplicons ) {
        $self->$method($amplicon)
            or return;
    }
    
    $self->status_message(
        sprintf(
            'Successfully assembled amplicon assembly with "%s."',
            $self->assembler,
        )
    );

    $self->amplicon_assembly->update_additional_properties(
        assembler => $self->assembler,
        assembler_params => $self->assembler_params,
    ) or return;
    
    return 1;
}

#< Assembling >#
sub _assemble_amplicon_with_phredphrap {
    my ($self, $amplicon) = @_;

    my $fasta_file = $amplicon->processed_fasta_file;
    return 1 unless -s $fasta_file; # ok
    my $phred_phrap_params = $self->_get_hashified_assembler_params;
    
    my $phrap = Genome::Model::Tools::PhredPhrap::Fasta->create(
        fasta_file => $fasta_file,
        %$phred_phrap_params,
    );

    unless ( $phrap ) { # bad
        $self->error_message("Can't create phrap command.");
        return;
    }

    $phrap->execute;

    return 1;
}

sub _get_hashified_assembler_params {
    my $self = shift;

    return {} unless $self->assembler_params; # ok
    
    unless ( $self->{_hashified_assembler_params} ) {
        my %params = Genome::Utility::Text::param_string_to_hash( $self->assembler_params );
        unless ( %params ) {
            $self->error_message('Malformed assembler params: '.$self->assembler_params);
            return;
        }
        $self->{_hashified_assembler_params} = \%params;

    }

    return $self->{_hashified_assembler_params};
}

1;

#$HeadURL$
#$Id$
