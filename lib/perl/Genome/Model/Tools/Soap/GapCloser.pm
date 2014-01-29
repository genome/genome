package Genome::Model::Tools::Soap::GapCloser;

use strict;
use warnings;

use Genome;
use File::Basename;

class  Genome::Model::Tools::Soap::GapCloser {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
        version => {
            is => 'String',
            doc => 'Version of GapCloser',
            valid_values => [qw/ 1.10 /], # only one, and it's deployed to /gsc/scripts/bin
        },
        assembly_directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'Assembly directory to derive input/output files from.',
        },
        a => {
            is => 'Text',
            is_optional => 1,
            doc => 'SOAP generated scafSeq file. Default is named "*.scafSeq" in the assembly directory',
        },
        b => {
            is => 'Text',
            is_optional => 1,
            doc => 'Config file used for the SOAP assembly. Default is named "config_file" in the assembly directory',
        },
        o => {
            is => 'Text',
            is_optional => 1,
            doc => 'GapCloser output fasta file name. It will be in the assembly directory. Default is named "gapfill".',
        },
        p => {
            is => 'Number',
            doc => 'Overlap length/K value. Typical default is 25. Max is 31.',
        },
    ],
};

sub help_brief {
    'GapCloser: close dem gaps!';
}

sub help_detail {
    return <<HELP;
HELP
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my @input_file_methods = (qw/ a b /);

    my $p = $self->p;
    if ( $p > 31 or $p < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ p /],
            desc => "The p ($p) must be an integer between 1 and 31!",
        );
    }

    if ( $self->assembly_directory ) {
        if ( not -d $self->assembly_directory ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => 'The assembly_directory is not a directory!',
            );
            return @errors;
        }
        if ( not -d $self->assembly_directory.'/edit_dir' ) {
            Genome::Sys->create_directory( $self->assembly_directory.'/edit_dir' );
        }
        $self->a( $self->_resolve_scaffold_sequence_file ) unless $self->a;
        $self->b( $self->_resolve_config_file ) unless $self->b;
        $self->o( $self->assembly_directory.'/edit_dir/gapfill' ) unless $self->o;
    }
    elsif ( not $self->o ) { 
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ o /],
            desc => 'No o given and no assembly_directory given to determine the output file!',
        );
    }

    for my $input_file_method ( @input_file_methods ) {
        my $file = $self->$input_file_method;
        if ( not $file ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [ $input_file_method ],
                desc => "Parameter ($input_file_method) is required! This can be resolved from the assmbly directory or passed in.",
            );
            return @errors;
        }
        if ( not -s $file ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [ $input_file_method ],
                desc => "File $file ($input_file_method) does not have any size!",
            );
            return @errors;
        }
    }

    return @errors;
}

sub execute {
    my $self = shift;

    $self->debug_message('SOAP GapCloser...');

    unlink $self->o;

    my $cmd = sprintf(
        'GapCloser -o %s -a %s -b %s -p %s',
        $self->o,
        $self->a,
        $self->b,
        $self->p,
    );
    $self->debug_message("Run GapCloser command: $cmd");
    my $rv = eval{ Genome::Sys->shellcmd( cmd => $cmd ) };
    if ( $rv ) {
        $self->error_message('GapCloser shell command failed!');
        return;
    }
    $self->debug_message('Run GapCloser command...OK');

    my $output = $self->o;
    if ( not -s $output ) {
        $self->error_message("GapCloaser ran ok, but output file ($output) was not created!");
        return;
    }
    $self->debug_message("Output file exists: $output");

    $self->debug_message('SOAP GapCloser...DONE');

    return 1;
}

1;

