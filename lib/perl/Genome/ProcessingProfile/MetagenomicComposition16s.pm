package Genome::ProcessingProfile::MetagenomicComposition16s;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::MetagenomicComposition16s {
    is => 'Genome::ProcessingProfile::Staged',
    has_param => [
        # About
        amplicon_processor => {
            is => 'Text',
            is_optional => 1,
            doc => 'A string of paramters to process amplicons by',
        },
        sequencing_center => {
            is => 'Text',
            doc => 'Place from whence the reads have come.',
            valid_values => [qw/ gsc broad /],
        },
        sequencing_platform => {
            is => 'Text',
            doc => 'Platform (machine) from whence the reads where created.',
            valid_values => [qw/ sanger 454 solexa /],
        },
        #< Assembler >#
        assembler => {
            is => 'Text',
            is_optional => 1,
            doc => 'Assembler name for assembling the reads.',
            valid_values => [qw/ phred_phrap /],
        },
        assembler_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'A string of parameters to pass to the assembler',
        },
        #< Classifier >#
        classifier => {
            is => 'Text',
            is_optional => 1,
            doc => 'Classifier name for classifing the amplicons.',
            valid_values => [qw/ rdp2-1 rdp2-2 rdp2-3 rdp2-5 /],
        },
        classifier_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'A string of parameters to pass to the classifier.',
        },
        #< chimera detector >#
        chimera_detector => {
            is => 'Text',
            is_optional => 1,
            doc => 'Chimera detector name, chimera slayer or nastier'
        },
        chimera_detector_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'A string of parameters to pass to chimera detector',
        },
    ],
};

#< Create >#
sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_)
        or return;

    # Validate params
    for my $type (qw/ assembler classifier /) { 
        my $method = $type.'_params_as_hash';
        $self->$method; # dies if error
    }
    
    # Validate amplicon processor
    if ( $self->amplicon_processor ) {
        $self->status_message('Validate amplicon processor: '.$self->amplicon_processor);
        if ( not eval { $self->amplicon_processor_commands } ) {
            $self->error_message('Failed to validate amplicon processor: '.$self->amplicon_processor);
            $self->delete;
            return;
        }
        $self->status_message('Validate amplicon processor OK');
    }
    # Validate classifier version
    # TODO

    my $chimera_detector_ok = $self->validate_chimera_detector;
    if ( not $chimera_detector_ok ) {
        $self->delete;
        return;
    }

    return $self;
}

#< BUILDING >#
sub stages {
    return (qw/ one /);
}

sub one_job_classes {
    my $self = shift;
    my @subclasses = (qw/ PrepareInstrumentData Classify Orient /);
    push @subclasses, (qw/ DetectAndRemoveChimeras /) if $self->chimera_detector;
    push @subclasses, (qw/ Reports /);
    return map { 'Genome::Model::Event::Build::MetagenomicComposition16s::'.$_ } @subclasses;
}

sub one_objects {
    return 1;
}

#< Hashify >#
sub _operation_params_as_hash {
    my ($self, $operation) = @_;

    my $method = $operation.'_params';
    my $params_string = $self->$method;
    return unless $params_string; # ok 

    my %params = Genome::Utility::Text::param_string_to_hash($params_string);

    unless ( %params ) { # not ok
        die $self->error_message("Malformed $operation params: $params_string");
    }

    return %params;
}

sub assembler_params_as_hash {
    return $_[0]->_operation_params_as_hash('assembler');
}

sub classifier_params_as_hash {
    return $_[0]->_operation_params_as_hash('classifier');
}

sub amplicon_processor_commands {
    my $self = shift;
    
    my $command_string = $self->amplicon_processor;
    return unless $command_string;
    my @valid_commands;
    my @commands = split( /\|/, $command_string );
    for my $command ( @commands ) {
        $command =~ s/^\s+//;
        $command =~ s/\s+$//;
        $command = "gmt sx $command";
        my $valid = Genome::Model::Tools::Sx::Validate->validate_command( $command );
        if ( not $valid ) {
            die $self->error_message("Invalid amplicon processor command: $command");
        }
        push @valid_commands, $command;
    }

    return @valid_commands;
}

sub validate_chimera_detector {
    my $self = shift;

    my $detector = $self->chimera_detector;
    my $params = $self->chimera_detector_params;
    return 1 if not $params and not $detector;
    $self->status_message('Validate chimera detector...');
    if ( not $params or not $detector ) {
        $self->error_message('Cannot give chimera detector without params or vice versa!');
        return;
    }
    $detector =~ s/_/\-/g;
    $self->chimera_detector($detector);

    my $class = 'Genome::Model::Tools::'.Genome::Utility::Text::string_to_camel_case(join(' ', split('-', $detector))).'::DetectChimeras';
    my $meta = eval{ $class->__meta__; };
    if ( not $meta ) {
        $self->error_message("Invalid chimera detector: $detector");
        return;
    }

    my $cmd = "gmt $detector detect-chimeras $params";
    $self->status_message('Chimera detector command: '.$cmd);
    $cmd .= ' -h 2>&1 > /dev/null'; # add help to check for invalid opts, redirect to dev/null

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd, print_status_to_stderr => 0); };
    if ( not $rv ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to validate chimera detector and params!');
        return;
    }

    $self->status_message('Validate chimera detector...OK');
    return 1;
}

1;

