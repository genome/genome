package Genome::Model::Tools::AmpliconAssembly::TrimAndScreen;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::AmpliconAssembly::TrimAndScreen {
    is => 'Genome::Model::Tools::AmpliconAssembly',
    has => [
    trimmer_and_screener => {
        is => 'Text',
        doc => 'The trimmer and screener to use.  Currently supported trimmers: '.join(', ', valid_trimmers_and_screeners()),
    },
    ],
    has_optional => [
    trimmer_and_screener_params => {
        is => 'Text',
        default_value => '',
        doc => 'String of parameters for the trimmer and screener',
    },
    ],
};

#< Valid Trimmers and Screeners >#
sub valid_trimmers_and_screeners {
    return (qw/ trim3_and_crossmatch lucy /);
}
#< Helps >#
sub help_detail {
    return <<EOS;
Runs quality trimming and vector screening. See above for valid trimmers and screeners.
EOS
}

sub help_synopsis {
}

#< Command >#
sub sub_command_sort_position { 20; }

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->trimmer_and_screener ) {
        $self->error_message("No trimmer and screener was given");
        $self->delete;
        return;
    }

   unless ( grep { $self->trimmer_and_screener eq $_ } valid_trimmers_and_screeners() ) {
        $self->error_message("Invalid trimmer and screener: ". $self->trimmer_and_screener);
        $self->delete;
        return;
    }

    unless ($self->_set_trim_and_screen_params ) {
        $self->error_message("Can't set the trim and screen parameters");
        $self->delete;
        return;
    }

    return $self;
}
        
sub execute {
    my $self = shift;

    my $amplicons = $self->get_amplicons
        or return;

    my $trim_and_screen_method = '_trim_and_screen_amplicon_by_'.$self->trimmer_and_screener;
    for my $amplicon ( @$amplicons ) {
        $self->$trim_and_screen_method($amplicon)
            or return;
    }
    
    $self->amplicon_assembly->update_additional_properties(
        trimmer_and_screener => $self->trimmer_and_screener,
        trimmer_and_screener_params => $self->trimmer_and_screener_params,
    ) or return;

    return 1;
}

#< Trim and SCreen Params >#
sub _set_trim_and_screen_params {
    my $self = shift;

    my ($trim_string, $screen_string) = split(
        m#;\s?#, $self->trimmer_and_screener_params
    );

    $self->{_trim_params} = {};
    if ( $trim_string ) {
        $self->{_trim_params} = Genome::Utility::Text::param_string_to_hash($trim_string)
            or return;
    }

    $self->{_screen_params} = {};
    if ( $screen_string ) {
        $self->{_screen_params} = Genome::Utility::Text::param_string_to_hash($screen_string)
            or return;
    }

    return 1;
}

sub _trim_params {
    my $self = shift;

    unless ( $self->{_trim_params} ) {
        $self->_set_trim_and_screen_params;
    }

    return %{$self->{_trim_params}};
}

sub _screen_params {
    my $self = shift;

    unless ( $self->{_screen_params} ) {
        $self->_set_trim_and_screen_params;
    }

    return %{$self->{_screen_params}};
}

#< Trim3 and CM >#
sub _trim_and_screen_amplicon_by_trim3_and_crossmatch {
    my ($self, $amplicon) = @_;

    return 1 unless -s $amplicon->fasta_file; # ok
    
    my $trim3 = Genome::Model::Tools::Fasta::Trim::Trim3->create(
        fasta_file => $amplicon->fasta_file,
        $self->_trim_params,
    );
    unless ( $trim3 ) { # not ok
        $self->error_message("Can't create trim3 command for amplicon: ".$amplicon->name);
        return;
    }
    $trim3->execute; # ok
   
    return 1 unless -s $amplicon->fasta_file; # ok

    my $screen = Genome::Model::Tools::Fasta::ScreenVector->create(
        fasta_file => $amplicon->fasta_file,
        $self->_screen_params,
    );
    unless ( $screen ) { # not ok
        $self->error_message("Can't create screen vector command for amplicon: ".$amplicon->name);
        return;
    }
    $screen->execute; # ok
    
    return 1;
}

#< Lucy >#
sub _trim_and_screen_amplicon_by_lucy {
    my ($self, $amplicon) = @_;

    return 1 unless -s $amplicon->fasta_file; # ok

    my $lucy = Genome::Model::Tools::Fasta::Trim::Lucy->create(
        fasta_file => $amplicon->fasta_file,
        $self->_trim_params,
    );
    unless ( $lucy ) {
        $self->error_message("Can't create lucy command for amplicon: ".$amplicon->name);
        return;
    }
    $lucy->execute;

    return 1;
}

1;

#$HeadURL$
#$Id$
