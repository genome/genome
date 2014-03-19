package Genome::Model::Build::Error;

use strict;
use warnings;

use Genome;

use Data::Dumper;
require Text::Wrap;

class Genome::Model::Build::Error {
    is => 'UR::Object',
    has => [
    # Build event
    build_event => {
        is => 'Genome::Model::Event',
        id_by => 'build_event_id',
        doc => 'The main build event.',
    },
    build_event_id => {
        is => 'Text',
        doc => 'The main build event id.',
    },
    # Stage
    stage_event_id => {
        is => 'Text',
        doc => 'The event id of the stage.',
    },
    stage => {
        is => 'Text',
        doc => 'The name of the stage.',
    },
    # Step event
    step_event => {
        is => 'Genome::Model::Event',
        id_by => 'step_event_id',
        doc => 'The step event.',
    },
    step_event_id => {
        is => 'Text',
        doc => 'The event id of the step.',
    },
    step => {
        is => 'Text',
        doc => 'The name of the step.',
    },
    # Error stuff
    error => {
        is => 'Text',
        doc => 'Error message text.',
    },
    error_wrapped => {
        is => 'Text',
        calculate_from => [qw/ error /],
        calculate => q{
        local $Text::Wrap::columns = $_[1] || 70;
        return Text::Wrap::wrap('', '', $error);
        },
    }
    ],
};

#< Error Log >#
sub error_log {
    my $self = shift;

    my $step_event = $self->step_event;
    unless ( $step_event ) {
        $self->error_message(
            sprintf(
                "Can't get step event for id (%s) to get error log file.",
                $self->step_event_id,
            )
        );
        return;
    }
    
    return $self->step_event->error_log_file;
}

sub error_log_for_web {
    my $self = shift;

    my $step_event = $self->step_event;
    unless ( $step_event ) {
        $self->error_message(
            sprintf(
                "Can't get step event for id (%s) to get error log file.",
                $self->step_event_id,
            )
        );
        return;
    }
    
    return 'file://'.$self->step_event->error_log_file;
}
#<>#

#< Create >#
sub create { 
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    for my $req (qw/ build_event_id stage stage_event_id step_event_id step error /) {
        next if defined $self->$req;
        $self->error_message("Property ($req) is required to create build error.");
        return;
    }
    
    return $self;
}

sub create_from_workflow_errors {
    my ($class, @wf_errors) = @_;

    # wf errors
    unless ( @wf_errors ) {
        $class->error_message("No workflow errors given to create build errors.");
        return;
    }

    my @errors;
    for my $wf_error ( @wf_errors ) {
        $class->_validate_error($wf_error)
            or return; #bad
        my %error_params = $class->_parse_error($wf_error)
            or next; # ok
        #print Dumper(\%params, $wf_error->path_name);
        
        my $error = $class->create(%error_params);
        unless ( $error ) {
            $class->error_message("Can't create build error from workflow error. See above.");
            return;
        }
        push @errors, $error;
    }

    return @errors;
}

sub _validate_error {
    my ($class, $error) = @_;

    unless ( $error->can('path_name') and $error->path_name ) {
        $class->error_message("No path name found in error: ".$error->id);
        $class->error_message( Dumper($error) );
        return;
    }

    #$class->status_message('Error path name: '.$error->path_name);

    return 1;
}

sub _parse_error {
    my ($class, $error) = @_;

    my %error;
    my @tokens = split(m#/#, $error->path_name);
    unless ( @tokens ) { # bad
        $class->error_message("Error parsing error path name: ".$error->path_name);
        print Dumper($error);
        return;
    }
    
    return unless @tokens == 3; # ok, only looking for errors with 3 parts

    # PATH NAME:
    # '%s all stages/%s %s/%s %s'
    # $build_event_id
    # $stage_event_id, (currently the build_event_id, but if stages get ids...)
    # $stage
    # $step
    # $step_event_id
    @error{qw/ build_event_id /} = split(/\s/, $tokens[0], 2);
    @error{qw/ stage_event_id stage /} = split(/\s/, $tokens[1]);
    @error{qw/ step step_event_id /} = split(/\s/, $tokens[2]);
    $error{error} = $error->error;

    unless($error{step_event_id} =~ /\d+/ and $error{stage_event_id} =~ /\d+/) {
        return; #even though there are three parts, this isn't what we're expecting
    }

    return %error;
}
#<>#

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2009 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$

