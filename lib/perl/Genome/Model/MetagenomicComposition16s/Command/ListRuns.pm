package Genome::Model::MetagenomicComposition16s::Command::ListRuns; 

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require Carp;

class Genome::Model::MetagenomicComposition16s::Command::ListRuns {
    is => 'Genome::Model::MetagenomicComposition16s::Command',
};

sub sub_command_category { return; }

sub help_brief { 
    return 'List the unique runs of traces in the build';
}

sub help_detail {
    return <<HELP;
    List the unique runs in a build's chromat_dir. This differs from listing instrument data
    which are tracked in the database.
     
    This command is backward compatible for amplicon assembly models/builds.
HELP
}

sub execute {
    my $self = shift;
    
    my $dbh = Genome::DataSource::GMSchema->get_default_handle or die "No dbh for genome"
        or Carp::confess("No db handle for Genome");

    my @builds = $self->_builds
        or return;;
    
    for my $build ( @builds ) {
        my $chromat_dir = $build->chromat_dir;
        unless ( -d $chromat_dir ) {
            $self->error_message('Build '.$build->id.' does not have a chromat_dir');
            next;
        }
        $self->_open_chromat_dir($chromat_dir);
        my $dh = Genome::Sys->open_directory($chromat_dir);
        unless ( $dh ) {
            $self->error_message("Can't open chomrat_dir ($chromat_dir) for build ".$build->id);
            next;
        }
        $dh->read; $dh->read; # . ..
        $self->{_dh} = $dh;
        my %runs;
        while ( my $traces = $self->_next_batch_or_traces ) {
            my $sth = $dbh->prepare(
                sprintf(
                    'select distinct(r.prep_group_id) from GSC.sequence_read r where r.trace_name in(%s)',
                    join(', ', map { "'$_'" } @$traces)
                )
            )
                or Carp::confess("Can't prepare sequence read statement");
            $sth->execute
                or Carp::confess("Can't execute sequence read statement");
            while ( my ($run) = $sth->fetchrow_array ) {
                $runs{$run} = 1;
            }
            $sth->finish;
        }
        if ( @builds > 1 ) {
            print '>'.join(
                ' ', $build->model->name, $build->model->id, $build->id
            )."\n";
        }
        print join("\n", sort keys %runs)."\n";
        $self->_close_chromat_dh;
    }

    return 1;
}

sub _open_chromat_dir{
    my ($self, $chromat_dir) = @_;

    $self->_close_chromat_dh;
    $self->{_dh} = undef;
    $self->{_dh} = Genome::Sys->open_directory($chromat_dir);
    $self->{_dh}->read; $self->{_dh}->read; # . ..

    return $self->{_dh};
}

sub _close_chromat_dh {
    return $_[0]->{_dh}->close if $_[0]->{_dh};
}

sub _next_batch_or_traces {
    my $self = shift;

    my $cnt = 0;
    my @traces;
    while ( my $trace = $self->{_dh}->read ) {
        $trace =~ s/\.gz//;
        push @traces, $trace;
        last if ++$cnt == 500;
    }

    return ( @traces ? \@traces : undef );
}

1;

#$HeadURL$
#$Id$
