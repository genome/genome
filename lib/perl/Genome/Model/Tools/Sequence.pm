package Genome::Model::Tools::Sequence;

use strict;
use warnings;
use Genome;
use IO::File;

class Genome::Model::Tools::Sequence{
    is => 'Command',
    has =>[
        chromosome => {
            is => 'Text',
            doc => "chromosome for desired sequence",
        },
        start => {
            is => 'Number',
            doc => 'sequence start position',
        },
        stop => {
            is => 'Number',
            doc => 'sequence stop position',
        },
        build_id => {
            is => 'Number',
            is_optional => 1,
            doc => "build id to grab sequence from, defaults to hb36",
        },
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            is_optional => 1,
        },
        species => {
            is => 'Text',
            is_optional => 1,
            doc => 'use reference for this species',
            default => 'human',
            valid_values => ['human', 'mouse'],
        },
        version => {
            is => 'Text',
            is_optional => 1,
            doc => 'reference version to be used',
            default => '36',
        },
        sequence =>{
            is => 'SCALAR',
            is_output => 1,
            doc => "This is populated with the sequence returned from running this command",
            is_optional => 1,
        },
    ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "get a seqeunce from an ImportedReferenceSequence model";
}

sub help_synopsis {
    return <<"EOS"
Given chromosome start and stop, returns a sequence from a ImportedReferenceSequence
EOS
}

sub help_detail {
    return <<"EOS"
Given chromosome start and stop, returns a sequence from a ImportedReferenceSequence
EOS
}

# For speed.  slightly unsafe, but there shouldn't ever be any problems, right?
sub __errors__() {
    return;
}

sub execute {
    my $self = shift;
    my $start = $self->start;
    $start = 1 if $start < 1;
    my $results =
        eval { 
            if ($self->build){
                lookup_sequence(
                        chromosome => $self->chromosome,
                        start => $start,
                        stop => $self->stop,
                        build => $self->build,
                        species => $self->species,
                        version => $self->version);
            }
            else{
                lookup_sequence(
                        chromosome => $self->chromosome,
                        start => $start,
                        stop => $self->stop,
                        build => '',
                        species => $self->species,
                        version => $self->version);
            }
        };
    if ($@) {
        $self->error_message($@);
        return 0;
    } else {
        $self->sequence($results);
        print $self->sequence, "\n"; #TODO: should this be here?
        return 1;
    }
}

sub lookup_sequence{
    my %args = @_;
    my($chromosome, $start, $stop, $build, $species, $version) = @args{'chromosome','start','stop','build','species','version'};

    unless ($build){
        my $model = Genome::Model->get(name => "NCBI-" . $species);
        unless ($model) {
            Carp::croak("Could not get imported reference model for " . $species);
            return 0;
        }

        $build = $model->build_by_version($version);
        unless ($build) {
            Carp::croak("Could not get imported reference version " . $version . " for " . $species);
            return 0;
        }
    }

    my $seq = $build->sequence($chromosome, $start, $stop);
    return $seq;
}

1;
