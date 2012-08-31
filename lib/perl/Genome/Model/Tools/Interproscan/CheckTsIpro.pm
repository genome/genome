package Genome::Model::Tools::Interproscan::CheckTsIpro;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use File::Slurp;

my $VERSION = '$Revision$';

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'list' => {
            is  => 'String',
            doc => "list of transcript names",
        },

    ],
);

sub help_brief
{
    "utility for checking if a transcript's protein has been run thru interpro and loaded into the MG database";
}

sub execute
{
    my $self = shift;

    my @transcript_names = read_file($self->list);
    chomp @transcript_names;
    $self->check_transcript_names(\@transcript_names);
    return 1;
}

sub check_transcript_names
{
    my $self = shift;
    my $names = shift;
    eval qq|
        use MPSampleData::Transcript;
        use MPSampleData::IproGeneTranscriptXref;
    |;
    die $@ if $@;

    foreach my $name (@$names)
    {
        my ($t) = MPSampleData::Transcript->search( transcript_name => $name );
        if ( !defined($t) )
        {
            print STDERR "transcript name ", $name, " not present\n";
            next;
        }
        my @xref =
            MPSampleData::IproGeneTranscriptXref->search( 
                       transcript_id => $t );
        if ( $#xref == -1 )
        {
            print STDERR $name, "not run or doesn't have any known domains\n";
        }
    }

    return 1;
}


1;

# $Id$
