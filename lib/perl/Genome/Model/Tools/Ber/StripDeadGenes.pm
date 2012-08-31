package Genome::Model::Tools::Ber::StripDeadGenes;

use strict;
use warnings;

use Genome;
use Command;

use Carp;
use English;

use BAP::DB::Sequence;
use BAP::DB::SequenceSet;
use BAP::DB::CodingGene;
use Genome::Sys;

use Data::Dumper;
use IPC::Run qw/ run timeout /;
use Time::HiRes qw(sleep);
use DateTime;
use MIME::Lite;

use File::Slurp;    # to replace IO::File access...

use Cwd;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'sqlitedatafile' => {
            is  => 'String',
            doc => "Name of sqlite output .dat file",
        },
        # do we need to adjust this?
        'sqliteoutfile' => {
            is  => 'String',
            doc => "Name of sqlite output .out file",
            is_optional => 1,
        },
        'keep_original' => {
            is => 'Boolean',
            doc => "do not write new file out over old one",
            is_optional => 1,
            default => 0,
        },
    ]
);

sub help_brief
{
    "remove dead genes from .dat file";
}

sub help_synopsis
{
    return <<"EOS"
EOS
}

sub help_detail
{
    return <<"EOS"
EOS
}

sub execute
{
    my $self          = shift;
    # get genes from db, find the ones that are dead
    # start cleaning out those that are dead
    my @lines = read_file($self->sqlitedatafile);
    print "lines ", scalar(@lines),"\n";
    my %dathash;
    my @non_dead_genes = ( );
GENE:    foreach my $line (@lines) {
        my ($gene_name) = (split(/\t/,$line))[0];
        $gene_name =~ /p(\d)_hybrid/;
        my $phase = $1;
        $gene_name =~ s/(\S+)p5_hybrid\.(\d{1,})$/$1$2/;       

        #print $gene_name,"\n";
        my $gi = BAP::DB::CodingGene->search({gene_name => $gene_name});
        my $gene = $gi->next;
        my @tags = $gene->tag_names;
        foreach my $t (@tags) {
            if($t->tag_name eq 'Dead') {
                #print $gene_name , " is dead\n";
                $self->status_message( $gene_name . " is dead");
                next GENE;
            }
        }
        push(@non_dead_genes, $line);
        # need this at all?
        #$dathash{$gene_name}{line} = $line;
        #$dathash{$gene_name}{phase} = $phase;
    }

    my $dead_count = scalar(@lines) - scalar(@non_dead_genes);
    if($dead_count < 0 ) {
        $self->error_message("negative number of dead genes");
        return 0;
    }

    $self->status_message("have $dead_count dead genes to ignore");
    
    # write out new .dat file.
    unless($self->keep_original) {
        # copy original .dat file to "file".orig
        Genome::Sys->copy_file($self->sqlitedatafile,
                                    $self->sqlitedatafile.".orig");
        $self->status_message("writing out new .dat file");
        write_file($self->sqlitedatafile,@non_dead_genes);
    }
    else {
        $self->status_message("not writing out file");
    }

    return 1;
}


1;
