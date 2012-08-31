package Genome::Model::Tools::Hgmi::Finish;


use strict;
use warnings;

use Genome;
use Command;
use Carp;
use DateTime;
use List::MoreUtils qw/ uniq /;
use IPC::Run qw/ run /;
use Data::Dumper;

use BAP::DB::Organism;
use BAP::DB::SequenceSet;
use BAP::DB::Sequence;

use Bio::SeqIO;
use Cwd;
use English;

class Genome::Model::Tools::Hgmi::Finish {
    is => 'Command',
    has => [
        organism_name => {
            is => 'String',
        },
        locus_tag => {
            is => 'String',
        },
        project_type => {
            is => 'String',
        },
    ],
    has_optional => [ 
        nr_db => {
            is => 'String',
            doc => "path to nr seq db",
            default => "/gscmnt/ams1102/info/annotation/blastdb/gsc_bacterial/bacterial_nr/bacterial_nr"
        },
        locus => {
            is => 'String',
        },
        gram_stain => {
            is => 'String',
        },
        ncbi_taxonomy_id => {
            is => 'String',
        },
        dev => {
            is => 'Boolean',
        },
        sequence_set_id => { 
            is => 'Integer',
            doc => "sequence set id" ,
        },
        acedb_version => { 
            is => 'String',
            doc => "Ace DB version (V1,V2,etc)"
        },
        sequence_set_name => {
            is => 'String',
        },
        work_directory => {
            is => 'String',
        }, 
        script_location => {
            is => 'String',
            doc => "path or name of bap finish project script",
            default => "bap_finish_project",
        },
        skip_acedb_parse => {
            is => 'Boolean',
            doc => "skip parsing into acedb for testing",
		    #default => 0, # this skips parsing into acedb...
        },
        assembly_name => {
            is  => 'String',
            doc => "assembly name from the config file",
        },
        org_dirname => {
            is  => 'String',
            doc => "organism directory name from config file",
        },
        assembly_version => {
            is  => 'String',
            doc => "assembly version from config file",
        },
        pipe_version => {
            is  => 'String',
            doc => "pipe version from config file",
        },
        path => {
            is  => 'String',
            doc => "path from config file",
        },
    ]
};

sub help_brief
{
    "tool for running the finish step (bap_finish_genes).  Not finished."
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
Runs the finish step (bap_finish_genes).
EOS

}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
Runs the finish step (bap_finish_genes).
EOS

}

sub execute
{
    my $self = shift;

    my %params = $self->gather_details();

    my $rv = Genome::Model::GenePrediction::Command::Bacterial::Finish->execute(%params);

    unless($rv) {
        $self->error_message("can't run finish project step");
        return 0;
    }
    return 1;
}

sub gather_details
{
    my $self = shift;

    my $organism_name    = $self->organism_name;
    my $locus_tag        = $self->locus_tag;
    my $project_type     = $self->project_type;
    my $assembly_name    = $self->assembly_name;
    my $org_dirname      = $self->org_dirname;
    my $assembly_version = $self->assembly_version;
    my $pipe_version     = $self->pipe_version;
    my $path             = $self->path;

    my ($ncbi_taxonomy_id, $gram_stain, $locus, $organism_id, );
    if (defined($self->dev)) { $BAP::DB::DBI::db_env = 'dev'; }
    my $organism_obj = BAP::DB::Organism->retrieve('organism_name'=> $organism_name);

    unless (defined($organism_obj)) 
    {
        croak " organism_obj is not set - are you running this before the predict and merge steps? Finish.pm\n\n";
    }

    $organism_name     = $organism_obj->organism_name();
    $organism_id       = $organism_obj->organism_id();
    $ncbi_taxonomy_id  = $organism_obj->ncbi_taxonomy_id();
    $gram_stain        = $organism_obj->gram_stain();
    $locus             = $organism_obj->locus();
    my @cols = ($organism_id, $organism_name, 
                $ncbi_taxonomy_id, $gram_stain, $locus);
    @cols = map { defined($_) ? $_ : 'NULL' } @cols;
    
    print join("\t", @cols), "\n\n";
    my $cwd;
    if(defined($self->work_directory))
    {
        $cwd = $self->work_directory;
    }
    else
    {
        $cwd = getcwd();
    }
    my $sequence_set_name = $assembly_name;
    my $sequence_set_name_obj;
    my $sequence_set_obj;
    my $sequence_set_id = $self->sequence_set_id;

    $sequence_set_name_obj = BAP::DB::SequenceSet->retrieve('sequence_set_name'=> $sequence_set_name);

    unless(defined($sequence_set_name_obj))
    {
        croak "nothing found for $sequence_set_name,\nperhaps you are running this before the predict and merge steps?Finish.pm\n\n ";
    }
    else
    {
        print "Sequence-set-name: '$sequence_set_name' already ready!! Here is your information:\n\n";
    }

    unless (defined($sequence_set_id))
    {
        $sequence_set_id = $sequence_set_name_obj->sequence_set_id();
    }

    my @list =($organism_id,$sequence_set_name, $sequence_set_id);
    print join("\t",@list),"\n\n";

    my $bapfinish_output = $cwd."/".$locus_tag."_bfp_BAP_screenoutput".$sequence_set_id.".txt";

    my %params = (
        'sequence_set_id' =>  $sequence_set_id,
        'locus_id' => $locus_tag,
        'project_type' => $self->project_type,
        'acedb_version' => $self->acedb_version,
		'assembly_name' => $assembly_name,
		'org_dirname' => $org_dirname,
		'assembly_version' => $assembly_version,
		'pipe_version' => $pipe_version,
		'path_base' => $path,
    );

    $params{dev} = $self->dev;
    $params{no_acedb} = $self->skip_acedb_parse;
    print Dumper(\%params),"\n";
    return %params;
}



1;
