package Genome::Model::Tools::Hgmi::Merge;


use strict;
use warnings;

use Genome;
use Command;
use Carp;
use DateTime;
use List::MoreUtils qw/ uniq /;
use IPC::Run qw/ run /;
use File::Slurp;
use Cwd;
require "pwd.pl";
use English;

use BAP::DB::Organism;
use BAP::DB::SequenceSet;
use BAP::DB::Sequence;
use Data::Dumper;
use YAML qw( DumpFile LoadFile );

use Bio::SeqIO;

class Genome::Model::Tools::Hgmi::Merge (
    is => 'Command',
    has => [
        organism_name => {
            is => 'String',
			doc => 'Organism\'s latin name',
		},
        locus_tag => {
            is => 'String',
			doc => 'Locus tag for project, containing DFT/FNL postpended',
        },
        project_type => {
            is => 'String',
			doc => 'Type of project',
		},
    ],
    has_optional => [
        nr_db => {
            is => 'String',
			doc => 'Path to the default non-redundant sequence database, only used if local copy unavailable',
			default => '/gscmnt/ams1102/info/annotation/blastdb/gsc_bacterial/bacterial_nr/bacterial_nr',
		},
        iprpath => {
            is => 'String',
            doc => "specify different version of iprscan",
            default => '/gsc/scripts/bin/iprscan',
        },
        locus_id => {
            is => 'String',
			doc => 'locus_tag with OUT DFT/FNL postpended',
        },
        gram_stain => {
            is => 'String',
			doc => 'gram stain for bacteria (positive or negative)',
		},
        ncbi_taxonomy_id => {
            is => 'String',
            doc => 'NCBI taxonomy id',
		},
        dev => {
            is => 'Boolean',
			doc => 'use development db',
		},
        work_directory => {
            is => 'String',
			doc => 'work directory',
		},
        sequence_set_id => { 
            is => 'Integer',
			doc => 'sequence set id',
		},
		runner_count => { 
            is => 'Integer',
			doc => 'Number of job runners created',
			default => 50,
		},
    ],
);

sub help_brief {
    return 'Runs the gene merging step of the gene prediction process';
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
Runs the gene merging step of the gene prediction process.
EOS
}

sub help_detail {
    my $self = shift;
    return <<"EOS"
Runs the gene merging step of the gene prediction process.
EOS
}


sub execute {
    my $self = shift;

    my %params = $self->gather_details();
    $self->status_message("Gathered params, now running merge:\n" . Data::Dumper::Dumper(\%params));
    my $rv = Genome::Model::GenePrediction::Command::Bacterial::Merge->execute(%params);

    unless($rv) {
        confess 'Could not run merge step!';
    }

    # executed successfully
    $self->validate(\%params);

    return 1;
}

sub gather_details
{
    my $self = shift;

    my $organism_name  = $self->organism_name;
    my $locus_tag      = $self->locus_tag;
    my $project_type   = $self->project_type;
    my $runner_count   = $self->runner_count;
    my $nr_db          = $self->nr_db;

    my ($ncbi_taxonomy_id, $gram_stain, $locus, $organism_id, );

    if (defined($self->dev)) { $BAP::DB::DBI::db_env = 'dev'; }
    my $organism_obj = BAP::DB::Organism->retrieve('organism_name'=> $organism_name);

    if (defined($organism_obj)) { 
        print "\n$organism_name, already exist!  Here is your information from Merge.pm:\n\n";
    }
    else
    {
        croak "$organism_name not in db, maybe you need to run bap_predict_genes.pl, first? Merge.pm\n";
    }

    $organism_name     = $organism_obj->organism_name();
    $organism_id       = $organism_obj->organism_id();
    $ncbi_taxonomy_id  = $organism_obj->ncbi_taxonomy_id();
    $gram_stain        = $organism_obj->gram_stain();
    $locus             = $organism_obj->locus();

    my @cols = ($organism_id, $organism_name, $ncbi_taxonomy_id, $gram_stain, $locus);
    @cols = map { defined($_) ? $_ : 'NULL' } @cols;

    print join("\t", @cols), "\n\n";

    my $cwd = $self->cwd();
    my @cwd = split(/\//x,$cwd);

    my ($sequence_set_name, $analysis_version_num, $hgmi_sequence_dir);

    # FIXME This REALLY needs to go away... a seemlingly minor change somewhere else to directory
    # structure could totally screw this up. There's no reason this information can't just be passed
    # to this tool instead of trying to deduce it in this manner.

    # these below are dangerous, should be altered.
    if ($project_type =~ /HGMI/x )
    {
        # these need to be based on the directory structure,
        # instead of just a 'raw' split.
        #unless($#cwd == 9)
        #{
        #    croak "directory structure is wrong in HGMI!?!? Merge.pm";
        #}
        #$sequence_set_name    = $cwd[6]; #HGMI projects
        $sequence_set_name    = $cwd[-4]; #HGMI projects
        #$analysis_version_num = $cwd[9]; #HGMI projects
        $analysis_version_num = $cwd[-1]; #HGMI projects
        $hgmi_sequence_dir    = join("\/", @cwd,'Sequence',$locus_tag); #HGMI projects
        
    }
    else # HMPP/Enterobacter
    {
        # this might be wrong and should be changed to something similar
        # to the above...
        unless($#cwd == 10)
        {
            croak "directory structure is wrong in HMPP !?!?!? Merge.pm";
        }
        $sequence_set_name = $cwd[7]; #HMPP and Enterobacter
        $analysis_version_num = $cwd[10]; #HMPP and Enterobacter
        $hgmi_sequence_dir = join("\/", @cwd[0..10],'Sequence',$locus_tag); #HMPP and Enterobacter
        
    }

    unless (defined($sequence_set_name)) 
    {
        croak " sequence_set_name is not set in Merge.pm! ";
    }

    my $sequence_set_name_obj;
    my $sequence_set_obj;
    my $sequence_set_id;

    $sequence_set_name_obj = BAP::DB::SequenceSet->retrieve('sequence_set_name'=> $sequence_set_name);

    if (defined($sequence_set_name_obj)) 
    {
        print "Sequence-set-name: '$sequence_set_name' already exist!! Here is your information from Merge.pm:\n\n";

    }
    else
    {
        croak "sequence set name not setup, maybe you need to run the bap_predict_genes.pl step? Merge.pm\n";
    }

    unless (defined($sequence_set_id))
    {
        $sequence_set_id = $sequence_set_name_obj->sequence_set_id();
    }

    $self->sequence_set_id($sequence_set_id);

    my @list =($organism_id,$sequence_set_name, $sequence_set_id);

    print join("\t",@list),"\n\n";
   
    if (!defined($nr_db)) 
    {
        croak "nr-db: \$nr_db is empty! Merge.pm";
    }
    elsif ( ! -f $nr_db )
    {
        croak "nr-db: ". $nr_db . " doesn't exist! Merge.pm";
    }
    
    my $bmg_job_stdout       = $cwd."/".$locus_tag."_bmg_BAP_job_".$sequence_set_id.".txt";
    my $bmg_job_stderr       = $cwd."/".$locus_tag."_bmg_BAP_job_err_".$sequence_set_id.".txt";
    my $debug_file           = $cwd."/".$locus_tag."_bmg_debug_file_".$sequence_set_id.".txt";
    my $bapmergegenes_output = $cwd."/".$locus_tag."_bmg_BAP_screenoutput_".$sequence_set_id.".txt";
    
    my %params = (
                   'sequence_set_id' => $sequence_set_id,
                   'job_stdout' => $bmg_job_stdout,
                   'job_stderr' => $bmg_job_stderr,
                   'runner_count' => $runner_count,
                   'debug_file' => $debug_file,
                   'nr_db' => $self->nr_db,
                   );

    if(defined($self->dev)) { $params{dev} = 1; }
    
    print Dumper(\%params),"\n";
    return %params;
}

sub _cwd {

    my ($self) = @_;

    my $cwd;
    
    if (defined($self->work_directory)) {
        $cwd = $self->work_directory;
    }
    else {
        $cwd = getcwd();
    }

    return $cwd;
    
}

sub validate {
    my $self = shift;
    my $params = shift;
    my %valid = ();
    %valid = $self->cmdline_to_hash();
    $valid{params} = %$params;
    my $cwd = getcwd();
    if(-f "$cwd/merge_valid") {
        unlink("$cwd/merge_valid");
    }
    #DumpFile("$cwd/merge_valid", $params);
    DumpFile("$cwd/merge_valid", \%valid);
    #write_file("$cwd/merge_valid","merge successful");
    return 1;
}


sub is_valid {
    my $self = shift;
    use Data::Compare;
    my $cwd = getcwd();
    if (-f "$cwd/merge_valid" ) {
        my %params = $self->gather_details();
        my %cmdline = $self->cmdline_to_hash();
        my %to_compare = %cmdline;
        $to_compare{params} = %params;
        my %validation = %{LoadFile("$cwd/merge_valid")};
        print Dumper(\%to_compare,\%validation),"\n";
        if(Compare(\%to_compare,\%validation)) {
            $self->status_message("previous merge run successful");
            return 1;
        }
        else
        {
            $self->status_message("$cwd/merge_valid exists, but doesn't match params");
            return 0;
        }

    }
    return 0;
}

sub cmdline_to_hash {
    my $self = shift;
    my %cmdline = ( );
    $cmdline{locus_tag} = $self->locus_tag;
    $cmdline{organism_name} = $self->organism_name;
    $cmdline{project_type} = $self->project_type;
    $cmdline{dev} = $self->dev;
    $cmdline{gram_stain} = $self->gram_stain;
    $cmdline{locus_id} = $self->locus_id;
    $cmdline{ncbi_taxonomy_id} = $self->ncbi_taxonomy_id;
    $cmdline{nr_db} = $self->nr_db;
    $cmdline{runner_count} = $self->runner_count;
    $cmdline{sequence_set_id} = $self->sequence_set_id;
    return %cmdline;
}

1;

# $Id$
