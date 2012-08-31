package Genome::Model::Tools::Hgmi::Predict;


use strict;
use warnings;

use Genome;
use Command;
use Carp;
use File::Slurp;
use IO::Dir;
use DateTime;
use List::MoreUtils qw/ uniq /;
use IPC::Run qw/ run /;
use Data::Dumper;

use BAP::DB::Organism;
use BAP::DB::SequenceSet;
use BAP::DB::Sequence;

use Bio::SeqIO;

use Cwd;
#require "pwd.pl";
use Data::Dumper;
use English;
use YAML qw( LoadFile DumpFile );


UR::Object::Type->define(
                         class_name => __PACKAGE__,
                         is  => 'Command',
                         has => [
                                 'organism_name'    => {is => 'String',
							doc => "Genius_species" ,
						       },
                                 'locus_tag'        => {is => 'String',
							doc => "Locus tag for project, containing DFT/FNL postpended",
						       },
                                 'project_type'     => {is => 'String',
							doc => "", 
						       },
                                 'locus_id'         => {is => 'String',
							doc => "locus_tag with OUT DFT/FNL postpended",
						       },
                                 'gram_stain'       => {is => 'String',
							doc => "",
							valid_values => ['positive','negative','none'],
						       },
                                 'ncbi_taxonomy_id' => {is => 'String',
                                                        doc => "",
                                                        is_optional =>1,
						       },
                                 'work_directory'   => { is => 'String',
							 doc => "",
							 is_optional => 1,
						       },
                                 'dev'              => {is => 'Boolean',
							doc => "use development db",
							is_optional => 1,
						       },
                                 'script_location'  => {is => 'String',
							doc =>"location for bap_predict_genes",
							is_optional => 1,
							default => 'bap_predict_genes',
						       },
				 'runner_count'     => { is => 'Number',
							 doc => "runner_count",
							 default => 10,
							 is_optional => 1,
						       },
                                 ]
                         );



sub help_brief
{
    "tool for running the predict step (bap_predict_genes). Not done yet.";
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
Runs predict step (bap_predict_genes).
Not done yet.
EOS

}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
Runs predict step (bap_predict_genes).
Gathers run time paramters for execution.

EOS

}

sub execute
{
    my $self = shift;

    #my @predict_command = $self->gather_details();
    my %params = $self->gather_details();
    #IPC::Run::run( @predict_command ) || croak "can't run Predict.pm : $CHILD_ERROR";

    #if($@)
    #{
    #    carp "uh, we should have quit by now";
    #    exit 1;
    #}
    my $rv = Genome::Model::GenePrediction::Command::Bacterial::Predict->execute(%params);
    unless($rv) {
        $self->error_message("can't run prediction step");
        return 0;
    }

    $self->validate(\%params);

    return 1;
}


sub gather_details
{
    my $self = shift;
    my $organism_name    = $self->organism_name;
    my $locus_tag        = $self->locus_tag;
    my $project_type     = $self->project_type;
    my $ncbi_taxonomy_id = $self->ncbi_taxonomy_id;
    my $gram_stain       = $self->gram_stain;
    my $locus_id         = $self->locus_id;
    my ($organism_id );

    if (defined($gram_stain)) {

      if ($gram_stain =~ /positive/){

	$gram_stain = "+";

      }
      else {

	$gram_stain = "-";

      }
    }

    if (defined($self->dev)) { $BAP::DB::DBI::db_env = 'dev'; }
    my $organism_obj = BAP::DB::Organism->retrieve('organism_name'=> $organism_name);

    if (defined($organism_obj)) { 
        print "\n$organism_name, already exist!  Here is your information:\n\n";
        
    }
    else
    {   
        $organism_obj = BAP::DB::Organism->insert({
                                                      'organism_name'    => $organism_name,
                                                      'ncbi_taxonomy_id' => $ncbi_taxonomy_id,
                                                      'gram_stain'       => $gram_stain,
                                                      'locus'            => $locus_id,
                                                  }
                                                  );
        
    }

    $organism_name     = $organism_obj->organism_name();
    $organism_id       = $organism_obj->organism_id();
    $ncbi_taxonomy_id  = $organism_obj->ncbi_taxonomy_id();
    $gram_stain        = $organism_obj->gram_stain();
    $locus_id          = $organism_obj->locus();
    my @cols = ($organism_id, $organism_name, $ncbi_taxonomy_id, $gram_stain, $locus_id);
    @cols = map { defined($_) ? $_ : 'NULL' } @cols;
    
    print join("\t", @cols), "\n\n";
    

    BAP::DB::DBI->dbi_commit();

    unless (defined($organism_obj)) 
    {
        croak " organism_obj is not set in Predict.pm ! ";
    }

    # cwd should look like:
    # /gscmnt/278/analysis/HGMI/B_catenulatum/Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb/Version_1.0/BAP/Version_1.0
    my $cwd;
    my @cwd;
    if(!defined($self->work_directory))
    {
        $cwd = getcwd();

    }
    else
    {
        $cwd = $self->work_directory;
    }

    @cwd = split(/\//x,$cwd);

    my ($sequence_set_name, $analysis_version_num, $hgmi_sequence_dir);

    # FIXME:
    # these below are dangerous
    # most of this looks like it is easy to index from the end of the string.
    if (($project_type =~ /HGMI/x)  )
    {
        # these need to be based on the directory structure,
        # instead of just a 'raw' split.
#        unless($#cwd == 9)
#        {
#            croak "directory structure is wrong or broken\n$cwd\n$#cwd  Predict.pm \n";
#        }
        #$sequence_set_name = $cwd[6]; #HGMI projects
        $sequence_set_name = $cwd[-4]; #HGMI projects
        #$analysis_version_num = $cwd[9]; #HGMI projects
        $analysis_version_num = $cwd[-1]; #HGMI projects
        $hgmi_sequence_dir = join("\/", $cwd, 'Sequence',$locus_tag);
#        $hgmi_sequence_dir = join("\/", @cwd[0..9],'Sequence',$locus_tag); #HGMI projects
        $self->status_message("hgmi seq dir: ". $hgmi_sequence_dir); 
    }
    else # HMPP/Enterobacter
    {
#        unless($#cwd == 10)
#        {
#            $self->error_message("directory structure is too long or short? ".$#cwd . " components");
#            croak "directory structure is wrong or broken in Predict.pm\n";
#        }
        $sequence_set_name = $cwd[-4];
        $analysis_version_num = $cwd[-1];
        #$hgmi_sequence_dir = join("\/", @cwd[0..10],'Sequence',$locus_tag); 
        $hgmi_sequence_dir = join("\/", $cwd,'Sequence',$locus_tag); 
        $self->status_message("sequence dir: ". $hgmi_sequence_dir); 
    }

    unless (defined($sequence_set_name)) 
    {
        croak " sequence_set_name is not set from Predict.pm! \n";
    }

    my $sequence_set_name_obj;
    my $sequence_set_obj;
    my $sequence_set_id;

    $sequence_set_name_obj = BAP::DB::SequenceSet->retrieve('sequence_set_name'=> $sequence_set_name);

    if (defined($sequence_set_name_obj)) 
    {
        print "Sequence-set-name: '$sequence_set_name' already exist!! Here is your information from Predict.pm:\n\n";

    }
    else 
    {
        my $short_ver_num;
        $short_ver_num = $analysis_version_num;
        $short_ver_num =~ s/Version_(\d)\.0/v$1/x;

        my $fasta_file = join(".",$hgmi_sequence_dir,$short_ver_num,'contigs','newname', 'fasta');
    
        unless (defined($fasta_file)) 
        {
            croak "fasta-file: '$fasta_file' is not set! Predict.pm\n";
        }
    
        if ((-z $fasta_file) ) 
        {
            croak "fasta-file: '$fasta_file 'is empty or non-existant! Predict.pm\n";
        }


        my $in = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');
    
        $sequence_set_obj = BAP::DB::SequenceSet->insert({
            sequence_set_name => $sequence_set_name,
            organism_id       => $organism_id,
        });
       
        $sequence_set_id = $sequence_set_obj->sequence_set_id();
    
        while (my $seq = $in->next_seq()) {
	
            my $sequence_obj = BAP::DB::Sequence->insert({
                sequence_set_id => $sequence_set_id,
                sequence_name   => $seq->display_id(),
                sequence_string => $seq->seq(),
            });
	
        }
	
        BAP::DB::DBI->dbi_commit();
    }

    unless (defined($sequence_set_id))
    {
        $sequence_set_id = $sequence_set_name_obj->sequence_set_id();
    }

    my @list =($organism_id,$sequence_set_name, $sequence_set_id);
    print join("\t",@list),"\n\n";

    my (
	$glimmer3_model, $glimmer3_pwm,
        $genemark_model,
	$runner_count,
       );

    $glimmer3_model = $cwd."/Sequence/".$locus_tag."_gl3.icm";

    if (-z $glimmer3_model) 
    {
        croak "glimmer3-model: '$glimmer3_model ' is empty! Predict.pm\n";
    }

    $glimmer3_pwm   = $cwd."/Sequence/".$locus_tag."_gl3.motif";

    if (-z $glimmer3_pwm) 
    {
        croak "glimmer3-pwm: '$glimmer3_pwm' is empty! Predict.pm \n";
    }

    my $model_file = undef;
    
    my $idir = IO::Dir->new($cwd."/Sequence");

    while(defined(my $fname = $idir->read))
    {
      if($fname =~ /heu_11_(\d+).mod/)
        {
	  $model_file = $fname;
        }
    }
    $idir->close;


    $genemark_model = $cwd."/Sequence/$model_file";


    if (-z $genemark_model) 
    {
        croak "genemark_model: '$genemark_model' is empty! Predict.pm \n";
    }

    unless (-f $genemark_model)
    {
        croak "genemark_model $genemark_model doesn't exist! Predict.pm \n";
    }
    # wow; were we hard coding this value all this time?
    #$runner_count   = 50;
    $runner_count   = $self->runner_count;

    my $bpg_job_stdout         = $cwd."/".$locus_tag."_bpg_BAP_job_".$sequence_set_id.".txt";
    my $bpg_job_stderr         = $cwd."/".$locus_tag."_bpg_BAP_job_err_".$sequence_set_id.".txt";
    my $bappredictgenes_output = $cwd."/".$locus_tag."_bpg_BAP_screenoutput_".$sequence_set_id.".txt";

    #print qq{\nbap_predict_genes.pl\n};
    # not execing script anymore

    my %params = (
                  'sequence_set_id' => $sequence_set_id,
                  'domain' => 'bacteria',
                  'glimmer3_model' => $glimmer3_model,
                  'glimmer3_pwm' => $glimmer3_pwm,
                  'genemark_model' => $genemark_model,
                  'runner_count' => $runner_count,
                  'job_stdout' => $bpg_job_stdout,
                  'job_stderr' => $bpg_job_stderr
                  );

    if(defined($self->dev)) { $params{dev} = 1; }
    print Dumper(\%params),"\n";

    return %params;
}

sub validate {
    my $self = shift;
    my $params = shift;
    my $cwd = getcwd();
    my %valid = ( );
    %valid = $self->cmdline_to_hash();
    $valid{params} = %{$params};
    # do we need to unlink any existing files?
    DumpFile("$cwd/prediction_valid",\%valid);
#    write_file("$cwd/prediction_valid","step valid\n"); 

    return 1
}

sub is_valid {
    my $self = shift;
    my %params = $self->gather_details();
    my %to_compare = $self->cmdline_to_hash();
    $to_compare{params} = %params;
    use Data::Compare;
    my $cwd = getcwd();
    # check for valid previous run file?
    if(-f "$cwd/prediction_valid") {
        my %validation = %{ LoadFile("$cwd/prediction_valid")} ;
        # compare %validation and %params
        print Dumper(\%to_compare,\%validation),"\n";
        if(Compare(\%to_compare,\%validation) ) {
        $self->status_message("previous valid prediction run exists");
        return 1;
        }
        else {
            $self->status_message("$cwd/prediction_valid exists, but contents do not match params");
            return 0;
        }
    }
    return 0;
}

sub cmdline_to_hash {
    my $self = shift;
    my %cmdline;
    $cmdline{gram_stain} = $self->gram_stain;
    $cmdline{locus_id} = $self->locus_id;
    $cmdline{locus_tag} = $self->locus_tag;
    $cmdline{organism_name} = $self->organism_name;
    $cmdline{project_type} = $self->project_type;
    $cmdline{dev} = $self->dev;
    $cmdline{ncbi_taxonomy_id} = $self->ncbi_taxonomy_id;
    $cmdline{runner_count} = $self->runner_count;
    return %cmdline;
}

1;

# $Id$
