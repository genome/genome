package Genome::Model::Tools::PooledBac::Run;

use strict;
use warnings;

use DateTime;
use Genome;
use Bio::SeqIO;
use PP::LSF;
use Data::Dumper;
use Cwd;

class Genome::Model::Tools::PooledBac::Run {
    is => 'Command',
    has => 
    [        
        project_dir =>
        {
            type => 'String',
            is_optional => 1,
            doc => "location of the finished pooled BAC projects"        
        },
        pooled_bac_dir =>
        {
            type => 'String',
            is_optional => 1,
            doc => "location of the input pooled BAC assembly"        
        },
        ref_seq_file =>
        {
            type => 'String',
            is_optional => 1,
            doc => "location of the reference sequence"        
        },
        ace_file_name =>
        {
            type => 'String',
            is_optional => 1,
            doc => "location of the finished pooled BAC projects"        
        },
        phd_ball_name =>
        {
            type => 'String',
            is_optional => 1,
            doc => "location of the finished pooled BAC projects"        
        },    
        percent_overlap => 
        {
            type => 'String',
            is_optional => 1,
            doc => "this is the percent overlap, default is 50%",
        },
        percent_identity =>
        {
            type => 'String',
            is_optional => 1,
            doc => "this is the percent identity, default is 85%",
        },
        params_file => 
        {
            type => 'String',
            is_optional => 1,
            doc => "Use this option toe specify the path to a params file, a file containing preset options to run the pipeline with",
        }, 

    ]
};

sub help_brief {
    "Run Pooled BAC Pipeline"
}

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Assemble Pooled BAC Reads
EOS
}

sub save_params
{
    my ($self, $file_name) = @_;
    my $params = $self->create_params_hash;    

    my $fh = IO::File->new(">$file_name");
    
    print $fh Dumper($params);

    return;
}

sub create_params_hash
{
    my ($self, $file_name) = @_;
    my $params = 
    { 
        project_dir => $self->project_dir,
        pooled_bac_dir => $self->pooled_bac_dir,
        ace_file_name => $self->ace_file_name,
        ref_seq_file => $self->ref_seq_file,
        phd_ball => $self->phd_ball_name,
        percent_overlap => $self->percent_overlap,
        percent_identity => $self->percent_identity,
    };

    return $params;
}

sub params_are_equal
{
    my ($self,$params1, $params2) = @_;
    foreach my $key (keys %{$params1})
    {
        return 0 if (defined $params1->{$key} && defined $params2->{$key} &&($params1->{$key} ne $params2->{$key}));
    }

    return 1;
}

sub get_params
{
    my ($self,$param_file) = @_;
    my $VAR1;
    my $params = eval `cat $param_file`;
    return $params;    
}
############################################################
sub execute { 
    my ($self) = @_;
    unless (`uname -m` =~ /64/) {
        $self->error_message('Pooled bac pipeline must be run from a 64-bit architecture');
        return;
    }
        
    my $orig_dir = cwd();

    my $params = {};
    $params = $self->get_params($self->params_file) if(defined $self->params_file && -e $self->params_file);
    my $project_dir = $self->project_dir || $self->project_dir($params->{project_dir});
    my $pooled_bac_dir = $self->pooled_bac_dir || $self->pooled_bac_dir($params->{pooled_bac_dir});
    my $ace_file_name = $self->ace_file_name || $self->ace_file_name($params->{ace_file_name} ||'Pcap.454Contigs.ace.1');
    my $ref_seq_file = $self->ref_seq_file || $self->ref_seq_file($params->{ref_seq_file});
    my $phd_ball = $self->phd_ball_name || $self->phd_ball_name($params->{phd_ball});
    my $percent_overlap = $self->percent_overlap || $self->percent_overlap($params->{percent_overlap});
    my $percent_identity = $self->percent_identity || $self->percent_identity($params->{percent_identity});


    $self->error_message("The pipeline needs for the project_dir to be specified in either the params file or on the command line in order to run.\n") and return if(!defined $project_dir);
    $self->error_message("The pipeline needs for the pooled_bac_dir to be specified in either the params file or on the command line in order to run.\n") and return if(!defined $pooled_bac_dir);
    $self->error_message("The pipeline needs for the ref_seq_file to be specified in either the params file or on the command line in order to run.\n") and return if(!defined $ref_seq_file);
        
    $self->error_message("Error creating directory $project_dir") and die unless Genome::Sys->create_directory($project_dir);
    my $dt = DateTime->now(time_zone => "America/Chicago");
    if($self->params_file && -e $self->params_file)
    {
        my $current_params = $self->create_params_hash;
        $self->save_params("$project_dir/params.".$dt->strftime("%y%m%d.%I%M\n")) unless ($self->params_are_equal($params,$current_params));
    }
    else
    {
        $self->save_params("$project_dir/params.".$dt->strftime("%y%m%d.%I%M\n"));
    }

    $self->error_message("Error running run-blast")  and die unless
    Genome::Model::Tools::PooledBac::RunBlast->execute(ref_seq_file=>$ref_seq_file, pooled_bac_dir=>$pooled_bac_dir,ace_file_name => $ace_file_name, project_dir => $project_dir);

    #run gmt map-contigs-to-assembly
    my %map_contigs_params = (
	pooled_bac_dir=>$pooled_bac_dir,
	ace_file_name => $ace_file_name,
	project_dir => $project_dir,
	percent_overlap => $percent_overlap,
	percent_identity => $percent_identity,
	);
    #following can be undefined if not $self->pecent_overlap and $self->percent_identity and not specified in params file
    #delete from params so MapContigsToAssembly can use it's own default param
    delete $map_contigs_params{'percent_overlap'} unless $percent_overlap;
    delete $map_contigs_params{'percent_identity'} unless $percent_identity;
    $self->error_message("Error running map-contigs-to-assembly")  and die unless
    Genome::Model::Tools::PooledBac::MapContigsToAssembly->execute( %map_contigs_params );

    $self->error_message("Error running add-linking-contigs")  and die unless
    Genome::Model::Tools::PooledBac::AddLinkingContigs->execute( project_dir => $project_dir);

    $self->error_message("Error generating reports")  and die unless
    Genome::Model::Tools::PooledBac::GenerateReports->execute( project_dir => $project_dir);

    $self->error_message("Error creating project directories")  and die unless
    Genome::Model::Tools::PooledBac::CreateProjectDirectoriesNew->execute(pooled_bac_dir=>$pooled_bac_dir,ace_file_name => $ace_file_name,phd_file_name_or_dir => $phd_ball, project_dir => $project_dir);

    $self->error_message("Error generating post assembly reports")  and die unless
    Genome::Model::Tools::PooledBac::GeneratePostAssemblyReports->execute( project_dir => $project_dir);

#    $self->error_message("Error updating seqmgr") unless
#    Genome::Model::Tools::Lims::UpdateSeqMgr->execute(project_dir => $project_dir);

    chdir( $orig_dir );

    return 1;
}


1;
