package Genome::Model::GenePrediction::Command::Pap::KEGGScan;

use strict;
use warnings;

use PAP;
use Carp 'confess';

use Bio::Annotation::DBLink;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use Compress::Bzip2;
use English;
use File::Basename;
use File::chdir;
use File::Temp;
use IO::File;
use IPC::Run;
use Cwd;
use File::Path qw(remove_tree);
use DateTime;

class Genome::Model::GenePrediction::Command::Pap::KEGGScan {
    is  => 'Command::V1',
    has => [
        fasta_file => { 
            is => 'Path',
            doc => 'fasta file name',            
            is_input => 1,
        },
    ],
    has_optional => [
        report_save_dir   => {
            is => 'Path',
            doc => 'directory to save a copy of the raw output to',
            is_input => 1,
        },
        bio_seq_feature   => { 
            is => 'ARRAY',  
            doc => 'array of Bio::Seq::Feature', 
            is_output => 1,
        },
        keggscan_top_output => {
            is => 'SCALAR',
            doc => 'instance of IO::File pointing to raw KEGGscan output',
        },
        keggscan_full_output => {
            is => 'SCALAR',
            doc => 'instance of IO::File pointing to raw KEGGscan output',
        },
        # These parameters are used when scheduling blast jobs, do NOT affect this lsf requirements
        # for this module are any modules it calls!
        blast_lsf_queue => { 
            is => 'Text',
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        blast_lsf_resource => { 
            is => 'Text',
            default_value => 'select[mem>4096] rusage[mem=4096]',
        }, 
        blast_lsf_max_memory => {
            is => 'Number',
            default => '4096000',
        },
        blast_lsf_job_limit => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of BLAST LSF jobs allowed to run at the same time',
        },
        fasta_chunk_size => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of sequences allowed in a fasta chunk',
        },
        _working_directory => {
            is => 'Path',
            doc => 'analysis program working directory',
        },
        version => {
                is => 'String',
                valid_values => ['50', '52', '56'],
                default => '56',
		},
		tar	=> {
				is	=> 'Boolean',
				is_input	=> 1,
				default	=> 1,
				doc	=> 'If set, tar/bz2 PAP_keggscan_* directory',
        },
    ],
    # These parameters tell workflow the requirements needed for this module 
    has_param => [
        lsf_resource => {
            default => "-R 'select[mem=8192,type==LINUX64] rusage[mem=8192,tmp=1024]' -M 8192000",
        },
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub help_brief {
    "Run KEGGscan";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;

    my $keggscan_faa_path = "/gscmnt/temp212/info/annotation/KEGG/Version_". $self->version . "/genes.v" . $self->version .".faa";

    # Status messages are not displayed by default...
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    # FIXME This directory needs to be a param with a default value
    my ($kegg_stdout, $kegg_stderr);    
    $self->_working_directory(
        File::Temp::tempdir(
            'PAP_keggscan_XXXXXXXX',
            DIR => '/gscmnt/temp212/info/annotation/PAP_tmp',
            CLEANUP => 0,
        )
    );
    chmod(0777, $self->_working_directory);
    my $fasta_file = $self->fasta_file();

    ## We're about to screw with the current working directory.
    ## Thusly, we must fixup the fasta_file property if it 
    ## contains a relative path. 
    #TODO bdericks: Is this necessary?
    unless ($fasta_file =~ /^\//) {
        $fasta_file = join('/', $CWD, $fasta_file);
        $self->fasta_file($fasta_file);
    }
    
    # FIXME The subject fasta path needs to be a param on this command object with a default value
    my ($fasta_name) = basename($self->fasta_file);
    my $kegg_command = Genome::Model::GenePrediction::Command::Pap::KEGGScan::RunKeggScan->create(
        species_name => "default", 
        query_fasta_path => $self->fasta_file,
        subject_fasta_path => $keggscan_faa_path,
        output_directory => $self->_working_directory . "/KS-OUTPUT." . $fasta_name,
        blast_lsf_queue => $self->blast_lsf_queue,
        blast_lsf_resource => $self->blast_lsf_resource,
        blast_lsf_max_memory => $self->blast_lsf_max_memory,
        fasta_chunk_size => $self->fasta_chunk_size,
        blast_lsf_job_limit => $self->blast_lsf_job_limit,
    );

    my $rv = $kegg_command->execute;
    unless ($rv) {
        $self->error_message("Problem running keggscan!");
        confess;
    }
        
    $self->parse_result();

    my $top_output_fh = $self->keggscan_top_output();
    my $full_output_fh = $self->keggscan_full_output();
    
    ## Be Kind, Rewind.  Somebody will surely assume we've done this,
    ## so let's not surprise them.
    # TODO bdericks: Is this necessary?
    $top_output_fh->seek(0, SEEK_SET);
    $full_output_fh->seek(0, SEEK_SET);
    
    $self->archive_result();

    ## Second verse, same as the first.
    $top_output_fh->seek(0, SEEK_SET);
    $full_output_fh->seek(0, SEEK_SET);
    
    return 1;
}

sub parse_result {

    my $self = shift;
 
    my $top_output_fh  = IO::File->new();
    my $full_output_fh = IO::File->new();
    
    my $top_output_fn = join('.', 'KS-OUTPUT', File::Basename::basename($self->fasta_file()));

    $top_output_fn = join('/', $self->_working_directory(), $top_output_fn, 'REPORT-top.ks');    
    
    $top_output_fh->open("$top_output_fn") or die "Can't open '$top_output_fn': $OS_ERROR";

    $self->keggscan_top_output($top_output_fh);

    my $full_output_fn = join('.', 'KS-OUTPUT', File::Basename::basename($self->fasta_file()));

    $full_output_fn = join('/', $self->_working_directory(), $full_output_fn, 'REPORT-full.ks');    
    
    $full_output_fh->open("$full_output_fn") or die "Can't open '$full_output_fn': $OS_ERROR";

    $self->keggscan_full_output($full_output_fh);
    
    my @features = ( );
    
    LINE: while (my $line = <$top_output_fh>) {

        chomp $line;

        my @fields = split /\t/, $line;
        
        my (
            $ec_number,
            $gene_name,
            $hit_name,
            $e_value,
            $description,
            $orthology
           ) = @fields[1,2,3,4,6,8];

        ## Prat filtered out lines with e-value >= 0.01 when 
        ## generating AceDB .ace files.
        if ($e_value >= 0.01) { next LINE; }
        
        ## The values in the third column should be in this form:
        ## gene_name (N letters; record M).
        ($gene_name) = split /\s+/, $gene_name; 
		unless ($main::locus_name) {
			$self->status_message("gene_name: ". $gene_name);
			($main::locus_name) = split(/_/, $gene_name);
			$self->status_message("locus_name: ". $main::locus_name);
		}

        ## Some descriptions have EC numbers embedded in them.
        ## Prat's removed them.
        $description =~ s/\(EC .+\)//;
        
        my $feature = new Bio::SeqFeature::Generic(-display_name => $gene_name);

        $feature->add_tag_value('kegg_evalue', $e_value);
        $feature->add_tag_value('kegg_description', $description);
        
        my $gene_dblink = Bio::Annotation::DBLink->new(
                                                       -database   => 'KEGG',
                                                       -primary_id => $hit_name,
                                                   );
        
        $feature->annotation->add_Annotation('dblink', $gene_dblink);

        ## Sometimes there is no orthology identifier (value is literally 'none').
        ## It is not unforseeable that it might also be blank/undefined.  
        if (defined($orthology) && ($orthology ne 'none')) {
        
            my $orthology_dblink = Bio::Annotation::DBLink->new(
                                                                -database   => 'KEGG',
                                                                -primary_id => $orthology,
                                                            );
            
            $feature->annotation->add_Annotation('dblink', $orthology_dblink);
                
        }
        
        push @features, $feature;
        
    }
    
    $self->bio_seq_feature( \@features );

    return;
    
}

sub archive_result {

    my $self = shift;
    
    
    my $report_save_dir = $self->report_save_dir();
    
	my $tar_bz2 = $self->tar;
    
    if (defined($report_save_dir)) {
        
        unless (-d $report_save_dir) {
            die "does not exist or is not a directory: '$report_save_dir'";
        }
        
        my $top_output_handle  = $self->keggscan_top_output();
        my $full_output_handle = $self->keggscan_full_output();
            
        my $top_target_file  = File::Spec->catfile($report_save_dir, 'REPORT-top.ks.bz2');
        my $full_target_file = File::Spec->catfile($report_save_dir, 'REPORT-full.ks.bz2');
        
        my $top_bz_file = bzopen($top_target_file, 'w') or
            die "Cannot open '$top_target_file': $bzerrno";
        
        while (my $line = <$top_output_handle>) {
            $top_bz_file->bzwrite($line) or die "error writing: $bzerrno";
        }
        
        $top_bz_file->bzclose();

        
        my $full_bz_file = bzopen($full_target_file, 'w') or
            die "Cannot open '$full_target_file': $bzerrno";
        
        while (my $line = <$full_output_handle>) {
            $full_bz_file->bzwrite($line) or die "error writing: $bzerrno";
        }
        
        $full_bz_file->bzclose();
        
    }

	## We will come here to tar ball PAP_keggscan_* directory
	if (defined($tar_bz2)) {
		my $parent_dir =  $self->_working_directory."/..";
		chdir($parent_dir) || confess "Unable to change directory to $parent_dir: $?";
		$self->status_message("Cwd: ". getcwd);

		my $user = $ENV{USER};

		my $now = DateTime->now(time_zone	=> 'America/Chicago');
		my $date_code	= $now->ymd('_');

#		my $tar_file_name = File::Basename::basename($self->fasta_file()). "-". $user. "-". $date_code. ".tar.bz2";
		my $tar_file_name = $main::locus_name. "-". $user. "-". $date_code. ".tar.bz2";
		my $tar_dir = $self->_working_directory;

		my $tar_bz_cmd = "tar -C $tar_dir -jcvf $tar_file_name .";
		$self->status_message("Cmd: ". $tar_bz_cmd);

		system($tar_bz_cmd) == 0 
			|| confess "system $tar_bz_cmd failed: $?";

		$self->status_message("Removing PAP_tmp dir: ". $tar_dir);
		remove_tree($tar_dir) || $self->status_message("Error removing $tar_dir: ". $?);	

	}

    return 1;
    
}

1;
