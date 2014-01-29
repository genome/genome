package Genome::Model::Tools::Lims::UpdateSeqMgr;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;

class Genome::Model::Tools::Lims::UpdateSeqMgr {
    is => 'Command',
    has => [        
        project_dir => {
            type => 'String',
            is_optional => 0,
            doc => "output dir for separate pooled bac projects"        
        },
    ],
};

sub help_brief {
    'Tool to copy assembly files into seqmgr linke pooled bac project directories';
}

sub help_synopsis {
    return <<"EOS"
gmt lims update-seq-mgr --project-dir /gscmnt/111/New_Pipeline_Pooled_BAC_Fosmid_Testing/Human_Pool_10_101119
EOS
}
sub help_detail {
    return <<"EOS"
This tool copies assembly files created by pooled-bac pipline into individual seqmgr
linked GSC project directories for finishing and manual sequencing improvement
EOS
}


sub execute { 
    my ($self) = @_;

    if ( not App::Init->initialized ) {
        App::DB->db_access_level('rw');
        App->init;
    }

    $self->debug_message ("Updating Seqmgr Projects ..");
    my $project_dir = $self->project_dir;

    chdir($project_dir);
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -file => 'ref_seq.fasta');
    $self->error_message("Erroring opening ref_seq_fasta") and die unless defined $seqio;

    while (my $seq = $seqio->next_seq)
    {    
        my $clone_name = $seq->display_id;

        my $p=GSC::Project->get(name => $clone_name);
        $self->error_message("Error retrieving project for $clone_name")  and die unless defined $p;
    
        unless ((-e "$project_dir/$clone_name") && (-d "$project_dir/$clone_name"))
        {
            $self->warning_message("Directory for $clone_name at $project_dir/$clone_name does not exist!");
            next;
        }
	#skip if project is not linked in seqmgr
	unless ( -l '/gscuser/seqmgr/'.$clone_name ) {
	    $self->warning_message( $clone_name." is not linked in seqmgr .. skipping" );
	    next;
	}

        chdir($project_dir."/$clone_name");
        next if (-e 'core'||!(-e 'edit_dir'));
        $p->set_project_status('pooled_bac_done'); 
        
        my $seqmgr_link = $p->seqmgr_link;
        $self->debug_message ("Updating $clone_name ..");

	unless (-d "$seqmgr_link/edit_dir" ) {
	    Genome::Sys->create_directory("$seqmgr_link/edit_dir");
	}
        foreach my $file (glob('edit_dir/*'))
        {
            my $new_file_name = File::Basename::basename($file);
            if ( $new_file_name =~ /$clone_name\.ace$/ ) {
                $new_file_name =~ s/\.ace$/\.fasta\.screen\.ace/;
            }
            system "/bin/cp -rfP $file $seqmgr_link/edit_dir/$new_file_name";
        }

	unless (-d "$seqmgr_link/phd_dir" ) {
	    Genome::Sys->create_directory("$seqmgr_link/phd_dir");
	}
        foreach my $phd_file (glob('phd_dir/*'))
        {        
            system "/bin/cp -rfP $phd_file $seqmgr_link/phd_dir/.";
        }

	unless ( -d "$seqmgr_link/chromat_dir" ) {
	    Genome::Sys->create_directory("$seqmgr_link/chromat_dir");
	}
        foreach my $chromat_file (glob('chromat_dir/*'))
        {
            system "/bin/cp -rfP $chromat_file $seqmgr_link/chromat_dir/.";
        }

	unless ( -d "$seqmgr_link/phdball_dir" ) {
	    Genome::Sys->create_directory("$seqmgr_link/phdball_dir");
	}
        foreach my $phdball_file (glob('phdball_dir/*'))
        {
            system "/bin/cp -rfP $phdball_file $seqmgr_link/phdball_dir/.";
        }

	unless ( -d "$seqmgr_link/sff_dir" ) {
	    Genome::Sys->create_directory("$seqmgr_link/sff_dir");
	}
        foreach my $sff_file (glob('sff_dir/*'))
        {
            system "/bin/cp -rfP $sff_file $seqmgr_link/sff_dir/.";
        }
    }
    return 1;
}

1;
