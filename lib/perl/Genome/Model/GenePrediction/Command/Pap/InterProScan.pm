#$Id$

package Genome::Model::GenePrediction::Command::Pap::InterProScan;

use strict;
use warnings;

use Workflow;

use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;

use Compress::Bzip2;
use English;
use File::Temp;
use IO::File;
use IPC::Run;
use File::Slurp;
use Data::UUID;
use Cwd;
use Carp 'confess';


my $debug_err = '';

class Genome::Model::GenePrediction::Command::Pap::InterProScan {
    is => ['Command::V1'],
        has => [
                fasta_file => {  
                               is  => 'SCALAR', 
                               doc => 'fasta file name',
                               is_input => 1,
                           },
                iprscan_output => {
                                   is          => 'SCALAR',
                                   is_optional => 1,
                                   doc         => 'instance of File::Temp pointing to raw iprscan output',
                               },
                iprscan_sorted_output => {
                                          is          => 'SCALAR',
                                          is_optional => 1,
                                          doc         => 'instance of File::Temp pointing to sorted iprscan output',
                                      },
                bio_seq_feature => { 
                                    is          => 'ARRAY',
                                    is_optional => 1,
                                    doc         => 'array of Bio::Seq::Feature',
                                    is_output => 1,
                                },
                report_save_dir => {
                                    is          => 'ARRAY',
                                    is_optional => 1,
                                    doc         => 'directory to save a copy of the raw output to',
                                    is_input => 1,
                                },
				version => {
									is => 'String',
									valid_values => ['4.5', '4.7', '4.8'],
									default => '4.8',
								},

				lsf_queue => { is_param => 1, 
                               default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},},
                lsf_resource => { is_param => 1,
                                  default_value => 'rusage[tmp=100]' },
                locus_tag  => {
                                    is          => 'SCALAR',
                                    doc         => 'locus tag for this genome',
									is_input => 1,
                                },
            ],
};

#operation Genome::Model::GenePrediction::Command::Pap::InterProScan {
#   input        => [ 'fasta_file', 'report_save_dir', 'locus_tag' ],
#  output       => [ 'bio_seq_feature' ],
# lsf_queue    => 'long',
#lsf_resource => 'rusage[tmp=100]',
#   
#};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Run iprscan";
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

	my $iprscan_path = "/gsc/scripts/pkg/bio/iprscan/iprscan-". $self->version ."/bin/iprscan";


    my $fasta_file  = $self->fasta_file();

    my $tmp_fh = File::Temp->new();
    my $tmp_fn = $tmp_fh->filename();

    my $sorted_tmp_fh = File::Temp->new();
    my $sorted_tmp_fn = $sorted_tmp_fh->filename();
    
    $self->iprscan_sorted_output($sorted_tmp_fh);
    $self->iprscan_output($tmp_fh);
    
    my @iprscan_command = (
                           # '/gsc/scripts/pkg/bio/iprscan/iprscan-4.7/bin/iprscan',
						   $iprscan_path,
                           '-cli',
                           '-appl', 'hmmpfam',
                           '-appl', 'hmmtigr',
                           '-goterms',
                           '-verbose',
                           '-iprlookup',
                           '-seqtype',
                           'p',
                           '-format',
                           'raw',
                           '-i',
                           $fasta_file,
                           '-o',
                           $tmp_fn,
                       );

    my ($iprscan_stdout, $iprscan_stderr);
   
    my $rv = IPC::Run::run(
                  \@iprscan_command, 
                  \undef, 
                  '>', 
                  \$iprscan_stdout, 
                  '2>', 
                  \$iprscan_stderr, 
                 );

    $debug_err = $self->report_save_dir."/debug.err.".$$;
    #write_file($self->report_save_dir."/debug.err.".$$,

    write_file($debug_err,
               $iprscan_stdout."\n====\n".$iprscan_stderr);
    unless($rv) {
        die "iprscan failed: $CHILD_ERROR";
    }

    # FIXME:
    # need to check thru the $iprscan_stdout and $iprscan_stderr for
    # SUBMITTED iprscan-<YYYYMMDD>-<HHMMSSSS>
    # and any other error values.

    my ($sort_stderr);
   
    ## Replicate the way Prat sorted the InterProScan
    ## output (sort -k 1,1 -k 2,2 -k 4,4), except
    ## that the field numbers are different here
    ## because we aren't selecting a subset of the
    ## columns via awk first.
    IPC::Run::run(
                  [
                   'sort',
                   '-k',
                   '1,1',
                   '-k',
                   '4,4',
                   '-k',
                   '12,12',
                  ],
                  '<',
                  $tmp_fn,
                  '>',
                  $sorted_tmp_fn,
                  '2>',
                  \$sort_stderr,
                 ) || die "sort failed: $CHILD_ERROR";
                 
    $self->parse_result();

    ## Be Kind, Rewind.  Somebody will surely assume we've done this, 
    ## so let's not surprise them.  
    $tmp_fh->seek(0, SEEK_SET);
    $sorted_tmp_fh->seek(0, SEEK_SET);

    $self->archive_result();

    ## Second verse, same as the first.
    $tmp_fh->seek(0, SEEK_SET);
    $sorted_tmp_fh->seek(0, SEEK_SET);

    return 1;

}

sub parse_result {

    my $self = shift;


    my $iprscan_fh = $self->iprscan_output();
    
    my @features = ( );

    my %least_evalue = ( );
    
    while (my $line = <$iprscan_fh>) {

        chomp $line;
        
        ## iprscan raw format is described here:
        ## ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/README.html#3
        my (
            $protein_name,
            $checksum,
            $length,
            $analysis_method,
            $database_entry,
            $database_member,
            $start,
            $end,
            $evalue,
            $status,
            $run_date,
            $ipr_number,
            $ipr_description,
            $go_description,
        ) = split /\t/, $line;

        ## Some selection criteria of Prat's
        ## refactored to be easier on the eyes
        ## and brain:
        if (
            ($ipr_number ne 'NULL') &&
            (
             ($evalue eq 'NA') ||
             ($evalue <= 0.01)
         )
        ){

            ## More Prat logic...this had me
            ## scratching my head, too...
            if ($evalue eq 'NA') {
                $evalue = 1e10;
            }

            my $feature = Bio::SeqFeature::Generic->new(
                                                        -display_name => $protein_name,
                                                        -primary      => $analysis_method,
                                                        -source_tag   => 'InterPro',
                                                        -start        => $start,
                                                        -end          => $end,
                                                        -score        => $evalue,
                                                    );
            
            $feature->add_tag_value('interpro_analysis',    $analysis_method);
            $feature->add_tag_value('interpro_evalue',      $evalue);
            $feature->add_tag_value('interpro_description', $ipr_description);
            
            my $dblink = Bio::Annotation::DBLink->new(
                                                      -database   => 'InterPro',
                                                      -primary_id => $ipr_number,
                                                  );

            $feature->annotation->add_Annotation('dblink', $dblink);

            ## Prat did thus, so thus do we (with more style).
            ## The <= may be a bug, but the current target is replication
            ## of results, bugs and all.
            if (exists($least_evalue{$protein_name}{$analysis_method}{$ipr_number})) {
                if ($evalue <= $least_evalue{$protein_name}{$analysis_method}{$ipr_number}->score()) {
                    $least_evalue{$protein_name}{$analysis_method}{$ipr_number} = $feature;
                }
            }
            else {
                $least_evalue{$protein_name}{$analysis_method}{$ipr_number} = $feature
            }
            
        }

    }

    foreach my $gene (keys %least_evalue) {

        foreach my $analysis (keys %{$least_evalue{$gene}}) {

            foreach my $ipr (keys %{$least_evalue{$gene}{$analysis}}) {

                push @features, $least_evalue{$gene}{$analysis}{$ipr};
                
            }

        }
        
    }
    
    $self->bio_seq_feature(\@features);
    
}

sub archive_result {

    my $self = shift;
    my $locus_tag = $self->locus_tag();


    my $report_save_dir = $self->report_save_dir();

    if (defined($report_save_dir)) {

        unless (-d $report_save_dir) {
            die "does not exist or is not a directory: '$report_save_dir'";
        }

        my $sorted_output_handle = $self->iprscan_sorted_output();

        my $sorted_target_file = File::Spec->catfile($report_save_dir, 'merged.raw.sorted.bz2');
        
        my $sorted_bz_file = bzopen($sorted_target_file, 'w') or
            die "Cannot open '$sorted_target_file': $bzerrno";

        while (my $line = <$sorted_output_handle>) {
            $sorted_bz_file->bzwrite($line) or die "error writing: $bzerrno"; 
        }

        $sorted_bz_file->bzclose();


        my $output_handle = $self->iprscan_output();

        my $target_file = File::Spec->catfile($report_save_dir, 'merged.raw.bz2');
        
        my $bz_file = bzopen($target_file, 'w') or
            die "Cannot open '$target_file': $bzerrno";
        
        while (my $line = <$output_handle>) {
            $bz_file->bzwrite($line) or die "error writing: $bzerrno"; 
        }
        
        $bz_file->bzclose();

		if ( -s $target_file ) {

			my $back_up_dir = "/gscmnt/gc2103/info/annotation/pap_tmp/iprscan";
		    my $user = $ENV{USER};

			my $fh = IO::File->new($debug_err, "r");
			confess "Could not get handle: $!"unless $fh;

			my $iprscan_dir = "/gsc/var/tmp/iprscan/";

			while (my $line = $fh->getline) {
				chomp $line;
				if ($line =~ m/SUBMITTED iprscan-(\d+)-(\d+)/) {
					$iprscan_dir .= "$1/iprscan-$1-$2";
					confess "Cannot build iprscan_dir" unless( $1 && $2 );

				}
			}
			
			my $cwd = getcwd();
			unless ( $cwd eq $back_up_dir )
			{
				chdir($back_up_dir)
				|| confess "Failed to change to $back_up_dir...\n\n";
			}

			#my $tar_file_name = $params->{locus_tag}. ".tar.bz2";
			my $tar_file_name = $locus_tag. "-". $user. ".tar.bz2";
			print $tar_file_name."\n";
			qx(tar -jcvf $tar_file_name $iprscan_dir);

		}
        
    }

    return 1;
    
}

1;
