#$Id$

package Genome::Model::GenePrediction::Command::Pap::RRnaScreen;

use strict;
use warnings;

use Workflow;

use Bio::Annotation::DBLink;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

use Compress::Bzip2;
use English;
use File::Basename;
use File::Spec;
use File::Temp;
use IPC::Run;


class Genome::Model::GenePrediction::Command::Pap::RRnaScreen {
    is  => ['Command::V1'],
    has => [
        fasta_file      => { 
                            is  => 'SCALAR', 
                            doc => 'fasta file name',
                           },
        blast_db => {
            is => 'SCALAR',
            doc => 'blast database for rRNA',
            is_optional => 1,
            default => "/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb",
        },
        blast_report => {
                         is          => 'SCALAR',
                         is_optional => 1,
                         doc         => 'instance of File::Temp pointing to raw blast output'
                        },
        bio_seq_feature => { 
                            is          => 'ARRAY',  
                            is_optional => 1,
                            doc         => 'array of Bio::Seq::Feature' 
                           },
        report_save_dir => {
                            is          => 'SCALAR',
                            is_optional => 1,
                            doc         => 'directory to save a copy of the blast report to',
                           },
        dead_genes => {
                        is          => 'ARRAY',
                        is_optional => 1,
                        doc         => 'array of sequence (query) names seen in the results, to mark as \'dead\'',
                       },
    ],
};

1;

__END__
operation Genome::Model::GenePrediction::Command::Pap::RRnaScreen {
    input        => [ 'fasta_file', 'report_save_dir' ],
    output       => [ 'bio_seq_feature'],
    lsf_queue    => 'long',
    lsf_resource => 'rusage[tmp=100]'
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Run blastp";
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

 
    my ($blastp_out, $blastp_err);

    ##FIXME:  This should not be hardcoded.  At least not here.
    my $rrna_db = $self->blast_db();

    my $fasta_file = $self->fasta_file();
    
    my $temp_fh    = File::Temp->new();
    my $temp_fn    = $temp_fh->filename();

    $self->blast_report($temp_fh);

    ## If 'blastp' invokes anything but WU-BLAST, stuff will probably
    ## go seriously foul in archive_result below
    my @blastp_command = (
                          'blastn',
                          $rrna_db,
                          $fasta_file,
                          '-o',
                          $temp_fn,
                          'V=1',
                          'B=1',
                         );

    IPC::Run::run(
                  \@blastp_command,
                  \undef,
                  '>',
                  \$blastp_out,
                  '2>',
                  \$blastp_err,
              ) || die "blastp failed: $CHILD_ERROR";
    
    $self->parse_result();

    ## Be Kind, Rewind.  Somebody will surely assume we've done this,
    ## so let's not surprise them.
    $temp_fh->seek(0, SEEK_SET);

    $self->archive_result();

    ## Second verse, same as the first.
    $temp_fh->seek(0, SEEK_SET);
    
    return 1;

}

sub parse_result {
    
    my $self = shift;
 
    my $fraction = 0.7;
    my $percentage = 70;
    ## According to the docs for Bio::Root::IO,
    ## -noclose should prevent the filehandle
    ## from being closed when $searchio gets
    ## garbage collected.  
    my $searchio = Bio::SearchIO->new(
                                      -format  => 'blast',
                                      -file      => $self->blast_report(),
                                      -noclose => 1, 
                                  );

    my @features    = ( );
    my @query_names = ( );
    
RESULT: while(my $result = $searchio->next_result())
    {
HIT:    my $hit = $result->next_hit;
        next RESULT unless defined($hit);
        my $hsp = $hit->next_hsp; 
        next HIT unless defined($hsp);
        if((int($result->query_length * $fraction) <= $hsp->length('total')) &&
           ($hsp->percent_identity >= $percentage))
        {
            print "========== ", $result->query_name," ==========\n";
            print join( "\t", $hsp->evalue, $hsp->bits, $hsp->percent_identity,
                              $hsp->start('query'), $hsp->end('query'),
                              $hsp->start('subject'), $hsp->end('subject'),
                              $hit->name ), "\n";
            push(@query_names,$result->query_name);
        }
    }
    
    $self->bio_seq_feature(\@features);
    $self->dead_genes(\@query_names);
}

sub archive_result {

    my $self = shift;


    my $report_save_dir = $self->report_save_dir();

    if (defined($report_save_dir)) {

        unless (-d $report_save_dir) {
            die "does not exist or is not a directory: '$report_save_dir'";
        }

        my $report_handle = $self->blast_report();

        ## This is an allegedly clever hack to split the
        ## blast report into multiple reports and end up with
        ## one per file, with the file named as the query
        ## sequence.  This relies on the assumption that
        ## we're dealing with WU-BLAST output, which separates
        ## reports by CTRL-L.

        my $bz_file;
        
        QUERY: foreach my $query_name (@{$self->query_names()}) {
        
            my $target_file = File::Spec->catfile($report_save_dir, "$query_name.bz2");
            
            $bz_file = bzopen($target_file, 'w') or
                die "Cannot open '$target_file': $bzerrno";
            
            while (my $line = <$report_handle>) {

                if ($line =~ qr/^\cL$/) {
                    $bz_file->bzclose();
                    next QUERY;
                }
                
                $bz_file->bzwrite($line) or die "error writing: $bzerrno";
                
            }
            
        }

        ## There will not be a CTRL-L at the end of the last report.
        $bz_file->bzclose();
        
    }

    return 1;
    
}

1;
