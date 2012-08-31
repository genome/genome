#$Id$

package PAP::Command::BlastP;

use strict;
use warnings;


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


class PAP::Command::BlastP {
    is  => 'PAP::Command',
    has => [
        fasta_file => { 
            is => 'SCALAR', 
            doc => 'fasta file name',
            is_input => 1,
        },
    ],
    has_optional => [
        blast_report => {
            is => 'SCALAR',
            doc => 'instance of File::Temp pointing to raw blast output',
        },
        bio_seq_feature => { 
            is => 'ARRAY',  
            doc => 'array of Bio::Seq::Feature' ,
            is_output => 1,
        },
        report_save_dir => {
            is => 'DirectoryPath',
            doc => 'directory to save a copy of the blast report to',
            is_input => 1,
        },
        query_names => {
            is => 'ARRAY',
            doc => 'array of sequence (query) names seen in the results',
        },
    ],
    has_param => [
        lsf_resource => { default_value => '-q long rusage[tmp=100]', },
    ],
};

sub execute {
    my $self = shift;
 
    my ($blastp_out, $blastp_err);

    ##FIXME:  This should not be hardcoded.  At least not here.
    my $bacterial_nr = '/gscmnt/ams1102/info/annotation/blastdb/gsc_bacterial/bacterial_nr/bacterial_nr';

    my $fasta_file = $self->fasta_file();
    
    my $temp_fh    = File::Temp->new();
    my $temp_fn    = $temp_fh->filename();

    $self->blast_report($temp_fh);

    ## If 'blastp' invokes anything but WU-BLAST, stuff will probably
    ## go seriously foul in archive_result below
    my @blastp_command = (
                          'blastp',
                          $bacterial_nr,
                          $fasta_file,
                          '-o',
                          $temp_fn,
                          'E=1e-10',
                          'V=1',
                          'B=50',
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
 
    
    ## According to the docs for Bio::Root::IO,
    ## -noclose should prevent the filehandle
    ## from being closed when $searchio gets
    ## garbage collected.  
    my $searchio = Bio::SearchIO->new(
                                      -format  => 'blast',
                                      -fh      => $self->blast_report(),
                                      -noclose => 1, 
                                  );

    my @features    = ( );
    my @query_names = ( );
    
    ## There should be one result per query sequence.
    ## The protein category is determined by the top (first) hit.
    ## The dblink is also set from the top hit, unless it's 'hypothetical'.
  RESULT: while (my $result = $searchio->next_result()) {

        my $query_name = $result->query_name();

        push @query_names, $query_name;
        
        my $feature = Bio::SeqFeature::Generic->new(-display_name => $query_name);

        my $protein_category = 'Predicted Protein';
        
        my $hit = $result->next_hit();

        if (defined($hit)) {
            
            my $hit_name        = $hit->name();
            my $hit_accession   = $hit->accession();
            my $hit_description = $hit->description();
            
            my $hsp = $hit->next_hsp();
            
            unless (defined($hsp)) { next RESULT; }
            
            if (defined($hsp)) {
                
                my $score       = $hsp->score;
                my $evalue      = $hsp->evalue;
                my $pctid       = sprintf("%.1f", $hsp->percent_identity());

                ## This is equivalent to the coverage reported by nrblast
                my $pctcov      = sprintf("%.1f", (($hsp->num_identical() / $hit->logical_length('query')) * 100));
                
                my $qstart      = $hsp->start('query');
                my $qend        = $hsp->end('query');
                my $sstart      = $hsp->start('subject');
                my $send        = $hsp->end('subject');
                
                if (($pctid >= 80) && ($pctcov >= 80)) {

                    my @hit_descriptions = split /\s+\>/, $hit_description;

                    my ($first_description) = @hit_descriptions;

                    my $all_magic_words = 1;
                    
                    foreach my $description (@hit_descriptions) {
                        unless ($hit_description =~ /fragment|homolog|hypothetical|like|predicted|probable|putative|related|similar|synthetic|unknown|unnamed/) {
                            $all_magic_words = 0;
                        }
                            
                    }
                    
                    if ($all_magic_words) {
                        
                        $protein_category = 'Conserved Hypothetical Protein';
                        
                    }
                    else {
                        
                        my $analogue = $first_description;
                        
                        $analogue =~ s/\[.*\]//x;
                        
                        $protein_category = "Hypothetical Protein similar to $analogue";
                        
                    }
                    
                }
                else {
                    $protein_category = 'Hypothetical Protein';
                }
                
                if($hit_description !~ /hypothetical/i) {
                    
                    $hit_description =~ s/>.*//;
                    
                    $feature->add_tag_value('blastp_bit_score', $score);
                    $feature->add_tag_value('blastp_evalue', $evalue);
                    $feature->add_tag_value('blastp_percent_identical', $pctid);
                    $feature->add_tag_value('blastp_query_start', $qstart);
                    $feature->add_tag_value('blastp_query_end', $qend);
                    $feature->add_tag_value('blastp_subject_start', $sstart);
                    $feature->add_tag_value('blastp_subject_end', $send);
                    $feature->add_tag_value('blastp_hit_name', $hit_name);
                    $feature->add_tag_value('blastp_hit_description', $hit_description);
                    
                    my $dblink = Bio::Annotation::DBLink->new(
                                                              -database   => 'GenBank',
                                                              -primary_id => $hit_accession,
                                                          );
                    
                    $feature->annotation->add_Annotation('dblink', $dblink);
                    
                }
                
            }
            
        }
        
        $feature->add_tag_value('blastp_category', => $protein_category);
        
        push @features, $feature;
        
    }
    
    $self->bio_seq_feature(\@features);
    $self->query_names(\@query_names);
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
        
            my $target_file = File::Spec->catfile($report_save_dir, "$query_name.blastp.bz2");
            
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
