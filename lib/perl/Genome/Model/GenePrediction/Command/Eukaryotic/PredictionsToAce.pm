
package Genome::Model::GenePrediction::Command::Eukaryotic::PredictionsToAce;

use strict;
use warnings;
use Genome;
use Carp 'confess';
use File::Temp; 
use Sort::Naturally qw/ ncmp nsort /;

use Bio::SeqIO;

class Genome::Model::GenePrediction::Command::Eukaryotic::PredictionsToAce {
    is => 'Genome::Command::Base',
    has_optional => [
        ace_file => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Path for output ace file',
        },
        protein_coding_only => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'If set, only genes that produce proteins are placed in ace file',
        },
        rna_only => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
            doc => 'If set, only rna genes are placed in the ace file',
        },
        build_id => {
            is => 'Number',
            is_input => 1,
            doc => 'Build ID of build that contains predictions',
        },
        prediction_directory => {
            is => 'DirectoryPath',
            is_input => 1,
            is_output => 1,
            doc => 'Path to directory containing predictions',
        },
        sequence_file => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Path to fasta file containing sequence on which prediction was run',
        },
    ],
};

sub help_brief {
    return "Generates an ace file from predictions creates by a eukaryotic gene prediction build";
}

sub help_synopsis {
    return "Generates an ace file from predictions creates by a eukaryotic gene prediction build";
}

sub help_detail {
    return "Generates an ace file from predictions creates by a eukaryotic gene prediction build";
}
    
sub execute {
    my $self = shift;

    # Either use provided directory or get a build and use it's data directory. Expect to find coding gene,
    # rna, protein, etc prediction files here
    my $prediction_directory = $self->prediction_directory;
    my $sequence_file = $self->sequence_file;
    unless (defined $prediction_directory and defined $sequence_file) {
        if (defined $self->build_id) {
            my $build = Genome::Model::Build::GenePrediction::Command::Eukaryotic->get($self->build_id);
            unless ($build) {
                $self->error_message("Could not get eukaryotic gene prediction with build ID " . $self->build_id);
                confess $self->error_message;
            }
            $prediction_directory = $build->data_directory;
            $sequence_file = $build->model->assembly_contigs_file;
        }
        else {
            $self->error_message("Must provide either a prediction directory or a build ID!");
            confess $self->error_message;
        }
    }

    confess "No sequence file found at $sequence_file!" unless -e $sequence_file;
    confess "No directory found at $prediction_directory!" unless -d $prediction_directory;

    if ($self->protein_coding_only and $self->rna_only) {
        confess "Both protein only and rna only are set. Only set one (or none)!";
    }

    # Now either use the supplied ace file or create a temp one in the predictions directory for output
    my $ace_fh;
    if (defined $self->ace_file) {
        if (-e $self->ace_file) {
            $self->warning_message("Removing existing output ace file at " . $self->ace_file);
            unlink $self->ace_file;
        }
        $ace_fh = IO::File->new(">" . $self->ace_file);
        confess "Could not get handle for " . $self->ace_file unless $ace_fh;
    }
    else {
        $ace_fh = File::Temp->new(
            TEMPLATE => 'predictions_XXXXXX',
            SUFFIX => '.ace',
            UNLINK => 0,
            CLEANUP => 0,
            DIR => $prediction_directory,
        );
        confess "Could not get handle for temp ace file in $prediction_directory" unless $ace_fh;
        $self->ace_file($ace_fh->filename);
    }
     
    $self->status_message("Generating predictions ace file at " . $self->ace_file . " using predictions in " . 
        $prediction_directory . " and sequence in $sequence_file");

    # Pre-fetching all genes now so only one file read is necessary
    my @coding_genes = Genome::Prediction::CodingGene->get(
        directory => $prediction_directory,
    ) unless $self->rna_only;
    my @rna_genes = Genome::Prediction::RNAGene->get(
        directory => $prediction_directory,
    ) unless $self->protein_coding_only;

    # Get list of sequences
    my @sequences = $self->_get_sequences_from_file($sequence_file);
    for my $sequence (nsort @sequences) { 
        my @seq_coding_genes;
        @seq_coding_genes = Genome::Prediction::CodingGene->get(
            directory => $prediction_directory,
            sequence_name => $sequence,
        ) unless $self->rna_only;

        my @seq_rna_genes;
        @seq_rna_genes = Genome::Prediction::RNAGene->get(
            directory => $prediction_directory,
            sequence_name => $sequence,
        ) unless $self->protein_coding_only;

        for my $gene (sort { ncmp($a->gene_name, $b->gene_name) } @seq_coding_genes) {
            my $gene_name = $gene->gene_name;
            my $start = $gene->start;
            my $end = $gene->end;
            ($start, $end) = ($end, $start) if $gene->strand eq '-1';
            my $source = $gene->source;
            my $strand = $gene->strand;

            $ace_fh->print("Sequence $sequence\n");
            $ace_fh->print("Subsequence $gene_name $start $end\n\n");

            my ($transcript) = $gene->transcript;
            my @exons = $transcript->exons;
            @exons = sort { $a->start <=> $b->start } @exons;
            @exons = reverse @exons if $transcript->strand eq '-1';

            my $spliced_length = 0;
            for my $exon (@exons) {
                $spliced_length += abs($exon->end - $exon->start) + 1;
            }

            # Coding gene names follow pattern <sequence_name>.<method>.<id>... need to extract method.
            my ($method) = $gene_name =~ /\Q$sequence\E.(\w+)\..*/;

            $ace_fh->print("Sequence : $gene_name\n");
            $ace_fh->print("Source $sequence\n");
            $ace_fh->print("Method $method\n");
            $ace_fh->print("CDS\t1 $spliced_length\n");
            $ace_fh->print("CDS_predicted_by $source\n");

            if ($gene->missing_start) {
                my $frame = $exons[0]->five_prime_overhang;
                $frame++;  # Predictors use frame 0 - 2, ace requires frame 1 - 3
                $ace_fh->print("Start_not_found $frame\n");
            }
            if ($gene->missing_stop) {
                $ace_fh->print("End_not_found\n");
            }

            my $transcript_start = $transcript->start;
            my $transcript_end = $transcript->end;

            for my $exon (@exons) {
                my $exon_start = $exon->start;
                my $exon_end = $exon->end;

                if ($exon_start > $exon_end) {
                    ($exon_start, $exon_end) = ($exon_end, $exon_start);
                }

                my ($exon_ace_start, $exon_ace_end);
            
                if ($gene->strand eq '+1') {    
                    $exon_ace_start = $exon_start - $transcript_start + 1;
                    $exon_ace_end = $exon_end - $transcript_start + 1;
                }
                elsif ($gene->strand eq '-1') {            
                    $exon_ace_start = $transcript_end - $exon_end + 1;
                    $exon_ace_end = $transcript_end - $exon_start + 1;
                }
                else {
                    die "Bad strand for coding gene " . $gene->gene_name . ": " . $gene->strand;
                }

                $ace_fh->print("Source_Exons $exon_ace_start $exon_ace_end\n");
            }
            $ace_fh->print("\n");
        }

        for my $gene (sort { ncmp($a->gene_name, $b->gene_name) } @seq_rna_genes) {
            my $gene_name = $gene->gene_name;
            my $accession = $gene->accession;
            my $codon = $gene->codon;
            my $amino_acid = $gene->amino_acid;
            my $amino_acid_code = substr($amino_acid, 0, 1);
            my $gene_score = $gene->score;
            my $source = $gene->source;

            my $start = $gene->start;
            my $end = $gene->end;
            ($start, $end) = ($end, $start) if $start > $end;
            
           # $ace_fh->print("Sequence $sequence\n");
           # $ace_fh->print("Subsequence $gene_name $start $end\n\n");
           # $ace_fh->print("Sequence : $gene_name\n");
           # $ace_fh->print("Source $sequence\n");

            my ($method, $remark, $transcript, $locus, $rfam_prod_line, @rfam_line);
            if ($source =~ /trnascan/i) {
                $method = 'tRNAscan';
                $remark = "\"tRNA-$amino_acid Sc=$gene_score\"";
                $transcript = "tRNA \"$codon $amino_acid $amino_acid_code\"";
                    $ace_fh->print("Sequence $sequence\n");
	            $ace_fh->print("Subsequence $gene_name $start $end\n\n");
        	    $ace_fh->print("Sequence : $gene_name\n");
        	    $ace_fh->print("Source $sequence\n");
                    
                    $ace_fh->print("Method $method\n") if defined $method;
                    $ace_fh->print("Remark $remark\n") if defined $remark;
                    $ace_fh->print("Locus $locus\n") if defined $locus;
                    $ace_fh->print("Transcript $transcript\n") if defined $transcript;
                    $ace_fh->print("\n");
            }
            elsif ($source =~ /rfam/i){
                
                if($gene_score >= 50) {
                $method = 'Rfam';
                $remark = "\"Predicted by Rfam ($accession), score $gene_score\"";
                $locus = $gene->description;
                $rfam_prod_line = `grep $locus /gscmnt/278/analysis/HGMI/RFAM/RFAM_v8.1.tbl | head -1`;
                       chomp $rfam_prod_line;
                       @rfam_line = split(/\t+/, $rfam_prod_line);
                    
 		    $ace_fh->print("Sequence $sequence\n");
	            $ace_fh->print("Subsequence $gene_name $start $end\n\n");
        	    $ace_fh->print("Sequence : $gene_name\n");
        	    $ace_fh->print("Source $sequence\n");                 

                    $ace_fh->print("Method $method\n") if defined $method;
                    $ace_fh->print("Remark $remark\n") if defined $remark;
                    $ace_fh->print("Locus $locus\n") if defined $locus;
                    $ace_fh->print("RFAM $gene_score\n") if defined $gene_score;
                    $ace_fh->print("Rfam_product \"$rfam_line[9]\" \n"); ##if defined $rfam_prod;
                    $ace_fh->print("Transcript $transcript\n") if defined $transcript;
                    $ace_fh->print("\n");
                  }
            }
            elsif ($source =~ /rnammer/i) {
                $method = 'RNAmmer';
                $remark = "\"Predicted by RNAmmer, score $gene_score\"";
                $locus = $gene->description;
               
		    $ace_fh->print("Sequence $sequence\n");
	            $ace_fh->print("Subsequence $gene_name $start $end\n\n");
        	    $ace_fh->print("Sequence : $gene_name\n");
        	    $ace_fh->print("Source $sequence\n");

                    $ace_fh->print("Method $method\n") if defined $method;
                    $ace_fh->print("Remark $remark\n") if defined $remark;
                    $ace_fh->print("Locus $locus\n") if defined $locus;
                    $ace_fh->print("Transcript $transcript\n") if defined $transcript;
                    $ace_fh->print("\n"); 
            }
            else {
                $method = $source;
                $remark = "\"Predicted by $method, score $gene_score\$";
                $locus = $gene->description;
               
                    $ace_fh->print("Sequence $sequence\n");
	            $ace_fh->print("Subsequence $gene_name $start $end\n\n");
        	    $ace_fh->print("Sequence : $gene_name\n");
        	    $ace_fh->print("Source $sequence\n");

                    $ace_fh->print("Method $method\n") if defined $method;
                    $ace_fh->print("Remark $remark\n") if defined $remark;
                    $ace_fh->print("Locus $locus\n") if defined $locus;
                    $ace_fh->print("Transcript $transcript\n") if defined $transcript;
                    $ace_fh->print("\n");
      
            }

#            $ace_fh->print("Method $method\n") if defined $method;
#            $ace_fh->print("Remark $remark\n") if defined $remark;
#            $ace_fh->print("Locus $locus\n") if defined $locus;
#            $ace_fh->print("Transcript $transcript\n") if defined $transcript;
#            $ace_fh->print("\n");
        }
    }

    $ace_fh->close;
    $self->status_message("Done, ace file is at " . $self->ace_file);
    return 1;
}

sub _get_sequences_from_file {
    my ($self, $file) = @_;
    my $seq_obj = Bio::SeqIO->new(
        -file => $file,
        -format => 'Fasta',
    );
    confess "Could not create Bio::SeqIO object for $file" unless $seq_obj;

    my @sequences;
    while (my $seq = $seq_obj->next_seq) {
        push @sequences, $seq->display_id;
    }

    return @sequences;
}

1;

