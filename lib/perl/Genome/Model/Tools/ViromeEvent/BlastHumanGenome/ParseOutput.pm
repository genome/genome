package Genome::Model::Tools::ViromeEvent::BlastHumanGenome::ParseOutput;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;
use Bio::SeqIO;
use Bio::SearchIO;

class Genome::Model::Tools::ViromeEvent::BlastHumanGenome::ParseOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "gzhao's Blast Human Genome parse output";
}

sub help_detail {
    return <<"EOS"
This script will check all .HGblast.parsed file in the 
given directory to make sure parsing blastn output file is finished 
for the given file.
EOS
}

sub execute {
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename ($dir);

    my $hg_blast_dir = $dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq_HGblast';
    unless (-d $hg_blast_dir) {
	$self->log_event("Failed to find HG blastN output dir for $sample_name");
	return;
    }
    my @fa_files = glob("$hg_blast_dir/*fa");

    if (@fa_files == 0) {
	$self->log_event("No more reads available for $sample_name");
	return 1;
    }

    foreach my $file (@fa_files) {
	next if $file =~ /HGfiltered\.fa$/; #ALREADY RAN BLAST OUT FILES IF PRESENT
	next if $file =~ /HGhits\.fa$/;  #FASTA OF BLAST HITS
	my $root = $file;
	$root =~ s/\.fa$//;
	my $blast_out = $root.'.HGblast.out';
	my $blast_parse = $root.'.HGblast.parsed';
	my $blast_filtered = $root.'.HGfiltered.fa';
	my $blast_hits = $root.'.HGhits.fa';
	my $blast_out_file_name = basename($blast_out);

	unless (-s $blast_out) { #SHOULD HAVE BEEN CREATED IN PREVIOUS STEP
	    $self->log_event("Failed to find HG blastN out file for $file");
	    return;
	}
	if (-s $blast_parse && -s $blast_filtered){
	    #PARSE BLAST OUT FILE ALREADY DONE
	    my $out = `tail -n 3 $blast_parse`;
	    if ($out =~ /\s+Summary:/) {
		$self->log_event("Parsing $blast_out_file_name already done");
		next;
	    }
	}

	unlink $blast_parse if -s $blast_parse;
	unlink $blast_filtered if -s $blast_filtered;
	unlink $blast_hits if -s $blast_hits;

	$self->log_event("Parsing $blast_out_file_name");

	#GET READ NAMES THAT WILL NOT BE KEPT
	my $filtered_reads = $self->parse_blast_file($blast_out);
	unless (scalar @$filtered_reads > 0) {
	    $self->log_event("Warning: No reads filtered for ".basename($blast_out));
	}

	my $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$blast_filtered");
	my $hits_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$blast_hits");
	my $in_io = Bio::SeqIO->new(-format => 'fasta', -file => $file);

	while (my $seq = $in_io->next_seq) {
	    my $read_name = $seq->primary_id;
	    unless (grep (/^$read_name$/, @$filtered_reads)) {
		$out_io->write_seq($seq);
	    }
            else
            {
                $hits_io->write_seq($seq);
            }
	}

	unless (-s $blast_filtered) {
            $self->log_event("No human reads found in $blast_out_file_name");
	    #return;
	}
    }

    return 1;
}

sub parse_blast_file {
    my ($self, $blast_file) = @_;

    my $parse_out_file = $blast_file;
    $parse_out_file =~ s/out$/parsed/;
    unlink $parse_out_file;
    my $parse_out_fh = Genome::Sys->open_file_for_writing( $parse_out_file );

    # get report from blast out file
    my %report_params = (
        blast_out_file => $blast_file,
        blast_type     => 'blastn',
    );
    my $report = $self->get_blast_report( %report_params );
    if ( not $report ) {
        $self->log_event('Failed get blastN blast report');
        return;
    } 

    my @filtered_reads; #READS THAT WILL NOT BE KEPT
    my $e_cutoff = 1e-10;
    my $read_count = 0;
    while (my $result = $report->next_result) {
	$read_count++;
	while (my $hit = $result->next_hit) {
	    if ($hit->significance <= $e_cutoff) {
		$parse_out_fh->print ( $result->query_name."\t".
				       $result->query_length."\tHomo\tHomo\t".
				       $hit->name."\t".$hit->significance."\n" );
		
		push @filtered_reads, $result->query_name; #READS THAT WILL NOT GO ON TO NEXT STAGES
		last; #ONLY LOOKING AT THE THE FIRST HIT
	    }
	}
    }
    my $keep_reads = $read_count - scalar @filtered_reads;
    $parse_out_fh->print("# Summary: $keep_reads out of $read_count is saved for BLASTN analysis\n");
    $parse_out_fh->close;
    return \@filtered_reads;
}

1;
