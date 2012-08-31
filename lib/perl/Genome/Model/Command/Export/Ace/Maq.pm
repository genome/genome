
# Rename the final word in the full class name to match the filename <---
package Genome::Model::Command::Export::Ace::Maq;

use strict;
use warnings;

use Genome;
use Command;

use Genome::Model::Command::Export::Ace;
use File::Basename;
use File::Temp;
use File::stat;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',                       
    has => [                                # Specify the command's properties (parameters) <--- 
        'maq'   => { type => 'String',      doc => "the maq map (input) file" },
        'ace'   => { type => 'String',      doc => "ace output file" },
        'refseq'   => { type => 'String',      doc => "the reference sequence fasta (input) file" },
        'chromosome'   => { type => 'String',      doc => "the chromosome to convert" },
        'start'   => { type => 'Integer',      doc => "the starting position to convert" },
        'end'   => { type => 'Integer',      doc => "the ending position to convert" },

        'phds'   => { type => 'String',      doc => "a directory to output phd files", is_optional => 1 },
        'endpad'   => { type => 'Integer',      doc => "number of extra bases to get at the end of the reference sequence", is_optional => 1 },

        'anychr'   => { type => 'Boolean',     doc => "output any chromosome", is_optional => 1 }
    ], 
);

sub help_brief {                            # Keep this to just a few words <---
    "converts from maq to ace (with phd files)"
}

sub help_synopsis {                         # Replace the text below with real examples <---
    return <<EOS
genome-model write ace maq --ace=aml.ace --maq=aml.map --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000
genome-model write ace maq --ace=aml.ace --maq=aml.map --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000 --endpad=0
genome-model write ace maq --ace=aml.ace --maq=aml.map --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000 --anychr
EOS
}

sub help_detail {                           # This is what the user will see with the longer version of help. <---
    return <<EOS 
This script converts from Maq to ace (with phd files).
EOS
}

sub execute {                               # Replace with real execution logic.
	my $self = shift;
	
	return unless (
								 defined($self->ace) && defined($self->refseq) && defined($self->maq) &&
								 defined($self->chromosome) && defined($self->start) && defined($self->end)
								);
	my $endpad = $self->endpad;
	$endpad ||= 33;
	
	my $refseq =
		Genome::Model::Command::Export::Ace::RefSeq(
																							 $self->refseq, $self->ace, $self->phds,
																							 $self->chromosome,$self->start,$self->end,$endpad);
	my ($ace_body_rd, $ace_body_qa, $num_reads) =	$self->MAQ_Ace();
	Genome::Model::Command::Export::Ace::WriteAce($self->ace,
																							 $refseq, $ace_body_rd, $ace_body_qa, $num_reads);
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub MAQ_Ace {
  my ($self) = @_;
  my ($maq_file, $d_phds, $chromosome,$start_pos,$end_pos) = 
		($self->maq, $self->phds,$self->chromosome,$self->start,$self->end);
	my($maq_prboffset) = 33;
  my $ace_body_rd_fh = new File::Temp( UNLINK => 0) or return 0; # $logger->fatal("Unable to open temporary ace read file");
  my $ace_body_rd = $ace_body_rd_fh->filename();
	
  my $ace_body_qa_fh = new File::Temp( UNLINK => 0) or return 0; # $logger->fatal("Unable to open temporary ace qa file");
  my $ace_body_qa = $ace_body_qa_fh->filename();
  my $num_reads = 0;
	
  my ($maq_mapview) = $ENV{GENOME_SW} . "/maq/maq-0.5.1_i686-linux/maq mapview $maq_file";
  open(MAQVIEW,"$maq_mapview |") or return 0; # $logger->fatal("Unable to run command: $maq_mapview");
  my ($bn) = basename($maq_file) =~ /(\S+)\.map/;
  my $rec_id = 0;
  while(<MAQVIEW>) {
    chomp;
    my ($read_id, $ref_id, $position, $strand, $relative_position_paired_read, $paired_flag,
				$mapping_quality, $single_end_mapping_quality, $alternative_mapping_quality,
				$number_of_mismatches_best_hit, $number_of_mismatches_second_hit,
				$read_size, $sequence, $quality) = split("\t");
    my ($read_chromosome);
    if ($ref_id =~ / NC_0000 ([^\.]+)/x ) {
			$read_chromosome = $1;
			$read_chromosome =~ s/^0+//;	# Strip leading zero
		} else {
			$read_chromosome = $ref_id;
		}
    if (($self->anychr || $read_chromosome eq $chromosome) &&
				$position >= $start_pos &&
				$position <= $end_pos) {
      my $nalignment = 1;
      my $now_string = localtime; # e.g. "Thu Oct 13 04:54:34 1994"
      
      ++$rec_id;
      ++$num_reads;
			#      my $rdname=$bn . '_' . $read_id;
      my $rdname=$read_id;
			$rdname =~ s/^Runs_//;
      my $f_phd=$rdname .'.phd.1';
      my $orientation = ($strand eq '+') ? 'U' : 'C';
      my $read_start = ($position - $start_pos) + 1;
      Genome::Model::Command::Export::Ace::WritePhd(rd=>
																									 {
																										seq=>$sequence,
																										qryseq=>$sequence,
																										orientation=> $orientation,
																										naln=>$nalignment,
																										saln=>$mapping_quality,
																										scpos=>$read_start,
																										prb=>$quality,
																										timestamp=>$now_string
																									 },
																									 f_phd=>$f_phd,d_phds=>$d_phds,prboffset=>$maq_prboffset);
      printf $ace_body_rd_fh "AF %s %s %d\n",$rdname,$orientation,$read_start;
      my $padded_rdlen=length($sequence);
      printf $ace_body_qa_fh "RD %s %d 0 0\n",$rdname,$padded_rdlen;
      Genome::Model::Command::Export::Ace::printSeq(fh=>$ace_body_qa_fh, seq=>$sequence,size=>50);
      printf $ace_body_qa_fh "\nQA 1 %d 1 %d\n",$padded_rdlen,$padded_rdlen;
      printf $ace_body_qa_fh "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TIME: %s\n\n",$rdname,$rdname, $now_string;
    }
  }
  $ace_body_rd_fh->close();
  $ace_body_qa_fh->close();
  close(MAQVIEW);
	
  return ($ace_body_rd, $ace_body_qa, $num_reads);
}

1;

