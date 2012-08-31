
# Rename the final word in the full class name to match the filename <---
package Genome::Model::Command::Export::Ace::Sxog;

use strict;
use warnings;

use Genome;
use Command;

use Genome::Model::Command::Export::Ace;
use File::Path;
use File::Basename;
use File::Temp;
use IO::File;
use PDL;
use PDL::IO::FlexRaw;
use File::stat;

do { no warnings; $PDL::BIGPDL = 1; };

my $LINKED_LIST_RECORD_SIZE = 8 + 8 + 8 + 1 + 1 + 2 + 60; 

my $LINKED_LIST_NODE_HEADER = [
                {Type=> 'longlong', NDims => 1, Dims=>1},
                {Type=> 'longlong', NDims => 1, Dims=>1},
                {Type=> 'double', NDims => 1, Dims=>1},
                {Type=> 'byte', NDims => 1, Dims=>1},
                {Type=> 'byte', NDims => 1, Dims=>1},
                {Type=> 'ushort', NDims => 1, Dims=>1},
                {Type=> 'byte', NDims => 1, Dims => [ 60 ]}
               ];

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',                       
    has => [                                # Specify the command's properties (parameters) <--- 
        'sxogs'   => { type => 'String',  doc => "the (input) directory for the sxog tile search binary files" },
        'prbfile' => { type => 'String',  doc => "the (input) filename for the (binary) prb values" },
        'ace'   => { type => 'String',      doc => "ace output file" },
        'refseq'   => { type => 'String',      doc => "the reference sequence fasta (input) file" },
        'chromosome'   => { type => 'String',      doc => "the chromosome to convert" },
        'start'   => { type => 'Integer',      doc => "the starting position to convert" },
        'end'   => { type => 'Integer',      doc => "the ending position to convert" },

        'phds'   => { type => 'String',      doc => "a directory to output phd files", is_optional => 1 },
        'endpad'   => { type => 'Integer',      doc => "number of extra bases to get at the end of the reference sequence", is_optional => 1 },

        'mosaik'   => { type => 'Boolean',     doc => "output mosaik compatible ace file", is_optional => 1 },
        'anychr'   => { type => 'Boolean',     doc => "output any chromosome", is_optional => 1 }
    ], 
);

sub help_brief {                            # Keep this to just a few words <---
    "converts from sxog to ace (with phd files)"
}

sub help_synopsis {                         # Replace the text below with real examples <---
    return <<EOS
cd edit_dir

genome-model write ace sxog -ace=edit_dir/aml.ace -sxogs=sxog_1 --prbfile=../../../../runs/solexa/rawprb.dat --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000
genome-model write ace sxog -ace=edit_dir/aml.ace -sxogs=sxog_1 --prbfile=../../../../runs/solexa/rawprb.dat --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000 --endpad=0
genome-model write ace sxog -ace=edit_dir/aml.ace -sxogs=sxog_1 --prbfile=../../../../runs/solexa/rawprb.dat --chromosome=7 --refseq=Hs_build36.fa --start=1 --end=500000 --mosaik
EOS
}

sub help_detail {                           # This is what the user will see with the longer version of help. <---
    return <<EOS 
This script converts from SynmatixOligoSearchG binary (gsc) format to ace (with phd files).
EOS
}

sub execute {                               # Replace with real execution logic.
	my $self = shift;
	
	return unless (
								 defined($self->ace) && defined($self->refseq) &&
								 defined($self->sxogs) && defined($self->prbfile) &&
								 defined($self->chromosome) && defined($self->start) && defined($self->end)
								);
	my $endpad = $self->endpad;
	$endpad ||= 33;
	
	my $refseq =
		Genome::Model::Command::Export::Ace::RefSeq(
																							 $self->refseq, $self->ace, $self->phds,
																							 $self->chromosome,$self->start,$self->end,$endpad);
	my ($ace_body_rd, $ace_body_qa, $num_reads) =	$self->SXOG_Ace();
	Genome::Model::Command::Export::Ace::WriteAce($self->ace,
																							 $refseq, $ace_body_rd, $ace_body_qa, $num_reads);
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub SXOG_Ace {
  my ($self) = @_;
  my ($sxog_base,$d_phds, $prb_filename, $chromosome,$start_pos,$end_pos,$output_mosaik) =
		($self->sxogs, $self->phds, $self->prbfile,
		 $self->chromosome,$self->start,$self->end, $self->mosaik);

  my $ace_body_rd_fh = new File::Temp( UNLINK => 0) or return 0; #$logger->fatal("Unable to open temporary ace read file");
  my $ace_body_rd = $ace_body_rd_fh->filename();

  my $ace_body_qa_fh = new File::Temp( UNLINK => 0) or return 0; #$logger->fatal("Unable to open temporary ace qa file");
  my $ace_body_qa = $ace_body_qa_fh->filename();
  my $num_reads = 0;

  my $prb = OpenPRB($prb_filename,$sxog_base);

  my $rds = SXOGReads($sxog_base,$chromosome,$start_pos,$end_pos,$output_mosaik);
	my %reads;
  foreach my $rd (@{$rds}) {
    my ($read_id, $strand, $position, $nalignment, $mapping_quality) = @{$rd};
    my ($qry_sequence, $quality, $ltxy_read_id) = GetPRBFQ($prb,$read_id);

    if ($output_mosaik) {
			$read_id = $ltxy_read_id;
		} else {
			my $read_id_suffix = (exists($reads{$read_id})) ?
				sprintf "_%04u",$reads{$read_id} : '';
			$reads{$read_id} += 1;
			$read_id .= $read_id_suffix;
		}

    my $rdname = $read_id;
    my $f_phd=$rdname .'.phd.1';
    my $orientation = ($strand) ? 'U' : 'C';
    # If the orientation is opposite, then reverse and complement
    # the sequence from the PRB file
    if ($strand == 0) {
      $qry_sequence =~ tr/ACGT/TGCA/;
      $qry_sequence = reverse($qry_sequence);
    }
    my $read_start = ($position - $start_pos) + 1;
    my $now_string = localtime; # e.g. "Thu Oct 13 04:54:34 1994"
    
    Genome::Model::Command::Export::Ace::WritePhd(rd=>
	      {
	       seq=>$qry_sequence,
	       qryseq=>$qry_sequence,
	       orientation=> $orientation,
	       naln=>$nalignment,
	       saln=>$mapping_quality,
	       scpos=>$read_start,
	       prb=>$quality,
	       timestamp=>$now_string
	      },
	      f_phd=>$f_phd,d_phds=>$d_phds);
    my %read=(name=>$rdname,
	      scpos=>$read_start,
	      ori=>$orientation,
	      readseq=>$qry_sequence,
	      timestamp=>$now_string
	     );
    $num_reads++;
    printf $ace_body_rd_fh "AF %s %s %d\n",$rdname,$orientation,$read_start;
    my $padded_rdlen=length($qry_sequence);
    printf $ace_body_qa_fh "RD %s %d 0 0\n",$rdname,$padded_rdlen;
    Genome::Model::Command::Export::Ace::printSeq(fh=>$ace_body_qa_fh, seq=>$qry_sequence,size=>50);
    printf $ace_body_qa_fh "\nQA 1 %d 1 %d\n",$padded_rdlen,$padded_rdlen;
    printf $ace_body_qa_fh "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TIME: %s\n\n",$rdname,$rdname, $now_string;
  }
  $ace_body_rd_fh->close();
  $ace_body_qa_fh->close();
  ClosePRB($prb);
  return ($ace_body_rd, $ace_body_qa, $num_reads);
}


sub SXOGReads {
  my ($sxog_base,$chromosome,$start_pos,$end_pos,$output_mosaik) = @_;
  my ($sxog_ndx_file) = $sxog_base . '_' . $chromosome . '_aln.ndx';
  my ($sxog_dat_file) = $sxog_base . '_' . $chromosome . '_aln.dat';
  my @reads;

  my $index_list_size = File::stat::stat($sxog_ndx_file)->size;
  my $max_end = int($index_list_size/8);
  $max_end = ($end_pos > $max_end) ? $max_end : $end_pos;
  my $position_slice = $end_pos - $start_pos;
  
  my $INDEX_NODE_HEADER = [ {Type=> 'longlong', NDims => 1, 
			     Dims => [ $position_slice+1 ]} ];
  
  my $alignment_list_fh = IO::File->new($sxog_dat_file);
  #$logger->fatal("Unable to open: $sxog_dat_file")
	return 0 unless $alignment_list_fh;
  
  my $alignment_list_size = File::stat::stat($sxog_dat_file)->size;
  
  my $sxog_ndx_fh = IO::File->new($sxog_ndx_file);
  seek($sxog_ndx_fh, $start_pos * 8,SEEK_SET);
  my $index_pdl = readflex($sxog_ndx_fh,$INDEX_NODE_HEADER);
  $sxog_ndx_fh->close();
  my $len = $index_pdl->nelem();
  
  for (my $cur_pos=0; $cur_pos < $len; $cur_pos++) {
    my $next_alignment_num = $index_pdl->at($cur_pos);
    if ($next_alignment_num) {
      do {
	seek($alignment_list_fh, $next_alignment_num * $LINKED_LIST_RECORD_SIZE,SEEK_SET);
	my ($last_alignment_num, $read_num, $prob, $read_len, $orient, $number_alignments, $pdl_mismatch_string) =
	  readflex( $alignment_list_fh, $LINKED_LIST_NODE_HEADER);
	
	# post-process the above PDL data
	# this might go into the get_alignment_node_for_alignment_num() call above
	($last_alignment_num, $read_num, $prob, $read_len, $orient, $number_alignments)
	  = map {$_->at(0)} ($last_alignment_num, $read_num, $prob, $read_len, $orient, $number_alignments);
	unless ($output_mosaik && $number_alignments > 1) {
		push @reads, [ $read_num, $orient, $cur_pos+$start_pos, $number_alignments, $prob ];
	}
	
	# This alignment record points to the previous one in the file, which is the next one we'll grab.
	# Since alignment positions are 1-based, alignment 0 is the end of the linked-list.
	$next_alignment_num = $last_alignment_num;
      } while($next_alignment_num);
    }
  }
  return \@reads;
}

sub GetPRBFQ {
  my ($prb,$read_id) = @_;
  my ($prb_fh, $ltxy_fh, $prb_hdr, $map_ref, $read_length, $max_read_number) = @{$prb};
  my ($prb_offset) = ( $read_id % 1000000000);
  $prb_offset = FindSXOGPRBOffset($map_ref,$prb_offset);
  $prb_fh->seek($read_length * 4 * ($prb_offset-1), SEEK_SET);
  my $pdlread = readflex($prb_fh,$prb_hdr);
  my @basetrans = ( 'A', 'C', 'G', 'T' );

  my ($sequence,$quality);
  foreach my $position (dog $pdlread) {
    my @base = map { $_ = ($_ > 127) ? $_ - 256 : $_; } list $position;
    my $maxscore = -100;
    my $maxbase = '';
    if (join('',@base) ne '0000') {	# Do not put out the blank values
      for(my $j=0;$j < 4;$j++) {
	if ($base[$j] > $maxscore) {
	  $maxscore = $base[$j];
	  $maxbase = $basetrans[$j];
	}
      }
      $sequence .= $maxbase;
      $quality .= chr($maxscore+64);
    }
  }

  $ltxy_fh->seek(8 * ($prb_offset-1), SEEK_SET);
  my $ltxy_hdr = [ { Type=>'short', NDims=>1,Dims=>[4]} ];
  my $pdlltxy = readflex($ltxy_fh,$ltxy_hdr);
	my ($ltxy_lane, $ltxy_tile, $ltxy_x, $ltxy_y) = list $pdlltxy;
	my $ltxy_read_id =
		$ltxy_lane . '_' . $ltxy_tile . '_' . $ltxy_x . '_' . $ltxy_y;

  return ($sequence, $quality, $ltxy_read_id);
}

sub ReadNdxMap {
  my ($prb_filename,$sxog_base) = @_;
  my ($sxog_reads_file) = $sxog_base . '_read_ndx.tsv';
  my ($sxog_reads_map_file) = $sxog_base . '_read_ndx.map';

  my %map_offset;
  if (-r $sxog_reads_map_file) {
    open(READMAP,"$sxog_reads_map_file") or return 0; #$logger->fatal("Unable to open read index map file for writing: $sxog_reads_file");
    while(<READMAP>) {
      chomp;
      my ($sxog_lane,$sxog_tile,$read_num) = split("\t");
      $map_offset{"$sxog_lane\t$sxog_tile"} = $read_num;
    }
    close(READMAP);
  } else {
    open(READMAP,">$sxog_reads_map_file") or return 0; #$logger->fatal("Unable to open read index map file for writing: $sxog_reads_file");
    open(READNDX,$sxog_reads_file) or return 0; #$logger->fatal("Unable to open read index file: $sxog_reads_file");
    while(<READNDX>) {
      chomp;
      my ($read_id, $prb_file) = split("\t");
      my ($sxog_lane, $sxog_tile) = $prb_file =~ /_(\d+)_(\d+)/;
      if (!exists($map_offset{"$sxog_lane\t$sxog_tile"})) {
	my ($read_num) = ( $read_id % 1000000000);
	$read_num -= 1;
	$map_offset{"$sxog_lane\t$sxog_tile"} = $read_num;
	print READMAP "$sxog_lane\t$sxog_tile\t$read_num\n";
      }
    }
    close(READMAP);
    close(READNDX);
  }
  my $map_filename = $prb_filename;
  $map_filename =~ s/dat$/map/;
  open(PRBMAP,$map_filename);
  my %map;
  while(<PRBMAP>) {
    chomp;
    my ($lane, $tile, $offset) = split("\t");
    if (!exists($map_offset{"$lane\t$tile"})) {
			return 0;
      #$logger->error("Missing lane: $lane and tile: $tile");
    } else {
      my $sxog_offset = $map_offset{"$lane\t$tile"};
      $map{$sxog_offset} = $offset - $sxog_offset;
    }
  }
  close(PRBMAP);

  return \%map;
}

sub FindSXOGPRBOffset {
  my ($map_ref, $read_num) = @_;
  my %map = %{$map_ref};
  my $last_offset;
  foreach my $sxog_offset (sort { $a <=> $b } (keys %map)) {
    $last_offset = $map{$sxog_offset};
    if ($sxog_offset >= $read_num) {
      return $last_offset + $read_num;
    }
  }
  return $read_num;
}

sub OpenPRB {
  my ($prb_filename,$sxog_base) = @_;
  my ($prb_fh) = new IO::File();
  my ($ltxy_fh) = new IO::File();

  my $map_ref = ReadNdxMap($prb_filename,$sxog_base);
  # Open the header and get the read length and max read number
  my ($read_length, $max_read_number);
  open(PRBHDR,$prb_filename . '.hdr');
  while(<PRBHDR>) {
    if ( / ^ \s* byte \s+ 3 \s+ 4 \s+ (\d+) \s+ (\d+) /x ) {
      ($read_length, $max_read_number) = ($1, $2);
    }
  }
  close(PRBHDR);

	my $ltxy_filename = $prb_filename;
	$ltxy_filename =~ s/\.dat/_ltxy.dat/;
  $ltxy_fh->open($ltxy_filename) or return 0;
	#$logger->fatal("Unable to open ltxy data file: $ltxy_filename");
  # Use the header file values
  my $prb_hdr = [ { Type=>'byte', NDims=>2,Dims=>[4,$read_length]} ];
  $prb_fh->open($prb_filename) or return 0; #$logger->fatal("Unable to open prb data file: $prb_filename");
  $prb_fh->binmode();
  return [ $prb_fh, $ltxy_fh, $prb_hdr, $map_ref, $read_length, $max_read_number ];
}

sub ClosePRB {
  my ($prb) = @_;
  my ($prb_fh, $ltxy_fh) = @{$prb};
  $prb_fh->close();
  $ltxy_fh->close();
}

1;

