
package Genome::Model::Command::Export::Ace;

use strict;
use warnings;

use Genome;
use Command; 

use Bio::Index::Fasta;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
);

sub help_brief {
    "export data into a Ace file format"
}

sub help_synopsis {
    return <<"EOS"

Write a subclass of this.  

Give it a name which is an extension of this class name.

Implement a new writing out for some part  of a genome model.

EOS
}

sub help_detail {
    return <<"EOS"

This module is an base class for commands that export data into the Ace file format.

Subclasses will implement different input formats.  This module
should handle common parameters, typically for handling the output file. 

EOS
}


sub _print_result {
    my ($pos,$coverage) = @_;

    print "$pos:$coverage\n";
}

sub RefSeq {
  my ($f_refseq, $f_ace, $d_phds, $chromosome,$start_pos,$end_pos,$end_pad) = @_;
  my $refseq=&GetRefSeq(file=>$f_refseq,
		       chromosome=>$chromosome,
		       start=>$start_pos,
		       end=>$end_pos+$end_pad);
  $refseq->{name} = 'MosaikAnchor.C1';
  my $refseq_phd=$refseq->{name}.'.scf.phd.1';
	my $co_d_phds = (defined($d_phds)) ? $d_phds :
		File::Basename::dirname($f_ace) . '/../phd_dir';
	unless (-e $co_d_phds) {
		mkpath $co_d_phds;
	}
  &WritePhd(rd=>$refseq,f_phd=>$refseq_phd,d_phds=>$co_d_phds);
  $refseq->{name} = '.MosaikAnchor.C1';
  return $refseq;
}

sub WriteAce {
  my ($f_ace, $refseq, $ace_body_rd, $ace_body_qa, $num_sreads) = @_;
  my $i;
  my $rd;
  open(ACE,">$f_ace")|| die "unable to open $f_ace\n";
  printf ACE "AS 1 %d\n\n",$num_sreads+1;  #1 represents the refseq read
  printf ACE "CO %s %d %d 1 U\n",$refseq->{name}, $refseq->{len}, $num_sreads+1;
  &printSeq(fh=>\*ACE, seq=>$refseq->{seq}, size=>50);
  printf ACE "\nBQ\n";
  for($i=0;$i<$refseq->{len};$i++){
    print ACE " 30";
    if(($i+1)%50==0){print ACE "\n";}
  }
  print ACE "\n";
  printf ACE "\nAF %s U 1\n", $refseq->{name};
  close(ACE);
  system("cat $ace_body_rd >> $f_ace");
  unlink $ace_body_rd;

  open(ACE,">>$f_ace")|| die "unable to reopen $f_ace\n";
  printf ACE "BS 1 %d %s\n\n", $refseq->{len},$refseq->{name};
  
  printf ACE "RD %s %d 0 0\n",$refseq->{name},$refseq->{len};
  &printSeq(fh=>\*ACE, seq=>$refseq->{seq}, size=>50);
  printf ACE "\nQA 1 %d 1 %d\n",$refseq->{len},$refseq->{len};
  printf ACE "DS CHROMAT_FILE: %s PHD_FILE: %s.phd.1 TIME: %s\n\n",$refseq->{name},$refseq->{name}, $refseq->{timestamp};

  close(ACE);
  system("cat $ace_body_qa >> $f_ace");
  unlink $ace_body_qa;
}


sub WritePhd{
  #&WritePhd(refseq=>$refseq,d_phds=>$d_phds);
  my %arg=@_;
  my $d_phds=$arg{d_phds} if defined $arg{d_phds};
  if (!defined($d_phds)) {
    return;
  }
  my $rd=$arg{rd} if defined $arg{rd};
  my $f_phd=$arg{f_phd} if defined $arg{f_phd};
  my $prboffset = (defined($arg{prboffset})) ? $arg{prboffset} : 64;

  open(PHD,">$d_phds/$f_phd") || die "unable to open $d_phds/$f_phd\n";
  printf PHD "BEGIN_SEQUENCE %s\n", $f_phd;
  print PHD "BEGIN_COMMENT\n";
  printf PHD "TIME: %s\n\n", $rd->{timestamp};
  print PHD "END_COMMENT\n\n";
  print PHD "BEGIN_DNA\n";
  my @seq=split //, $rd->{seq};
  my @prb;
  my $fakeprb=0;
  if(defined $rd->{prb}){
    @prb=split //, $rd->{prb};
  }
  else{
#    $self->error_message("$f_phd does not have prb values, using: 30");
    $fakeprb=30;
  }

  for(my $i=0;$i<=$#seq;$i++){
    my $qual;
    if($fakeprb>0){
      $qual=$fakeprb;
    }
    else{
      $qual=ord($prb[$i])-$prboffset;
    }
    printf PHD "%s %d %d\n", uc($seq[$i]),$qual,$i+1;
  }


  print PHD "END_DNA\n";
  print PHD "END_SEQUENCE\n";
  close(PHD);
}

sub printSeq{
  my %arg=@_;
  my $fh=$arg{fh} if defined $arg{fh};
  my $seq=$arg{seq} if defined $arg{seq};
  my $size=$arg{size} if defined $arg{size};

  my $len=length($seq);
  for(my $i=0;$i<$len;$i+=$size){
    my $slen=$len-$i;
    $slen=($slen>$size)?$size:$slen;
    my $str=unpack("x$i a$slen", $seq);
    printf $fh "%s\n",$str;
  }

  #return substr($self->sequence_base_string, $start, $length);
  #return unpack("x$start a$length", $self->sequence_base_string);
}


# here is where the retrieval key is specified for Bio::Index::Fasta in GetRefSeq
sub get_chromosome_id {
  my $line = shift;
	my ($chromosome);
	if ($line =~ / NC_0000 ([^\.]+)/x) {
		$chromosome = $1;
		$chromosome =~ s/^0+//;	# Strip leading zero
	} else {
		$chromosome = $line;
	}
  return $chromosome;
}

sub GetRefSeq{
  my %arg=@_;
  my $f_refseq=$arg{file} if defined $arg{file};
  $f_refseq=~/(\S+)\.fa/;
  my $name=$1;

  my $index_filename = $f_refseq . '.dbm';
  my $fasta_index;
  if (-r $index_filename) {
    $fasta_index = Bio::Index::Fasta->new(
					     '-filename' => $index_filename,
					    );
  } else {
    $fasta_index = Bio::Index::Fasta->new(
					     '-filename' => $index_filename,
					     '-write_flag' => 1,
					    );
    $fasta_index->id_parser(\&get_chromosome_id);
    # make the index
    $fasta_index->make_index($f_refseq);
  }

  my $ref_seq_name = $arg{chromosome};
  $ref_seq_name = ($fasta_index->get_all_primary_ids)[0];
  my $ref_seq_obj = $fasta_index->fetch($ref_seq_name);

  unless (defined($ref_seq_obj)) {
#    $self->error_message("Failed to fetch sequence for $ref_seq_name");
		return 0;
  }

  my $seq = $ref_seq_obj->trunc($arg{start}, $arg{end})->seq();

  my $now_string = localtime; # e.g. "Thu Oct 13 04:54:34 1994"
  my $len=length($seq);
  my %refseq=(name=>$name,
	      seq=>$seq,
	      len=>$len,
	      start=>$arg{start},
	      end=>$arg{end},
	      timestamp=>$now_string
	     );
  return \%refseq;
}



1;

