package Genome::Model::Tools::Hgmi::CollectSequence;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use Bio::Seq;
use Bio::SeqIO;
use English;

use File::Copy;
use Cwd;

use Compress::Bzip2;
use File::Basename;
use File::Spec;
use File::Temp;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
    has => [
	    'seq_file_name'  => { is => 'String',
				  doc => "Fasta file of contigs", 
				},
	    'seq_file_dir'   => { is => 'String',
				  doc => "Where to retrieve the fasta file from",
				},
	    'minimum_length' => { is => 'Integer',
				  doc => "Minimum contig length", 
				},
	    'new_ctgs_out'   => { is => 'String',
				  doc => "Contigs greater than minimum_length unless part of a scaffold ",
				  is_optional => 1,
				},

	   ]
			);

sub help_brief
{
    "tool for picking out only the sequences of specified length unless part of a scaffold";
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
need to put help synopsis here
EOS
}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
need to put help detail here.
EOS
}


sub execute
{
    my $self = shift;

    my $seq_file_name = $self->seq_file_name;
    my $seq_file_dir  = $self->seq_file_dir;

    my $cwd = getcwd();

    my $old_file = qq{$seq_file_dir/$seq_file_name};

    my $new_file = qq{$cwd/$seq_file_name};

    if (-e $old_file) {

      copy($old_file, $new_file);

    }
    else {

      die qq{ \n\n No $old_file found! CollectSequence.pm: $OS_ERROR \n\n };

    }

    my (

	$quals,
	$qualfile,
	$newqualfile,
	$agp,
	$agp_file,
	$newagp_file,

       );

    $quals       = qq{contigs.quals};
    $qualfile    = qq{$seq_file_dir/$quals};
    $newqualfile = qq{$cwd/contigs.quals};

    if (( -e $qualfile ) && ( ! -z $qualfile )) {

      copy($qualfile, $newqualfile);

    }
    else {

      warn qq{ \n\n No $qualfile found! CollectSequence.pm\n\n };

    }

    $agp         = qq{supercontigs.agp};
    $agp_file    = qq{$seq_file_dir/$agp};
    $newagp_file = qq{$cwd/supercontigs.agp};

    if (( -e $agp_file ) && (! -z $agp_file )) {

      copy($agp_file, $newagp_file);

    }
    else {

      warn qq{\n\n No $agp_file found! CollectSequence.pm\n\n};

    }

    my %lookup    = ( );
    my %sequences = ( );
    my @Counter   = ( );

    my $in_file   = $seq_file_name;
    my $out_file  = qq{contigs_gt_th_200.fasta};
    my $out_short = qq{contigs_ls_th_200.fasta};

    $self->new_ctgs_out($out_file);

    my $seqin    = new Bio::SeqIO(-file => $in_file, -format => "fasta");
    my $seqout   = new Bio::SeqIO(-file => ">$out_file", -format => "fasta");
    my $seqshort = new Bio::SeqIO(-file => ">$out_short", -format => "fasta");

    while ( my $seq = $seqin->next_seq() ) {

      my ($superctg, $chunk)  = split(/\./, $seq->primary_id());
      my $ctglength           = length($seq->seq());
      my $primecontig         = $seq->primary_id();
      my $sequence            = $seq->seq();

      $lookup{$superctg}{$primecontig} += $ctglength;

      push @Counter, $superctg;

      $sequences{$primecontig} = $sequence;

    }

    foreach my $supercontigs ( sort keys %lookup ) {

      ##perhaps I should make getting these counts a subroutine later...

      my %count = ( );

      foreach $supercontigs ( @Counter ) {

	$count{$supercontigs}++;

      }

      foreach my $contigs ( sort keys %{$lookup{$supercontigs}} ) {

	my $lengthofctg = $lookup{$supercontigs}{$contigs};
	my $sequence    = $sequences{$contigs};
	my $seq         = Bio::Seq->new ( -id => $contigs, -seq => $sequence);

	if ( $lengthofctg >= $self->minimum_length ) {

	  $seqout->write_seq($seq);

	}
	else {

	  my $ctgcount = 0;

	  $ctgcount = $count{$supercontigs};

	  if ( $ctgcount >= 2 ) {

	    $seqout->write_seq($seq);

	    unless ( $lengthofctg > 41 ) {

	      warn qq{\n\nWARNING: $contigs is $lengthofctg bps long and is equal to or less than 40bps! DO not run contig! CollectSequence.pm\n\n};

	    }

	  }
	  else {

	    $seqshort->write_seq($seq);

	  }

	}

      }

    }

    my @clean_up = ($seq_file_name, $quals);

    foreach my $file (@clean_up) {

      unless (defined($file)) {

	warn "$file is not defined in CooectSequence.pm for bz2 operations!";
	next;

      }

      if ( ( -e $file ) && ( ! -z $file)) {

	my $bz_file;
	my $end_file = File::Spec->catfile($cwd, "$file.bz2");

	open(FILE, $file) or
	  die "Can not open $file, CollectSequence.pm:  $OS_ERROR \n";

	$bz_file = bzopen($end_file, 'w') or
	  die qq "Can not open '$end_file' from CollectSequence.pm : $bzerrno\n";

	while (my $line = <FILE>) {

	  $bz_file->bzwrite($line) or
	    die "Error writing bzip2 file from CollectSequence.pm: $bzerrno! \n\n";

	}

	$bz_file->bzclose();
	close(FILE);

	unlink($file) or
	  die "Can't delete $file from CollectSequence.pm : $OS_ERROR \n";
      }
      else {

	warn "$file is not present or Zero file size in CollectSequence.pm for bz2 operations!";
	next;

      }

    }

    return 1;
  }

1;
