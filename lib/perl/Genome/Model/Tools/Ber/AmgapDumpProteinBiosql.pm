package Genome::Model::Tools::Ber::AmgapDumpProteinBiosql;

use strict;
use warnings;

use Genome;
use Command;

use Carp;
use English;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;

use Data::Dumper;
use IO::File;


UR::Object::Type->define(
			 class_name => __PACKAGE__,
			 is         => 'Command',
			 has        => [
					'locus_tag'       => {
							      is => 'String',
							      doc => "Locus tag for project, containing DFT/FNL postpended",
							     },
					'proteinnamefile' => { is  => 'String',
							       doc => "Name of proteins dumped out from biosql file",
							       is_optional => 1,
							   },
				       ],
             has_optional => [
                   'outputdir' => { is => 'String',
                                    doc => 'output directory files are to be written in',
                                  },
                   'dev' => { is => 'Boolean',
                              doc => 'dump from development',
                              default => 0,
                    },
                   ],
			);

sub help_brief
  {
    "Tool for dumping the protein files from BioSQL for BER product naming pipeline";
  }

sub help_synopsis
  {
    return <<"EOS"
      Tool for dumping the protein files from BioSQL for BER product naming pipeline.
EOS
  }

sub help_detail
  {
    return <<"EOS"
Tool for dumping the protein files from BioSQL for BER product naming pipeline.
EOS
  }


sub execute
  {
    my $self = shift;
    my @basket = ( );
    my $locus_tag = $self->locus_tag;

    my $dbadp;
    unless($self->dev) {
        $dbadp = Bio::DB::BioDB->new(
				    -database => 'biosql',
				    -user     => 'sg_user',
				    -pass     => 'sg_us3r',
				    -dbname   => 'DWRAC',
				    -driver   => 'Oracle'
				   );
    }
    else
    {
        $self->status_message("using development biosql");
        $dbadp = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sgus3r',
            -dbname   => 'DWDEV',
            -driver   => 'Oracle',
        );

    }

    my $adp = $dbadp->get_object_adaptor("Bio::SeqI");

    $adp->dbh->{'LongTruncOk'} = 0;
    $adp->dbh->{'LongReadLen'} = 1000000000;

    my $query = Bio::DB::Query::BioQuery->new();
    $query->datacollections( [ "Bio::PrimarySeqI s", ] );

    $query->where( ["s.display_id like '$locus_tag%'"] );

    my $res = $adp->find_by_query($query);

  #GENE: while ( my $seq = $res->next_object() ) {
  while ( my $seq = $res->next_object() ) {

      my $gene_name = $seq->display_name();
      #print $gene_name, "\n";

      my @feat = $seq->get_SeqFeatures();
GENE:      foreach my $f (@feat) {

	my $display_name = $f->display_name();
	#print STDERR $display_name," ", $f->primary_tag,"\n";

	next GENE if $f->primary_tag ne 'gene';
	next GENE if $f->has_tag('Dead');

	my $new_display_name;

        if (
            ($display_name =~ /^(\w+\d+)\.(\w+)\.(\d+)$/) 
            ||
            ($display_name =~ /^(\w+\d+)\.(\w+)\.(\d+)\.InterPro\.\d+$/) 
            ||
            ($display_name =~ /^(\w+\d+\.\d+)\.(\w+)\.(\d+)$/)
            ||
            ($display_name =~ /^(\w+\d+\.\d+)\.(\w+)\.(\d+)\.InterPro\.\d+$/)
           ) {
	  my ($seq_id, $source, $number) = ($1, $2, $3);
	  $new_display_name = join('.', $seq_id, $source, 'p5_hybrid', $number, 'fasta');
        }
        else {
	  die "failed to parse '$display_name':\n - does not match expected format (seqid.predictor.sequence)";
        }

	my $ss;

	$ss = $seq->subseq( $f->start, $f->end );

	unless(defined($ss)) { 

	  die "failed to fetch sequence for '$new_display_name' from BioSQL";

	}

	my $protseq = new Bio::Seq(
				   -display_id => $new_display_name,
				   -seq        => $ss
				  );
	my $newseq;

	if ( $f->strand < 0 ) {
	  $newseq = $protseq->revcom->translate->seq;
	}
	else {
	  $newseq = $protseq->translate->seq;
	}
	my $seqobj = new Bio::Seq(
				  -display_id => $new_display_name,
				  -seq        => $newseq
				 );

    my $outputfile;
    my $outputdir = $self->outputdir;
    if( defined($outputdir)  && (-d $outputdir) )
    {
        $outputfile = $outputdir."/".$new_display_name;
    }
    else
    {
        $outputfile = $new_display_name;
    }
	my $seqout = new Bio::SeqIO(
				    #-file   => ">$new_display_name",
				    -file   => ">$outputfile",
				    -format => "fasta"
				   );

	$seqout->write_seq($seqobj);

	#print STDERR "sequence should be written out\n";

	push @basket, $new_display_name;

      }
    }

    my $protname_fh = IO::File->new();

    my $basketfile = qq{$locus_tag.proteinname.fof};
    if($self->outputdir)
    {
        $basketfile = $self->outputdir."/".$basketfile;
    }

    $self->proteinnamefile($basketfile);

    $protname_fh->open(">$basketfile")
      or die "Can't open '$basketfile': $OS_ERROR";

    foreach my $name (@basket) {

      print $protname_fh qq{$name\n};

    }

    $protname_fh->close;




    return 1;

}
