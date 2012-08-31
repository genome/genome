package Genome::Model::Tools::Bmr::CalculateFdr;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::Bmr::CalculateFdr {
  is => 'Genome::Command::Base',
  has_input => [
  working_dir => {
    is => 'String',
    is_optional => 0,
    doc => 'A valid working directory for BMR calculations and SMG test',
  },
  ]
};

sub help_brief {
  "Combines parallel SMG tests and calculates FDR."
}

sub help_detail {
  "Combines the files containing p-values per gene, and then calculates FDR for each using R."
}

sub execute {
  my $self = shift;
  my $rlibrary = "SMG_test.R";

  my $work_dir = $self->working_dir;
  $work_dir =~ s/\/$//;
  my $genesum_dir = $work_dir . "/gene_summary_results";
  my $smgtest_dir = $work_dir . "/smg_test_results";
  my $class_file = $work_dir . "/combined.class_summary";
  my $all_genesums = $work_dir . "/combined.gene_summary";
  my $all_pvalues = $work_dir . "/combined.pvalues";
  my $all_fdrs = $work_dir . "/combined.fdr";
  my $final_results = $work_dir . "/combined.final";

  #Combine the separate gene_summary results
  opendir( GSUM_DIR, $genesum_dir );
  my @files = readdir( GSUM_DIR );
  closedir( GSUM_DIR );
  @files = grep { /\.gene_summary/ } @files;
  @files = map {$_ = "$genesum_dir/" . $_ } @files;
  my $combineGSumFh = IO::File->new( $all_genesums, ">" );
  print $combineGSumFh "Gene\tClass\tBases_Covered\tNon_Syn_Mutations\tBMR\n";
  foreach my $file ( @files )
  {
    my $fh = IO::File->new( $file );
    $fh->getline; #discard the header
    while( my $line = $fh->getline )
    {
      print $combineGSumFh $line;
    }
    $fh->close;
  }
  $combineGSumFh->close;

  #Fetch the names of the separate SMG test outputs
  opendir( SMG_DIR, $smgtest_dir );
  @files = readdir( SMG_DIR );
  closedir( SMG_DIR );
  @files = grep { /\.pvalues/ } @files;
  @files = map {$_ = "$smgtest_dir/" . $_ } @files;

  #Combine the separate SMG test outputs
  my $combineSmgFh = IO::File->new( "$all_pvalues\_full", ">" );
  $combineSmgFh->print( "Gene\tClass\tBases_Covered\tNon_Syn_Mutations\tBMR\tp\tlh0\tlh1\ttn\ttx\tp.fisher\tp.lr\tp.convol\tqc\n" );
  my %mutsPerGene = (); #This will store the total mutations per gene
  foreach my $file ( @files )
  {
    my $fh = IO::File->new( $file );
    $fh->getline; #discard the header
    while( my $line = $fh->getline )
    {
      my @segs = split( /\t/, $line );
      my ( $gene, $muts ) = ( $segs[0], $segs[3] );
      $mutsPerGene{$gene} += $muts;
      $combineSmgFh->print( $line );
    }
    $fh->close;
  }
  $combineSmgFh->close;

  $combineSmgFh = IO::File->new( "$all_pvalues\_full" );
  $combineSmgFh->getline; #Discard header
  my @pvalueLines = ();
  push( @pvalueLines, "Gene\tp.fisher\tp.lr\tp.convol\n" );
  my %skipLine = ();
  while( my $line = $combineSmgFh->getline )
  {
    my @segs = split( "\t", $line );
    next if( defined $skipLine{$segs[0]} ); #Grab only 1 line per gene
    $skipLine{$segs[0]} = 1;
    if( $mutsPerGene{$segs[0]} > 0 )
    {
      $line = join( "\t", $segs[0], $segs[10], $segs[11], $segs[12] ) . "\n";
      push( @pvalueLines, $line );
    }
    #If the total number of non-syn mutations in a gene is zero, make its pvalues = 1
    else
    {
      $line = join( "\t", $segs[0], '1', '1', '1' ) . "\n";
      push( @pvalueLines, $line );
    }
  }
  $combineSmgFh->close;

  #Write these new pvalues for input to the FDR calculations
  $combineSmgFh = IO::File->new( $all_pvalues, ">" );
  $combineSmgFh->print( @pvalueLines );
  $combineSmgFh->close;

  #Call R for FDR calculation
  my $fdr_calc_cmd = "smg_fdr(in.file='$all_pvalues',fdr.file='$all_fdrs');";
  my $fdr_calc_rcall = Genome::Model::Tools::R::CallR->create(command=>$fdr_calc_cmd,library=>$rlibrary);
  $fdr_calc_rcall->execute;

  #Gather FDR stats
  my $fdrfh = IO::File->new( $all_fdrs );
  my %FDR;
  my $fdr = \%FDR;
  $fdrfh->getline; #Discard header
  while( my $line = $fdrfh->getline )
  {
    next if $line =~ /gene/i;
    chomp $line;
    my ( $gene, $pfisher, $plr, $pconvol, $fdrfisher, $fdrlr, $fdrconvol ) = split( /\t/, $line );
    $FDR{$gene}{'pfisher'} = $pfisher;
    $FDR{$gene}{'plr'} = $plr;
    $FDR{$gene}{'pconvol'} = $pconvol;
    $FDR{$gene}{'fdrfisher'} = $fdrfisher;
    $FDR{$gene}{'fdrlr'} = $fdrlr;
    $FDR{$gene}{'fdrconvol'} = $fdrconvol;
  }
  $fdrfh->close;

  #Grab the classes from the bmr file
  my @classes;
  my $classesref = \@classes;
  my $classfh = IO::File->new( $class_file );
  while( my $line = $classfh->getline )
  {
    next if $line =~ /class/i;
    my ( $class ) = split( /\t/, $line );
    push( @classes, $class );
  }
  $classfh->close;
  @classes = sort @classes;

  #Print final output file with pvalues and FDRs per gene
  my $outfh = IO::File->new( $final_results, ">" );
  print $outfh "Gene\tTotal_Muts\t";
  print $outfh "$_\_muts\t" foreach( @classes );
  print $outfh "$_\_cov\t" foreach( @classes );
  print $outfh "p.fisher\tp.lr\tp.convol\tfdr.fisher\tfdr.lr\tfdr.convol\n";

  #Loop through BMR file and gather and print mutation information
  my %COVMUTS;
  my $covmuts = \%COVMUTS;

  my $genefh = IO::File->new( $all_genesums );
  $genefh->getline; #Discard header
  while( my $line = $genefh->getline )
  {
    my ( $nextgene, $class, $cov, $muts ) = split( /\t/, $line );
    if( !defined $COVMUTS{'gene'} )
    {
      $COVMUTS{'gene'} = $nextgene;
    }
    if ($nextgene eq $COVMUTS{'gene'})
    {
      $COVMUTS{$class}{'cov'} = $cov;
      $COVMUTS{$class}{'muts'} = $muts;
    }
    if ($nextgene ne $COVMUTS{'gene'})
    {
      #we are at next gene in file, so print output from the last gene
      my $gene2print = $COVMUTS{'gene'};
      $self->print_gene( $gene2print, $classesref, $covmuts, $fdr, $outfh );
      undef( %COVMUTS );
      $COVMUTS{'gene'} = $nextgene;
      $COVMUTS{$class}{'cov'} = $cov;
      $COVMUTS{$class}{'muts'} = $muts;
    }
  }
  #Print stuff from the last gene in the file
  $self->print_gene( $COVMUTS{'gene'}, $classesref, $covmuts, $fdr, $outfh );
  $genefh->close;
  $outfh->close;

  return 1;
}

sub print_gene {
  my ($self,$gene2print,$classesref,$covmuts,$fdr,$outfh) = @_;
  my $total_muts;

  for my $class (@$classesref)
  {
    next if( $class =~ m/gene/i );
    $total_muts += $covmuts->{$class}->{'muts'} if( defined $covmuts->{$class}->{'muts'} );
  }

  #print mutation info
  print $outfh "$covmuts->{'gene'}\t$total_muts\t";
  for my $class (@$classesref)
  {
    if( defined $covmuts->{$class}->{'muts'} )
    {
      print $outfh "$covmuts->{$class}->{'muts'}\t";
    }
    else
    {
      print STDERR "Number of variants undefined for $class in $gene2print\n";
    }
  }

  #print coverage info
  for my $class (@$classesref)
  {
    if( defined $covmuts->{$class}->{'cov'} )
    {
      print $outfh "$covmuts->{$class}->{'cov'}\t";
    }
    else
    {
      print STDERR "Amount of coverage undefined for $class in $gene2print\n";
    }
  }

  #print pvalue and fdr info. If undefined, then it had no coverage
  if( defined $fdr->{$gene2print}->{'pfisher'} )
  {
    print $outfh "$fdr->{$gene2print}->{'pfisher'}\t$fdr->{$gene2print}->{'plr'}\t$fdr->{$gene2print}->{'pconvol'}\t";
    print $outfh "$fdr->{$gene2print}->{'fdrfisher'}\t$fdr->{$gene2print}->{'fdrlr'}\t$fdr->{$gene2print}->{'fdrconvol'}\n";
  }
  else
  {
    print $outfh "NC\tNC\tNC\tNC\tNC\tNC\n";
  }
}

1;
