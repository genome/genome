package Genome::Model::Tools::Vcf::Convert::Snv::Strelka;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Snv::Strelka {
    is => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from strelka output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from strelka snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "Strelka";
}

#Override most of the print header code so that Strelka's header is concatenated to the TGI default header
sub print_header{
    my $self = shift;
    my $public_reference = $self->_get_public_ref();

    my $input_fh = $self->_input_fh;
    my $output_fh = $self->_output_fh;

    my @header_columns = $self->_get_header_columns();
    _print_header($public_reference, $input_fh, $output_fh, \@header_columns);
    return 1;
}

sub _print_header {
    my ($public_reference, $input_fh, $output_fh, $header_columns) = @_;
    my @header_columns = @{$header_columns};

    $output_fh->print("##reference=$public_reference" . "\n");
    $output_fh->print("##phasing=none" . "\n");
    $output_fh->print("##Original Strelka header follows:" . "\n");

    while(my $line = <$input_fh>) {
      chomp($line);
      if ($line =~ /^\#\#/){
        # catch issue with Strelka DP format field having different description to other Vcfs in DV2
        # incompatible_line = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">';
        # FIXME Once joinx has ability to ignore description fields we should remove this and start using that.
        if ($line =~ /FORMAT..ID.DP,/) {
            my $compatible_line = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">';
            $output_fh->print($compatible_line, "\n");
            next;
        }
        $output_fh->print($line, "\n");
      }elsif($line =~ /^\#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+NORMAL\s+TUMOR/){
        last;
      }else{
        die "Bad header: $line";
      }
    }

    #column header:
    $output_fh->print( "#" . join("\t", @header_columns) . "\n");
    return 1;
}


#Override the entire conversion process so that Strelka VCF lines are passed through unaltered
# Loop through each input line, parse it, and print it to output
sub convert_file {
  my $self = shift;
  my $input_fh = $self->_input_fh;

  #Skip comments and check header line and die if we find a problem
  while(my $line = <$input_fh>) {
    chomp($line);
    if ($line =~ /^\#\#/){
      next;
    }elsif($line =~ /^\#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+NORMAL\s+TUMOR/){
      last;
    }else{
      last;
    }
  }

  #Simply send the data line as is without any parsing
  while(my $line = $self->get_record($input_fh)) {
    chomp $line;
    $self->write_line($line);
  }

  return 1;
}


#EXAMPLE VCF FILE PUT OUT BY STRELKA NATIVELY
##fileformat=VCFv4.1
##fileDate=20120710
##source=strelka
##startTime=Tue Jul 10 19:35:13 2012
##content=strelka somatic snv calls
##germlineSnvTheta=0.001
##priorSomaticSnvRate=1e-06
##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
##FILTER=<ID=DP,Description="Greater than 3.0x chromosomal mean depth in Normal sample">
##FILTER=<ID=BCNoise,Description="Fraction of basecalls filtered at this site in either sample is at or above 0.4">
##FILTER=<ID=SpanDel,Description="Fraction of reads crossing site with spanning deletions in either sample exceeeds 0.75">
##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
##cmdline=/gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/strelka/scripts/consolidateResults.pl --config=/gscmnt/gc2142/techd/analysis/strelka/all_chrs/results/config/run.config.ini
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
#1       10231   .       C       A       .       QSS_ref NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1      84:6:69:0:7,21:71,192:0,0:0,1
#1       10333   .       C       T       .       QSS_ref NT=ref;QSS=4;QSS_NT=4;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    13:1:0:0:0,0:12,33:0,0:0,1      49:8:2:0:0,0:37,92:0,0:4,10
#1       10440   .       C       A       .       QSS_ref NT=ref;QSS=5;QSS_NT=5;SGT=CC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    16:1:27:0:0,1:15,56:0,0:0,0     107:2:69:0:11,18:94,204:0,0:0,0
#1       10473   .       G       A       .       QSS_ref NT=ref;QSS=6;QSS_NT=6;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    9:1:0:0:0,0:0,1:8,13:0,0        62:1:0:0:18,23:0,2:43,52:0,0


#EXAMPLE VCF RESULT AFTER CONVERTING FROM SOMATIC SNIPER
##fileformat=VCFv4.1
##fileDate=20120530
##source=Sniper
##reference=ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/special_requests/assembly_variants/NCBI36_WUGSC_variant.fa.gz
##phasing=none
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth corresponding to alternate alleles 1/2/3... after software and quality filtering">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=BQ,Number=A,Type=Integer,Description="Average Base Quality corresponding to alternate alleles 1/2/3... after software and quality filtering">
##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Somatic status relative to normal counterpart: 0(wildtype), 1(germline), 2(somatic), 3(loh), 4(unknown)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Average Mapping Quality">
##FORMAT=<ID=FA,Number=1,Type=Float,Description="Fraction of reads supporting ALT">
##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CONTROL_SAMPLE_123	TUMOR_SAMPLE_123
#1	10002	.	A	C	.	PASS	.	GT:GQ:DP:BQ:MQ:AD:FA:VAQ	.:.:21:.:.:.:.:.	0/1:63:50:23:63:.:.:60
#1	10003	.	A	T	.	PASS	.	GT:GQ:DP:BQ:MQ:AD:FA:VAQ	.:.:25:.:.:.:.:.	0/1:33:57:22:33:.:.:30
#1	10248	.	A	T	.	PASS	.	GT:GQ:DP:BQ:MQ:AD:FA:VAQ	.:.:69:.:.:.:.:.	0/1:36:116:26:36:.:.:23
#1	10469	.	C	A	.	PASS	.	GT:GQ:DP:BQ:MQ:AD:FA:VAQ	.:.:16:.:.:.:.:.	0/1:19:34:34:19:.:.:16


#DIRECT COMPARISON OF DATA LINES ONLY
#CHROM  POS     ID  REF  ALT  QUAL  FILTER   INFO                                                            FORMAT                          NORMAL                         TUMOR                                  
#1       10231   .   C    A    .     QSS_ref  NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1      84:6:69:0:7,21:71,192:0,0:0,1 #Strelka
#1	      10002   .   A    C    .     PASS     .	                                                             GT:GQ:DP:BQ:MQ:AD:FA:VAQ	       .:.:21:.:.:.:.:.	               0/1:63:50:23:63:.:.:60        #Sniper



