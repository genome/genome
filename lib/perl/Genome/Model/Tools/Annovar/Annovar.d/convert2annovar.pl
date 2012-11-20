#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: 508 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $';

our ($verbose, $help, $man);
our ($variantfile);
our ($outfile, $format, $includeinfo, $snpqual, $snppvalue, $coverage, $maxcoverage, $chr, $chrmt, $altcov, $allelicfrac, $fraction, $species, 
	$filterword, $confraction, $allallele, $withzyg);

our %iupac = (R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-'); ### <<< FOR 5500SOLiD LifeScope ( S=>'GC' is replaced by S=>'CG')
our %iupacrev = reverse %iupac; ### <<< FOR 5500SOLiD LifeScope

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'includeinfo'=>\$includeinfo,
	'snpqual=f'=>\$snpqual, 'snppvalue=f'=>\$snppvalue, 'coverage=i'=>\$coverage, 'maxcoverage=i'=>\$maxcoverage, 'chr=s'=>\$chr, 'chrmt=s'=>\$chrmt, 
	'fraction=f'=>\$fraction, 'altcov=i'=>\$altcov, 'allelicfrac'=>\$allelicfrac,
	'species'=>\$species, 'filter=s'=>\$filterword, 'confraction=f'=>\$confraction, 'allallele!'=>\$allallele, 'withzyg'=>\$withzyg) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($variantfile) = @ARGV;

$chrmt ||= 'M';

if (not $format) {
	$format = 'pileup';
	print STDERR "NOTICE: the default --format argument is set as 'pileup'\n";
}

if (defined $outfile) {
	open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
}

defined $snpqual and $format eq 'pileup' || $format eq 'vcf4' || pod2usage ("Error in argument: the --snpqual is supported only for the 'pileup' or 'vcf4' format");
defined $snppvalue and $format eq 'gff3-solid' || pod2usage ("Error in argument: the --snppvalue is supported only for the 'gff3-solid' format");
if (not defined $snpqual and $format eq 'pileup') {
	$snpqual = 20;
	print STDERR "NOTICE: the default --snpqual argument for pileup format is set as 20\n";
}

if (not defined $snppvalue) {
	$snppvalue = 1;		#by default, do not use any of the P-value cutoff in filtering out GFF3-SOLID files (this is differnt from handling pileup files)
}

if (not defined $coverage) {
	$coverage = 0;
}

if (defined $fraction) {
	$format eq 'pileup' or $format eq 'vcf4' or pod2usage ("Error in argument: the '--fraction' argument is supported for the pileup or vcf4 format only");
	$format eq 'vcf4' and print STDERR "NOTICE: the --fraction argument works ONLY on indels for vcf4 format\n";
	$fraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --fraction argument must be between 0 and 1 inclusive");
} else {
	$fraction = 0;
}

if (defined $confraction) {
	$format eq 'vcf4' and print STDERR "NOTICE: the --confraction argument works ONLY on indels for vcf4 format\n";
	$confraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --confraction argument must be between 0 and 1 inclusive");
} else {
	$confraction = 0;
}

if (defined $altcov) {
	$format eq 'pileup' or pod2usage ("Error in argument: the '--altcov' argument is supported for the '--format pileup' only");
	$altcov < $coverage or pod2suage ("Error in argument: the --altcov argument must be less than --coverage");
	$altcov > 0 or pod2suage ("Error in argument: the --altcov argument must be a positive integer");
}

if (defined $species) {
	$format eq 'gff3-solid' or pod2usage ("Error in argument: the '--species' argument is only necessary for the '--format gff3-solid'");
}

if ($allallele) {
	$format eq 'vcf4' or pod2usage ("Error in argument: the '--allallele' argument is only supported for the '--format vcf4'");
}

if ($format eq 'pileup') {
	convertPileup ($variantfile);
} elsif ($format eq 'cg') {
	convertCG ($variantfile);
} elsif ($format eq 'cgmastervar') {
	convertCGMasterVar ($variantfile);
} elsif ($format eq 'gff3-solid') {
	convertGFF3SolidSNP ($variantfile);
} elsif ($format eq 'soap') {
	print STDERR "WARNING: the support for '--format soap' is not well developed yet and may contain bugs for indel analysis.\n";
	convertSOAP ($variantfile);
} elsif ($format eq 'maq') {
	print STDERR "WARNING: the support for '--format maq' is not well developed yet and may contain bugs.\n";
	convertMAQSNP ($variantfile);
} elsif ($format eq 'casava') {
	if (not defined $chr) {
		pod2usage ("Error in argument: please supply --chr argument for the '--format casava'");
	}
	convertCASAVA ($variantfile, $chr);
} elsif ($format eq 'vcf4') {
	convertVCF4 ($variantfile);
} elsif ($format eq 'annovar') {
	convertANNOVAR ($variantfile);
} else {
	pod2usage ("Error in argument: the --format $format is not currently supported. Please contact ANNOVAR developer for adding the support");
}


sub convertPileup {
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0/;
	
	if ($variantfile eq 'stdin') {
		*VAR = *STDIN;
	} else {
		open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	}
	print STDERR "NOTICE: Column 6-9 in output are heterozygosity status, SNP quality, total reads, reads with mutation\n";

	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my $hom = 'hom';
		my @field = split (/\t/, $_);
		@field >= 10 or die "Error: invalid record found in pileupfile $variantfile (at least 10 fields expected): <$_>\n";
		my ($chr, $pos, $wt, $call, @other) = @field;
		my ($cons_qual, $snp_quality, $readcount, $readallele) = @field[4,5,7,8];
		$chr =~ s/^chr//;
		$wt = uc $wt;					#wt may or may not in upper case, it depends on the input FASTA file
		$call = uc $call;				#usually call is in upper case
		$readallele = uc $readallele;			#lower case indicate the opposite strand
		
		$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
		
		$snp_quality >= $snpqual or next;		#quality of the variant call
		$readcount >= $coverage or next;		#read count of the variant call
		$maxcoverage and $readcount <= $maxcoverage || next;	#maximum coverage of the variant call
		
		if ($wt eq '*') {				#indel
			#example:
			#1       970271  *       +C/+C   39      106     44      5       +C      *       1       4       0       0       0
			#1       1548977 *       */+CCG  29      29      42      3       *       +CCG    2       1       0       0       0
			#1       1674810 *       */+C    24      24      42      6       *       +C      5       1       0       0       0
			#1       968466  *       -CT/-CT 53      339     55      5       -CT     *       5       0       0       0       0
			#1       1093600 *       -GAAA/* 29      29      53      3       -GAAA   *       1       2       0       0       0
			#1       1110101 *       */-A    41      41      17      6       *       -A      5       1       0       0       0
			#1       1215395 *       */-TC   26      26      32      4       *       -TC     3       1       0       0       0
			my @obs = split (/\//, $call);		#example: '+AG/+AG' as homozygotes, '*/-TA' or '*/+T' as heterozygotes
			@obs == 2 or die "Error: pileup record contains more than two alternative alleles: <$_>\n";
			my ($end, $ref, $alt);
			my ($indelreadcount);			#number of reads supporting the indel
			
			
			if ($obs[0] eq $obs[1]) {
				#something weird may occur in SamTools output: 22      15231121        *       */*     360     0       32      158     *       +GA     156     2       0       0       0
				$obs[0] eq '*' and next;	
	
				#for deletions, SAMTOOLS represent deletion using a location before the first deleted base in the reference sequence coordinate system
				#for example, a deletion in Samtools output is "1       109266688       *       */-CTT  1429    1429    58      43      *       -CTT    24      19      0       0       0"
				#the correct ANNOVAR input (for rs35029887) should be "1       109266689       109266691       CTT     -       het     1429"
				#insertions are fine without change; for example, check rs5745466 in Genome Browser; samtools report "1       76119508        *       +AT/+AT"
				#for this insertion, ANNOVAR input file (for rs5745466) becomes "1       76119508        76119508        -       AT      hom     1601"

				if ($obs[0] =~ m/^\-/) {
					$pos++;			#add 1 to position in deletion
				}
				
				$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
				$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $indelreadcount >= $altcov || next;
				
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = sprintf ("%.3f", $indelreadcount/$readcount);
				}
				($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
				print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
				$counthom++;
			} else {
				$hom = 'het';
				if ($obs[0] =~ m/^[\-\+]/) {
					$obs[0] =~ m/^\-/ and $pos++;
					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
					$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
					defined $altcov and $indelreadcount >= $altcov || next;
					
					if ($chr eq $chrmt or $allelicfrac) {
						$hom = sprintf ("%.3f", $indelreadcount/$readcount);
					}
					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
					$counthet++;
				}
				if ($obs[1] =~ m/^[\-\+]/) {
					$obs[1] =~ m/^\-/ and $pos++;
					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[1]);
					$indelreadcount = calculateIndelReadCount ($obs[1], \@field);
					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
					defined $altcov and $indelreadcount >= $altcov || next;
					
					if ($chr eq $chrmt or $allelicfrac) {
						$hom = sprintf ("%.3f", $indelreadcount/$readcount);
					}
					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
					$counthet++;
				}
			}
			$countindel++;
		} else {
			#1       798494  G       A       36      36      58      3       AAA     bbb
			#1       798677  T       K       33      33      52      26      ,$.,,G.GG,.,......,..G,,...     b`bbBaJIbFbZWaTNQbb_VZcbbb
			#1       856182  G       A       42      42      50      5       AAAAA   B\bbb
			#1       861034  A       M       48      48      49      14      ,$,.,..,cc.c.,c bbBbb`]BFbHbBB
			#1       864289  T       K       22      22      56      6       .g,,g,  BH^_BB
			
			$wt eq $call and next;			#this is not a SNP
			my $obs = $iupac{$call} or die "Error: invalid best call ($call) in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] ne $obs[1]) {
				$hom = 'het';
			}
				
			
			if ($obs[0] eq $wt) {			#obs[0] is guaranteed to be an alternative allele
				@obs = @obs[1,0];
			}
			if ($wt eq 'A' and $obs[0] eq 'G' or $wt eq 'G' and $obs[0] eq 'A' or $wt eq 'C' and $obs[0] eq 'T' or $wt eq 'T' and $obs[0] eq 'C') {
				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
					$countti++;
				}
				
			} else {
				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
					$counttv++;
				}
			}
			
			my $mutallelecount;
			
			if ($obs[1] eq $wt) {			#het SNP
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthet++;
			} elsif ($obs[1] ne $obs[0]) {		#het SNP but both differ from reference allele
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[1], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[1], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[1], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				if ($chr eq $chrmt) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthet++;
				$counthet++;
			} else {				#homo SNP
				if ($chr eq $chrmt or $allelicfrac) {
					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
				}
				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
				
				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
				$counthom++;
			}
			$countsnp++;
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
}

sub calculateIndelReadCount {
	my ($obs, $field) = @_;
	#make sure to use upper case in the comparison, for example:
	#chr10   130133  *       */-ca   189     189     59      31      *       -ca     27      4       0       0       0
	if ($obs eq uc $field->[8]) {
		return $field->[10];
	} elsif ($obs eq uc $field->[9]) {
		return $field->[11];
	} else {
		die "Error: invalid record in pileup file (indel counts cannot be inferred): <$obs> vs <@$field>\n";
	}
}

sub calculateMutAlleleCount {
	my ($allele, $string) = @_;	#they should be already both in upper case
	$string =~ s/\^.//g;		#^ is followed by mapping quality
	$string =~ s/\$//g;
	$string =~ s/[+-]1[^\d]//g;	#1 followed by a non-number
	$string =~ s/[+-]2..//g;
	$string =~ s/[+-]3...//g;
	$string =~ s/[+-]4....//g;
	$string =~ s/[+-]5.....//g;
	$string =~ s/[+-]6......//g;
	$string =~ s/[+-]7.......//g;
	$string =~ s/[+-]8........//g;
	$string =~ s/[+-]9.........//g;
	$string =~ s/[+-]10..........//g;
	
	#make sure to use upper case letter
	my @string = split (//, uc $string);
	my $count = 0;
	for my $i (0 .. @string-1) {
		$allele eq $string[$i] and $count++;
	}
	return $count;
}

sub calculateAllelicFraction {
	my ($obs, $readbase, $readcount) = @_;
	my @readbase = split (//, $readbase);
	my $count=0;
	for my $i (0 .. @readbase-1) {
		uc $obs eq uc $readbase[$i] and $count++;
	}
	my $hom = $count/$readcount;
	length ($hom) > 5 and $hom > 0.001 and $hom = sprintf ("%.3f", $hom);
	return $hom;
}

sub recalculateEndRefObs {		#recalculate end position, reference allele and observed allele
	my ($end, $ref, $obs) = @_;
	if ($obs =~ m/^\-(\w+)/) {	#deletion
		$end += (length ($1)-1);
		$ref = $1;
		$obs = '-';
	} elsif ($obs =~ m/^\+(\w+)/) {	#insertion
		$ref = '-';
		$obs = $1;
	} else {
		die "Error: cannot handle $end, $ref, $obs\n";
	}
	return ($end, $ref, $obs);
}

sub convertCG {
	my ($variantfile) = @_;
	
	my ($foundheader, $countline, @field);
	my ($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = qw/0 0 0 0 0 0 0 0/;
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	print STDERR "NOTICE: Converting variants from $variantfile\n";
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		if (m/^>locus/) {
			$foundheader++;
		}
		if (not $foundheader) {
			$countline > 50 and die "Error: invalid CG-var file format for $variantfile (>locus record is not found within the first 50 lines)\n";
			next;
		}
		my ($locus, $ploidy, $haplo, $chr, $start, $end, $vartype, $ref, $obs, $score, $haplolink, $xref) = split (/\t/, $_);
		$chr =~ s/^chr//;
		$vartype eq 'ins' or $start++;		#CG use zero-start, half-open coordinate. Insertion does not need to be processed (example, "16476   2       2       chr1    751820  751820  ins             T       49              dbsnp:rs59038458")
		$obs eq '' and $obs = '-';
		$ref eq '' and $ref = '-';

		if ($vartype =~ m/^snp|ins|del|delins|sub$/) {		#in new versions of the files, "sub" is used instead of "delins".
			#$chr eq 'M' and next;			#ignore chrM markers as they are not diploid
			if ($chr eq $prechr and $start eq $prestart and $end eq $preend and $obs eq $preobs) {		#homozygous mutation
				print $chr, "\t", $start, "\t", $end, "\t", $ref, "\t", $obs, "\t", $vartype, "\t", ($score+$prescore)/2, "\t", "hom\t", $xref, "\n";
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = qw/0 0 0 0 0 0 0 0/;
			} else {
				if ($prestart and $preend) {
					print $prechr, "\t", $prestart, "\t", $preend, "\t", $preref, "\t", $preobs, "\t", $prevartype, "\t", $prescore, "\thet\t", $prexref, "\n";
				}
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = ($chr, $start, $end, $vartype, $ref, $obs, $score, $xref);
			}
		}
	}
	if ($prestart and $preend) {
		print $prechr, "\t", $prestart, "\t", $preend, "\t", $preref, "\t", $preobs, "\t", $prevartype, "\t", $prescore, "\thet\t", $prexref, "\n";
	}
	print STDERR "NOTICE: Done with $countline lines\n";
}

sub convertCGMasterVar {
	#this subroutine converts CG masterVar format into ANNOVAR input format
	#example input file is below:
	
	#SEGDUP_GENERATED_AT    2010-Dec-01 13:40
	#SOFTWARE_VERSION       2.0.2.26
	#TYPE   VAR-OLPL
	#>locus  ploidy  chromosome      begin   end     zygosity        varType reference       allele1Seq      allele2Seq      allele1VarScoreVAF      allele2VarScoreVAF      allele1VarScoreEAF      allele2VarScoreEAF      allele1VarQuality       allele2VarQuality       allele1HapLink  allele2HapLink  allele1XRef     allele2XRef     evidenceIntervalId      allele1ReadCount        allele2ReadCount        referenceAlleleReadCount        totalReadCount  allele1Gene     allele2Gene     pfam    miRBaseId       repeatMasker    segDupOverlap   relativeCoverageDiploid calledPloidy    relativeCoverageNondiploid      calledLevel
	#1       2       chr1    0       10000   no-call no-ref  =       ?       ?                                                                                                                                                                                                       
	#2       2       chr1    10000   11038   no-call complex =       ?       ?                                                                                                                                                                               1.13    N       1.02    1.006
	#3       2       chr1    11038   11055   hom     ref     =       =       =                                                                                                                                                                               1.13    N       1.02    1.006
	#4       2       chr1    11055   11082   no-call complex =       ?       ?                                                                                                                                                                               1.13    N       1.02    1.006
	#5       2       chr1    11082   11109   hom     ref     =       =       =                                                                                                                                                                               1.13    N       1.02    1.006

	my ($variantfile) = @_;
	
	my ($foundheader, $countline);
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	print STDERR "NOTICE: Converting variants from $variantfile\n";
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		if (m/^>locus/) {
			$foundheader++;
		}
		if (not $foundheader) {
			$countline > 50 and die "Error: invalid CG-var file format for $variantfile (>locus record is not found within the first 50 lines)\n";
			next;
		}
		my ($locus, $ploidy, $chr, $start, $end, $zygosity, $vartype, $ref, $obs1, $obs2, @otherinfo) = split (/\t/, $_, -1);
		my ($a1scorevaf, $a2scorevaf, $a1scoreeaf, $a2scoreeaf, $a1varqual, $a2varqual, $a1haplink, $a2haplink, $a1xref, $a2xref, $interval, $a1read, $a2read, $refread, $totalread) = @otherinfo;

		$chr =~ s/^chr//;
		$vartype eq 'ins' or $start++;		#CG use zero-start, half-open coordinate. Insertion does not need to be processed (example, "16476   2       2       chr1    751820  751820  ins             T       49              dbsnp:rs59038458")
		$obs1 eq '' and $obs1 = '-';
		$obs2 eq '' and $obs2 = '-';
		$ref eq '' and $ref = '-';

		#zygosity explanation:
		##no-call: All alleles are partially or fully no-called.
		##hap: Haploid, fully called locus.
		##half: Diploid locus where one of the alleles is fully called and the other contains no-calls.
		##hom: Diploid, homozygous, fully called locus.
		##het-ref: Diploid, heterozygous, fully called locus where one of the alleles is identical to the reference.
		##het-alt: Diploid, heterozygous, fully called locus where both alleles differ from the reference.

		#vartype explanation:
		##snp, ins, del, or sub: Fully called or half-called locus that contains only a single isolated variation.
		##ref: Fully called or half-called locus that contains only reference calls and no calls and at least one allele is fully called.
		##complex: Locus that contains multiple variations or has no-calls in all alleles. This is also the value for all loci where the reference itself is ambiguous.
		##no-ref: Locus where the reference genome is N.
		##PAR-called-in-X: Locus on the pseudo-autosomal region of the Y chromosomes in males.
                
                $zygosity eq 'no-call' and next;		#ignore locus without calls
                

                
                
		if ($vartype =~ m/^snp|ins|del|delins|sub$/) {		#in new versions of the files, "sub" is used instead of "delins".
	                if ($coverage) {
	                	$totalread >= $coverage or next;
	                }
	                if ($maxcoverage) {
	                	$totalread <= $maxcoverage or next;
	                }
			if ($zygosity eq 'hom' or $zygosity eq 'hap') {
				print $chr, "\t", $start, "\t", $end, "\t", $ref, "\t", $obs1, "\t", 'hom', "\t", $vartype, "\t", $totalread, "\t", $includeinfo?join("\t", "\t", @otherinfo):'', "\n";
			} else {
				print $chr, "\t", $start, "\t", $end, "\t", $ref, "\t", $obs1, "\t", 'het', "\t", $vartype, "\t", $totalread, "\t", $includeinfo?join("\t", "\t", @otherinfo):'', "\n";    
			}
		}
	}
	print STDERR "NOTICE: Done with $countline lines\n";
}

sub convertGFF3SolidSNP {
	my ($variantfile) = @_;
	my ($countline, $countvar, $countallvar, @other) = (0, 0, 0);
	my ($unknown_count);		#count of variants with 'unknown' variation type
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	$_ = <VAR>;
	s/[\r\n]+$//;
	m/^##gff-version\s+3/ or die "Error: invalid first line in GFF3 file ('##gff-version 3' expected): <$_>\n";
	$_ = <VAR>;
	s/[\r\n]+$//;
	(m/^##solid-gff-version/ || m/^##source-version/) or print STDERR "WARNING: problematic second line in GFF3-SOLiD file ('##solid-gff-version' or '##source-version' expected): <$_>\n"; ### <<< FOR 5500SOLiD LifeScope

	print STDERR "NOTICE: Column 6-9 in output are heterozygosity status, variant score (P-value), total clipped normal coverage reads, total reads with mutated allele\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		m/^##/ and next;		#header of comment lines
		m/^#/ and next;			#header of results lines
		
		my @field = split (/\t/, $_);
		@field == 9 or die "Error: invalid record found in $variantfile (10 fields expected): <$_>\n";
		my ($chr, $program, $type, $pos, $end, $score, $attribute) = @field[0,1,2,3,4,5,8];		#score is the P-value for the SNP calls
		$chr eq 'chr_name' and next;	#header line
		
		if ($score ne '.') {
			$score >=0 and $score <=1 or die "Error: invalid score record found in file (0-1 range expected): <$_>\n";
			$score <= $snppvalue or next;
		}
		
		if ($species and $species eq 'human') {
			$chr eq '23' and $chr = 'X';
			$chr eq '24' and $chr = 'Y';
			$chr eq '25' and $chr = 'M';
		}

		$includeinfo and @other = ($attribute);			#unless -includeinfo is set, the other will not be printed

		my ($readcount, $mutallelecount) = ('.', '.');		#total coverage, coverage for mutated alleles
		
		if ($type eq 'unknown') {
			#SOLiD GDD3 may have unknown variation types
			#chr1    AB_SOLiD Small Indel Tool       unknown 3833062 3833153 1       .       .       ID=5483;len=no_call;allele-call-pos=3833062;allele-call=/CCAC;allele-pos=3833057;alleles=atccatccacccatc/aTCCATCCACCCACCCATC/NO_CALL;allele-counts=REF,2,2;tight_chrom_pos=none;loose_chrom_pos=3833058-3833069;no_nonred_reads=3;coverage_ratio=8.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_50_10_r,L1_1_50_15_r,L1_1_50_15_r,L1_1_50_12_r;bead_ids=1018_196_970,699_1263_465,220_513_923,2022_1532_1071;overall_qvs=4,6,2,50;no_mismatches=5,4,2,0;read_pos=27,29,31,13;from_end_pos=23,21,19,37;strands=+,+,+,+;tags=R3,F3,F3,F3;indel_sizes=-92,-112,4,4;non_indel_no_mismatches=0,0,8,0;unmatched-lengths=50,50,50,50;ave-unmatched=50.0000;anchor-match-lengths=48,49,49,49;ave-anchor-length=48.7500;read_seqs=G23223321322112233223100132013201320110011001322332,T33223321322112233223100132013201320110013021322332,T33223321322112233223100132013201320110011001322332,T31001320132013201100110013223322113030332233113032;base_qvs=;non_indel_seqs=T21322332211221121322332230321212121223322332233221,G12020202202020012001200213022002130012332310122030,G12020202202020012001000210022012110312331331122030,G22111012101031010100002002321020002202121121313021;non_indel_qvs=
			$unknown_count++;
			next;		#do not count this one!
		}
		
		if ($program eq 'SOLiD_diBayes' or $program eq 'AB_SOLiD SNP caller') {		#SNP variants
			#detailed description can be found at http://solidsoftwaretools.com/gf/download/docmanfileversion/208/866/DiBayes_Documentation_v1.2.pdf
			#chr1    SOLiD_diBayes   SNP     559817  559817  0.094413        .       .       genotype=Y;reference=T;coverage=9;refAlleleCounts=5;refAlleleStarts=4;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=14;diColor1=11;diColor2=33;het=1;flag= 
			#chr1    SOLiD_diBayes   SNP     714068  714068  0.000000        .       .       genotype=M;reference=C;coverage=13;refAlleleCounts=7;refAlleleStarts=6;refAlleleMeanQV=25;novelAlleleCounts=6;novelAlleleStarts=4;novelAlleleMeanQV=22;diColor1=00;diColor2=11;het=1;flag= 
			#chr1    SOLiD_diBayes   SNP     714835  714835  0.041579        .       .       genotype=R;reference=A;coverage=5;refAlleleCounts=3;refAlleleStarts=3;refAlleleMeanQV=18;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=20;diColor1=02;diColor2=20;het=1;flag= 

			$pos == $end or die "Error: invalid record found in GFF3-SOLiD file: start and end discordant: <$_>\n";
	
			my ($wt, $call);
			my ($hit); ### <<< FOR 5500SOLiD LifeScope

			if ($attribute =~ m/ref_base=(\w)/) {
				$wt = $1;
			} elsif ($attribute =~ m/reference=(\w)/) {
				$wt = $1;
			} else {
				die "Error: invalid record found in GFF3-SOLiD file (ref_base/reference was not found): <$_>\n";
			}
			
			if ($attribute =~ m/consen_base=(\w)/) {
				$call = $1;
			} elsif ($attribute =~ m/genotype=(\w)/) {
				$call = $1;
			} elsif ($attribute =~ m/allele-call=([\w\/]+)/) { ### <<< FOR 5500SOLiD LifeScope
			        $hit = $1;
			        if ($hit =~ m/\//) {
			            $call = $iupacrev{join("",sort(split(/\//,$hit)))}; 
			        } else {
				    $call = $hit;
			        }
			} else {
				die "Error: invalid record found in GFF3-SOLiD file (consen_base was not found): <$_>\n";
			}
						
			if ($attribute =~ m/coverage=(\d+)/) {
				$readcount = $1;
				$readcount >= $coverage or next;		#read count of the variant call
				$maxcoverage and $readcount <= $maxcoverage || next;
			}
			if ($attribute =~ m/novelAlleleCounts=(\d+)/) {
				$mutallelecount = $1;
				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				defined $altcov and $mutallelecount >= $altcov || next;
			}
			
			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
				$countallvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
			}
		} elsif ($program eq 'AB_CNV_PIPELINE') {	#CNV
			if ($attribute =~ m/copynum=(\d+)/ or $attribute =~ m/copynumber=(\d+)/) {
				if ($1 < 2) {
					print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", '-', "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
				} elsif ($1 > 2) {
					print $chr, "\t", $end, "\t", $end, "\t", '-', "\t", 0, "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
				}
			} else {
				print $chr, "\t", $end, "\t", $end, "\t", '-', "\t", 0, "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
			}
		} elsif ($program eq 'AB_SOLiD Large Indel Tool') {	#CNV
			#http://solidsoftwaretools.com/gf/download/docmanfileversion/182/780/Large_Indel_Documentation_v1.0.0.pdf
			## [FIELDS] (1) chromosome (2) version (3) indel type (4) breakpoint start (5) breakpoint end (6) p-value (7) NA (8) NA (9) attributes
			#chr10   AB_SOLiD Large Indel Tool       insertion       151910  151969  2.77548e-11     .       .       dev=-71;avgDev=-1.63884;zygosity=HOMOZYGOUS;nRef=0;nDev=14;refDev=0;devDev=-1.60924;refVar=0;devVar=0.0159438;beadIds=1750_720_1641,649_1680_794,503_1756_608,1726_174_1362,1208_1772_353,872_594_1604,1924_990_858,1899_961_1848,901_1226_378,323_1750_1017,1185_180_1237,1519_490_1074,1291_94_324,285_758_922,1135_95_1594,1055_218_1279,
			#chr10   AB_SOLiD Large Indel Tool       insertion       154109  154729  2.1559e-11      .       .       dev=-66;avgDev=-1.51253;zygosity=HOMOZYGOUS;nRef=0;nDev=15;refDev=0;devDev=-1.02864;refVar=0;devVar=0.133236;beadIds=1728_1671_1739,1250_231_25,811_783_1090,1035_908_491,649_1680_794,503_1756_608,1726_174_1362,1208_1772_353,872_594_1604,1924_990_858,1899_961_1848,901_1226_378,323_1750_1017,1185_180_1237,1519_490_1074,1291_94_324,285_758_922,1135_95_1594,1055_218_1279,
			my ($call, @call, $zygosity);
			if ($attribute =~ m#zygosity=HEMIZYGOUS#) {
				$zygosity = 'het';
			} elsif ($attribute =~ m#zygosity=HOMOZYGOUS#) {
				$zygosity = 'hom';
			} else {
				$zygosity = 'unk';
			}
			if ($type eq 'insertion') {
				#the true boundary is unknown (start is different from end) so we cannot use "-" to represent reference allele.
				print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", 0, "\t", $zygosity, "\t", "$score\t.\t.\t", join ("\t", @other), "\n";
			} elsif ($type eq 'deletion') {
				print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", '-', "\t", $zygosity, "\t", "$score\t.\t.\t", join ("\t", @other), "\n";
			}
		} elsif ($program eq 'AB_SOLiD Small Indel Tool') {		#small indels
			#typical simple insertion and deletions
			#chr1    AB_SOLiD Small Indel Tool       deletion        1352612 1352614 1       .       .       ID=1290;del_len=3;allele-call-pos=1352612;allele-call=cct/;allele-pos=1352610;alleles=cccctccat/cCCCAT;allele-counts=REF,2;tight_chrom_pos=1352612-1352614;loose_chrom_pos=1352612-1352614;no_nonred_reads=2;coverage_ratio=11.5000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_3_r,L1_1_25_8_r;bead_ids=1470_2000_506,822_1710_1767;overall_qvs=18,19;no_mismatches=3,3;read_pos=6,13;from_end_pos=19,12;strands=-,+;tags=R3,R3;indel_sizes=-3,-3;non_indel_no_mismatches=1,-1;unmatched-lengths=25,25;ave-unmatched=25.0000;anchor-match-lengths=24,99;ave-anchor-length=61.5000;read_seqs=G0112310001100003120031200,G0300213000011000132110021;base_qvs=;non_indel_seqs=T2120033002022200220000002,;non_indel_qvs=
			#chr1    AB_SOLiD Small Indel Tool       insertion_site  1311162 1311162 1       .       .       ID=1249;ins_len=1;allele-call-pos=1311162;allele-call=/G;allele-pos=1311161;alleles=gaggggggg/GAGGGGGGGG/NO_CALL;allele-counts=REF,3,1;tight_chrom_pos=none;loose_chrom_pos=1311160-1311169;no_nonred_reads=3;coverage_ratio=4.6667;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_6_r,L1_1_50_10_r,L1_1_25_2_r,L1_1_25_3_r;bead_ids=850_837_429,1160_878_181,404_1050_1881,1084_64_1343;overall_qvs=20,56,25,25;no_mismatches=3,2,2,1;read_pos=11,22,11,11;from_end_pos=14,28,14,14;strands=+,-,-,-;tags=R3,F3,F3,F3;indel_sizes=1,1,1,1;non_indel_no_mismatches=-1,1,0,1;unmatched-lengths=25,50,25,25;ave-unmatched=31.2500;anchor-match-lengths=99,49,24,24;ave-anchor-length=49.0000;read_seqs=G1020001130221020000000020,T03223323210110021000000022122030100020221222222122,T0102210000000221223301000,T0102210000000221220301000;base_qvs=;non_indel_seqs=,G21202030032202013220021321131212021000122300013132,G1331133120001221220120120,G1331133120001221220120220;non_indel_qvs=
			
			#sometimes, allele-call is ambiguous that requires a "block substitution" representation (although they were annotated as insertion or deletion by SOLiD, they should be treated as block substitution by ANNOVAR)
			#sometimes, mutiple allele calls may be present at the same site
			#chr1    AB_SOLiD Small Indel Tool       deletion        1209357 1209360 1       .       .       ID=1101;del_len=4;allele-call-pos=1209357;allele-call=ggtggg/TT;allele-pos=1209355;alleles=ggggtgggggggtt/gGTTGGGGTT/gGTGTTTTGCCTT/NO_CALL;allele-counts=REF,3,1,1;tight_chrom_pos=none;loose_chrom_pos=1209357-1209363;no_nonred_reads=4;coverage_ratio=3.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=0.9888;run_names=L1_1_25_1_r,L1_1_25_2_r,L1_1_25_4_r,L1_1_25_3_r,L1_1_25_7_r;bead_ids=1017_1024_53,1493_1896_615,1794_647_1473,307_1116_687,773_1492_1671;overall_qvs=24,24,28,24,8;no_mismatches=2,3,2,3,2;read_pos=14,9,14,9,15;from_end_pos=11,16,11,16,10;strands=-,+,-,+,+;tags=F3,R3,F3,F3,F3;indel_sizes=-4,-4,-4,-4,3;non_indel_no_mismatches=1,0,0,0,0;unmatched-lengths=25,25,25,25,25;ave-unmatched=25.0000;anchor-match-lengths=24,24,24,24,24;ave-anchor-length=24.0000;read_seqs=T2221100101000101000221100,G0001122000100000101001020,T2221100101000101000221100,T1112200010100010100112000,T1011220000111000130200001;base_qvs=;non_indel_seqs=G0312033221312111022200300,T0111113210210112100001130,G0312133221312111022200300,G0231003132222112000012020,G3121331033101113122312020;non_indel_qvs=
			#chr1    AB_SOLiD Small Indel Tool       deletion        1209436 1209436 1       .       .       ID=1103;del_len=1;allele-call-pos=1209436;allele-call=ag/A/G;allele-pos=1209434;alleles=tgagggggtt/tGAGGGGTT/tGGGGGGTT;allele-counts=REF,1,1;tight_chrom_pos=none;loose_chrom_pos=1209436-1209441;no_nonred_reads=2;coverage_ratio=5.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_6_r,L1_1_25_2_r;bead_ids=1315_1584_2005,1706_194_437;overall_qvs=28,21;no_mismatches=0,3;read_pos=9,7;from_end_pos=16,18;strands=-,-;tags=R3,R3;indel_sizes=-1,-1;non_indel_no_mismatches=-1,0;unmatched-lengths=25,25;ave-unmatched=25.0000;anchor-match-lengths=99,24;ave-anchor-length=61.5000;read_seqs=G3001010000011001010000001,G3010100022110010111000110;base_qvs=;non_indel_seqs=,T1112003220020013202122300;non_indel_qvs=
			#chr1    AB_SOLiD Small Indel Tool       insertion_site  1424540 1424540 1       .       .       ID=1376;ins_len=3;allele-call-pos=1424540;allele-call=tt/CCCAC;allele-pos=1424537;alleles=ttttttg/TTTCCCACTG/NO_CALL;allele-counts=REF,1,1;tight_chrom_pos=none;loose_chrom_pos=1424536-1424543;no_nonred_reads=2;coverage_ratio=11.5000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_7_r,L1_1_50_16_r;bead_ids=703_437_370,1089_878_1744;overall_qvs=1,9;no_mismatches=3,4;read_pos=5,35;from_end_pos=20,15;strands=-,-;tags=R3,F3;indel_sizes=3,3;non_indel_no_mismatches=2,0;unmatched-lengths=25,50;ave-unmatched=37.5000;anchor-match-lengths=24,47;ave-anchor-length=35.5000;read_seqs=G2032002200200000000000020,T30100102220312202103112130230322210121100200002100;base_qvs=;non_indel_seqs=T2121120003012303000000000,G22213300221101011121030022002222300220322213303102;non_indel_qvs=
			my ($call, @call, $zygosity);
			my ($refcall, $gapnonred, %temphash); ### <<< FOR 5500SOLiD LifeScope
			#if ($attribute =~ m#experimental-zygosity=HEMIZYGOUS# ||$attribute =~ m#zygosity=HEMIZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
			#	$zygosity = 'het';
			#} elsif ($attribute =~ m#experimental-zygosity=HOMOZYGOUS# || $attribute =~ m#zygosity=HOMOZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
			#the above 3 lines are replaced by the following 3 lines on 20120618
			if ($attribute =~ m#zygosity=(MULTI-)?HEMIZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
			   $zygosity = 'het';
			} elsif ($attribute =~ m#zygosity=(MULTI-)?HOMOZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
				$zygosity = 'hom';
			} else {
				$zygosity = 'unk';
			}
			$score = '.';			#by default, score=1 in the output
			
			#no_nonred_reads: Number of reads with unique start positions (non-redundant reads).
			#coverage_ratio: Clipped normal coverage/number of non-redundant reads.Clipped coverage is where the parts of the read that are within 5 bp at either end are not counted as a part of coverage.
			if ($attribute =~ m/no_nonred_reads=(\d+);coverage_ratio=([\d\.]+)/) {
				$readcount = int ($1*$2);	
				$readcount >= $coverage or next;		#clipped normal read count of the variant call (this could even be lower than mut allele count)
				$maxcoverage and $readcount <= $maxcoverage || next;
			} elsif ($attribute =~ m/gap-nonred-reads=(\d+)/) { ### <<< FOR 5500SOLiD LifeScope
				$gapnonred = $1;
				$attribute =~ m/coverage_ratio=(\d+)/;
				$readcount = int($gapnonred*$1);
				$readcount >= $coverage or next;
				$maxcoverage and $readcount <= $maxcoverage || next;
			} else {
				$readcount = '.';
			}
			if ($attribute =~ m/allele-counts=REF,(\d+)/) {
				$mutallelecount = $1;
			} elsif ($attribute =~ m/context-variant-reads=(\d+)/) { ### <<< FOR 5500SOLiD LifeScope
			    	$mutallelecount = $1;
			}
			if ($attribute =~ m#reference=([\w\-]+)#) { ### <<< FOR 5500SOLiD LifeScope (using "reference" tag for the reference allele) 
			       	$refcall = $1;
				$attribute =~ m#;allele\-call=([\w\-\/]+)#;
				foreach my $item(split(/\//, $1)) { $temphash{$item}++; } # collecting unique alleles
				delete $temphash{"possibleOthers"}; # ingore the "possibleOthers" allele
				@call = keys %temphash;

				if ($1 eq '-/-') { # a simple deletion ["allele-call=-/-"] (for the end position, "$end" is already not used)
					print $chr, "\t", $pos, "\t", $pos+length($refcall)-1, "\t", $refcall, "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
				} elsif ($refcall eq '-') { # a simple insertion (single or multiple allele) ["reference=-"]
					for my $i (0 .. @call-1) {					    
					    	next if ($refcall eq $call[$i]);
						print $chr, "\t", $pos, "\t", $pos, "\t", '-', "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
						$i > 0 and $countallvar++;
					}
				} else { # an indel that may have several alleles, or may require a block substitution representation
					for my $i (0 .. @call-1) {
					    	next if ($refcall eq $call[$i]);
					    	# for the end position, "$pos+length($call[0])-1" is already not used.
						print $chr, "\t", $pos, "\t", $pos+length($refcall)-1, "\t", $refcall, "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
						$i > 0 and $countallvar++;
					}
				} 
			} elsif ($attribute =~ m#allele\-call=([\w\/]+)#) {
				@call = split (/\//, $1);
				
				if (@call == 1) {		#a simple deletion
					print $chr, "\t", $pos, "\t", $end, "\t", $call[0], "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
				} elsif ($call[0] eq '') {	#a simple insertion (single or multiple allele)
					for my $i (1 .. @call-1) {
						print $chr, "\t", $pos, "\t", $pos, "\t", '-', "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
						$i > 1 and $countallvar++;
					}
				} else {			#an indel that may have several alleles, or may require a block substitution representation
					for my $i (1 .. @call-1) {
						print $chr, "\t", $pos, "\t", $pos+length($call[0])-1, "\t", $call[0], "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
						$i > 1 and $countallvar++;
					}
				}
			} else {
				$call = '0';
				print $chr, "\t", $pos, "\t", $end, "\t", $call, "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
			}
		} else {
			die "Error: unrecognizable genotype calling program encountered (valid types are SOLiD_diBayes, AB_CNV_PIPELINE, AB_SOLiD Large Indel Tool, AB_SOLiD Small Indel Tool): <$_>\n";
		}
			
		$countvar++;		#variation positions
		$countallvar++;		#all variants (maybe several at one variation position)
	}
	print STDERR "NOTICE: Finished processing $variantfile with $countline input lines\n";
	print STDERR "NOTICE: Wrote variants in $countvar variation positions ($countallvar variants considering multi-allelic ones)\n";
	$unknown_count and print STDERR "WARNING: $unknown_count variants with 'unknown' variation type were skipped\n";
}



sub convertSOAP {
	my ($variantfile) = @_;
	my ($countline, $countvar, @other);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		
		my @field = split (/\t/, $_);
		if (@field == 18) {		#snp file
			my ($chr, $pos, $wt, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
	
			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
		} elsif (@field == 6) {		#indel file
			my ($chr, $pos, $strand, $indellen, $call, $homo) = @field;
			$homo eq 'Homo' and $homo = 'hom';
			$homo eq 'Hete' and $homo = 'het';
			$chr =~ s/^chr//;
	
			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
	
			if ($indellen =~ m/^\+(\d+)$/) {		#insertion
				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
				print join("\t", $chr, $pos, $pos, '-', $call, $homo), "\n";
			} elsif ($indellen =~ m/^\-(\d+)$/) {		#deletion
				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
				print join("\t", $chr, $pos, $pos+$1-1, $call, '-', $homo), "\n";
			} else {
				die "Error: invalid record found in SOAPindel file: <$_>\n";
			}
		} else {
			die "Error: invalid record found in $variantfile (18 or 6 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
		$countvar++;
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}


sub convertANNOVAR {
	my ($variantfile) = @_;
	my ($countline, $countvar, $countsnp, $invalid) = (0, 0, 0, 0);
	my ($countti, $counttv) = (0, 0);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	while (<VAR>) {
		$countline++;
		s/[\r\n]+$//; 
		my @field = split (/\s+/, $_);
		@field >= 5 or die "Error: invalid record found in annovar input file (at least 5 tab or space delimited fields expected): <$_>\n";
		my ($chr, $start, $end, $ref, $obs) = @field;
		if ($ref =~ m/^[ACGT]$/ and $obs =~ m/^[ACGT]$/) {
			if ($ref eq 'A' and $obs eq 'G' or $ref eq 'G' and $obs eq 'A' or $ref eq 'C' and $obs eq 'T' or $ref eq 'T' and $obs eq 'C') {
				$countti++;
			} else {
				$counttv++;
			}
			$countsnp++;
		}
		
		($ref, $obs) = (uc $ref, uc $obs);
		$chr =~ s/^chr//;
		if ($chr =~ m/[^\w\.]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {		#chr name could contain . (example: GL000212.1)
			$invalid++;
		} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
			or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
			or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
			or $start =~ m/[^\d]/ 			#start is not a number
			or $end =~ m/[^\d]/ 			#end is not a number
			or $start > $end			#start is more than end
			or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
			or $ref eq '-' and $start != $end	#length mismatch for insertion
			) {
			$invalid++;
		}
		print "$_\n";
		$countvar++;
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
	$invalid and print STDERR "WARNING: $invalid input lines have invalid formats\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions\n";
}

sub convertMAQSNP {
	my ($variantfile) = @_;
	my ($countline, $countvar, @other);
	
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		
		my @field = split (/\t/, $_);
		my @other = ();
		if (@field == 12) {					#SNPs
			my ($chr, $pos, $wt, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
	
			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
			$countvar++;
		} elsif (@field == 13) {				#indels; the deletion start site do not need changes; the duplication start site need additional investigation by ANNOVAR developers
			my ($chr, $pos, $type, $numread, $call, @other) = @field;
			$chr =~ s/^chr//;
	
			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
	
			my @obs = split (/:/, $call);
			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] =~ m/^\-(\d+)/) {		#deletion
				my $len = $1;
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[1], "\t", '-', "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[0] =~ m/^(\d+)/) {	#insertion
				my $len = $1;
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";	#2011jul12: changed pos to pos-1 for insertions
			}
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (12 or 13 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}

sub convertCASAVA {
	my ($variantfile, $chr) = @_;
	my ($countline, $countvar, @other);
	
	my ($intype);
	my ($pos_index, $call_index, $reference_index, $type_index, $score_index, $total_index, $used_index);
	my ($ref_indel_index, $quality_index, $maxgtype_index, $bp1_reads_index, $ref_reads_index, $indel_reads_index, $other_reads_index);
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my @field;

		if (m/^#/) {
			s/^#//;
			if (s/^\$\sCOLUMNS\s//) {
				@field = split (/\s+/, $_);
			} else {
				@field = split (/\t/, $_);
			}
			if (m/\bposition\b/ or m/\bpos\b/) {
				for my $i (0 .. @field-1) {
					if ($field[$i] eq 'position' or $field[$i] eq 'pos') {
						$pos_index = $i;
					} elsif ($field[$i] eq 'max_gt|poly_site') {		#this has priority over max_gt per se
						$intype = 'snp';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$call_index = $i;
					} elsif ($field[$i] eq 'modified_call' or $field[$i] eq 'max_gt') {
						$intype = 'snp';
						$call_index = $i;
					} elsif ($field[$i] eq 'reference' or $field[$i] eq 'ref') {
						$reference_index = $i;
					} elsif ($field[$i] eq 'type') {
						$type_index = $i;
					} elsif ($field[$i] eq 'score') {
						$score_index = $i;
					} elsif ($field[$i] eq 'total') {
						$total_index = $i;
					} elsif ($field[$i] eq 'used') {
						$used_index = $i;
					} elsif ($field[$i] eq 'ref/indel') {
						$intype = 'indel';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$ref_indel_index = $i;
					} elsif ($field[$i] eq 'Q(max_gt|poly_site)') {	#this has priority over Q(max)
						$quality_index = $i;
					} elsif ($field[$i] eq 'Q(max_gt)') {	#this has priority over Q(max)
						$quality_index = $i;
					} elsif ($field[$i] eq 'Q(indel)') {
						$quality_index = $i;
					} elsif ($field[$i] eq 'max_gtype') {
						$maxgtype_index = $i;
					} elsif ($field[$i] eq 'bp1_reads') {
						$bp1_reads_index = $i;
					} elsif ($field[$i] eq 'ref_reads') {
						$ref_reads_index = $i;
					} elsif ($field[$i] eq 'indel_reads') {
						$indel_reads_index = $i;
					} elsif ($field[$i] eq 'other_reads') {
						$other_reads_index = $i;
					}
				}
			}
			next;
		}
		
		##$ COLUMNS seq_name pos bcalls_used bcalls_filt ref Q(snp) max_gt Q(max_gt) max_gt|poly_site Q(max_gt|poly_site) A_used C_used G_used T_used
		#chr21.fa	9411785	1	0	G	11	GT	3	GT	3	0	0	0	1
		#chr21.fa	9414658	1	0	T	10	CT	3	CT	3	0	1	0	0
		#chr21.fa	9415181	2	0	C	52	TT	5	TT	5	0	0	0	2
		#chr21.fa	9415317	2	0	C	6	CT	6	CT	34	0	1	0	1

		
		$intype or die "Error: unable to recognize the correct type of the input file (make sure that header line is present in $variantfile)\n";
		@field = split (/\t/, $_);
		
		if ($intype eq 'snp') {					#SNPs
			defined $pos_index and defined $reference_index and defined $call_index or die "Error: unalbe to find the position, reference and modified_call column header in $variantfile\n";
			my ($pos, $wt, $obs) = @field[$pos_index, $reference_index, $call_index];
			my (@other);
			defined $pos and defined $wt and defined $obs or die;
			$includeinfo and @other = @field;
			
			length ($obs) == 1 and $obs .= $obs;
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed allele $obs should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
			$countvar++;
		} elsif ($intype eq 'indel') {				#indels
			defined $pos_index and defined $ref_indel_index and defined $maxgtype_index or die "Error: unable to find the pos, ref_indel and max_gtype column header in $variantfile\n";
			my ($pos, $call, $hom, @other) = @field[$pos_index, $ref_indel_index, $maxgtype_index];
			$includeinfo and @other = @field;

			#hg19 coordinate below; insertion needs position adjustment!!! deletion is fine
			#948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
			#978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
			#1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
			#1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0
			#185137455     11I10M2I        TATGTGTCCT      -----------TTTTTTATTT--/AAATGATAGACTTTTTTTTTTAA ATTTCAGAAA      1126    het     988     hom    45       20      24      7       N/A     0       0
			#1276931 2D41M4I CACACACATG      CACACACACGCACACACGTGCAATGTGAAAACACCTCATGCAG----/--CACACACGCACACACGTGCAATGTGAAAACACCTCATGCAGACAC ACACATGCAC      548     hom     16      het     8       0       11      11      N/A     0       0
			
			my @obs = split (/\//, $call);
			@obs == 2 or die "Error: observed indel allele $call should correspond to two alleles: <$_>\n";
			if ($obs[0] =~ m/^\-+$/) {		#insertion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif ($obs[1] =~ m/^\-+$/) {		#deletion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[0], "\t", '-', "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif (length ($obs[0]) eq length ($obs[1])) {	#block substitution
				$obs[0] =~ s/\-//g;
				$obs[1] =~ s/\-//g;
				print $chr, "\t", $pos, "\t", $pos+length($obs[0])-1, "\t", $obs[0], "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} else {
				die "Error: invalid record found in indel line: <$_>\n";
			}
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (11 or 15 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}

sub convertVCF4 {
	my ($variantfile) = @_;
	
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;
	
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		$countline++;
		
		if (m/^##fileformat=VCFv(\d+\.)/) {
			$1<4 and print STDERR "ERROR: Input file is not in VCF version 4 format but is $_" and exit;
		}
		if (m/^##UnifiedGenotyper/) {
			$source_program = 'gatksnp';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK UnifiedGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth\n";
			$fraction and print STDERR "WARNING: the --fraction argument will be ignored for GATK SNP calls!!!\n";
			$confraction and print STDERR "WARNING: the --confraction argument will be ignored for GATK SNP calls!!!\n";
		}
		if (m/^##IndelGenotyper/) {
			$source_program = 'gatkindel';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK IndelGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality\n";
		}
		
		if (not m/^#/ and not $source_program) {	#finished reading header line but did not detect the source program
			$includeinfo or print STDERR "NOTICE: for SNPs, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth, if these information can be recognized automatically\n";
			$includeinfo or print STDERR "NOTICE: for indels, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality, if these information can be recognized automatically\n";
			$source_program = 'unknown';
		}
		
		m/^#/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		my $otherinfo = $_;	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
	
		#format description: http://www.1000genomes.org/wiki/Analysis/vcf4.0
		#standard VCF4 should have 8 columns, but some software may produce more columns (for example, for genotype calls). The first 8 columns should follow the specification
		
		#example of VCF4 generated by GATK SNP caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       55      .       T       G       34.82   .       DP=2;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=14.16;MQ0=0;QD=17.41;SB=-10.00     GT:DP:GL:GQ     0/1:1:-6.66,-0.30,-0.00:1.76
		#1       2646    .       G       A       40.91   .       DP=4;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=7.50;MQ0=3;QD=10.23;SB=-10.00      GT:DP:GL:GQ     0/1:1:-7.27,-0.30,-0.00:1.76
		
		#example of VCF4 generated by GATK indel caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       2525324 .       G       GC      .       PASS    AC=5,5;DP=12;MM=4.8,3.7142856;MQ=29.0,42.285713;NQSBQ=33.0,46.463768;NQSMM=0.24,0.20289855;SC=0,5,1,6  GT       0/1
		#1       3553372 .       GC      G       .       PASS    AC=6,6;DP=6;MM=0.8333333,0.0;MQ=60.0,0.0;NQSBQ=63.533333,0.0;NQSMM=0.0,0.0;SC=0,6,0,0   GT      1/0
		#1       6093011 .       CG      C       .       PASS    AC=31,31;DP=32;MM=0.7096774,2.0;MQ=59.64516,60.0;NQSBQ=64.192184,39.666668;NQSMM=0.0,0.11111111;SC=23,8,0,1     GT      1/0
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       533     .       G       C       .       PASS    AA=.;AC=6;AN=120;DP=423
		#1       41342   .       T       A       .       PASS    AA=.;AC=29;AN=120;DP=188
		#1       41791   .       G       A       .       PASS    AA=.;AC=5;AN=120;DP=192
		#1       44449   .       T       C       .       PASS    AA=C;AC=2;AN=120;DP=166
		#1       44539   rs2462492       C       T       .       PASS    AA=T;AC=2;AN=120;DP=131    
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       1000153 .       TCACA   T       100     PASS    AF=0.115095;HP=1;NF=16;NR=13;NS=52;CA=0;DP=615
		#1       1000906 .       CA      C       48      PASS    AF=0.0772696;HP=1;NF=2;NR=9;NS=51;CA=0;DP=281
		#1       1000950 rs60561655;-/G  CG      C       100     PASS    AF=0.447771;HP=5;DB;NF=10;NR=20;NS=50;CA=M;DP=291
		#1       1010786 rs36095298;-/G,mills,venter     A       AG      100     PASS    AF=0.774334;HP=1;DB;NF=21;NR=27;NS=51;CA=0;DP=306
		#1       1026158 .       T       TGGGGG  100     PASS    AF=0.115637;HP=1;NF=5;NR=2;NS=52;CA=0;DP=591
                
                #example of VCF4 generated by SamTools mpileup (Note that GT was not the first field in the FORMAT string)
                ##fileformat=VCFv4.0
		##samtoolsVersion=0.1.16 (r963:234)
		##fileformat=VCFv4.0
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  1247MFL0003.NOVO.srt.bam
		#chr1    14574   .       A       G       3.54    .       "DP=3;AF1=0.4999;CI95=0.5,0.5;DP4=1,0,2,0;MQ=21;PV4=1,1,1,1"    PL:GT:GQ        "31,0,34:0/1:32"
		#chr1    14930   .       A       G       37      .       "DP=19;AF1=0.5;CI95=0.5,0.5;DP4=7,5,5,1;MQ=25;PV4=0.6,6.3e-05,1,0.23"   PL:GT:GQ        "67,0,103:0/1:70"
		#chr1    16495   .       G       C       28      .       "DP=4;AF1=0.5;CI95=0.5,0.5;DP4=0,0,4,0;MQ=32"   PL:GT:GQ        "70,12,70:0/1:58"
		#chr1    59040   .       T       C       4.77    .       "DP=4;AF1=0.4999;CI95=0.5,0.5;DP4=0,2,2,0;MQ=22;PV4=0.33,0.21,1,1"      PL:GT:GQ        "33,0,39:0/1:35"
		#chr1    69270   .       A       G       46      .       "DP=20;AF1=0.5;CI95=0.5,0.5;DP4=2,0,18,0;MQ=24;PV4=1,1,1,0.28"  PL:GT:GQ        "94,18,100:0/1:78"
		#chr1    69511   .       A       G       24      .       "DP=5;AF1=0.5;CI95=0.5,0.5;DP4=1,0,2,1;MQ=25;PV4=1,0.46,1,0.039"        PL:GT:GQ        "54,0,57:0/1:55"
                
		#reserved VCF4 sub-fields in the INFO field
		#    * AA ancestral allele
		#    * AC allele count in genotypes, for each ALT allele, in the same order as listed
		#    * AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
		#    * AN total number of alleles in called genotypes
		#    * BQ RMS base quality at this position
		#    * CIGAR cigar string describing how to align an alternate allele to the reference allele
		#    * DB dbSNP membership
		#    * DP combined depth across samples, e.g. DP=154
		#    * END end position of the variant described in this record (esp. for CNVs)
		#    * H2 membership in hapmap2
		#    * MQ RMS mapping quality, e.g. MQ=52
		#    * MQ0 Number of MAPQ == 0 reads covering this record
		#    * NS Number of samples with data
		#    * SB strand bias at this position
		#    * SOMATIC indicates that the record is a somatic mutation, for cancer genomics
		#    * VALIDATED validated by follow-up experiment


		#SAMtools/BCFtools specific information
		#SAMtools/BCFtools may write the following tags in the INFO field in VCF/BCF.
		#Tag	Description
		#I16	16 integers:
		#1	#reference Q13 bases on the forward strand 	2	#reference Q13 bases on the reverse strand
		#3	#non-ref Q13 bases on the forward strand 	4	#non-ref Q13 bases on the reverse strand
		#5	sum of reference base qualities 	6	sum of squares of reference base qualities
		#7	sum of non-ref base qualities 	8	sum of squares of non-ref base qualities
		#9	sum of ref mapping qualities 	10	sum of squares of ref mapping qualities
		#11	sum of non-ref mapping qualities 	12	sum of squares of non-ref mapping qualities
		#13	sum of tail distance for ref bases 	14	sum of squares of tail distance for ref bases
		#15	sum of tail distance for non-ref bases 	16	sum of squares of tail distance for non-ref
		#INDEL	Indicating the variant is an INDEL.
		#DP	The number of reads covering or bridging POS.
		#DP4	Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling. Sum can be smaller than DP because low-quality bases are not counted.
		#PV4	P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)
		#FQ	Consensus quality. If positive, FQ equals the phred-scaled probability of there being two or more different alleles. If negative, FQ equals the minus phred-scaled probability of all chromosomes being identical. Notably, given one sample, FQ is positive at hets and negative at homs.
		#AF1	EM estimate of the site allele frequency of the strongest non-reference allele.
		#CI95	Equal-tail (Bayesian) credible interval of the site allele frequency at the 95% level.
		#PC2	Phred-scaled probability of the alternate allele frequency of group1 samples being larger (,smaller) than of group2 samples.
		#PCHI2	Posterior weighted chi^2 P-value between group1 and group2 samples. This P-value is conservative.
		#QCHI2	Phred-scaled PCHI2
		#RP	Number of permutations yeilding a smaller PCHI2

		#example of triallelic variants generated by mpileup/bcftools
		#1       156706559       .       A       C,G     114     .	DP=20;AF1=1;CI95=1,1;DP4=0,0,1,19;MQ=60;FQ=-63  GT:PL:GQ	1/2:237,126,90,162,0,138:99
		#6       31129642        .       A       G,C     76      .	DP=31;AF1=1;CI95=1,1;DP4=0,0,28,3;MQ=60;FQ=-75  GT:PL:GQ	1/2:255,194,146,164,0,119:99
		#1       11297762        .       T       C,A     98      .	DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ	1/1:131,51,0,120,28,117:99
		
		my @field=split(/\t/,$_);
		@field >=8 or die "Error: invalid record found in VCF4 file (at least 8 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $ID, $ref_allele, $mut_allele, $quality_score, $filter, $info, $format, $sample) = @field;
		my ($end);
		my ($mut_allele2, $zygosity);
		
		if ($filterword) {		#ensure that the filter field contains the filterword
			$filter =~ m/\b$filterword\b/i or next;
		}
		
		$info =~ s/^"//; $info =~ s/"$//;
		
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		#if ($ID eq '.' || $ID =~ /^rs/) {		#per MISHIMA, Hiroyuki suggestion (vcf4's third column (ID column) are not always ".")
		#	$end = $start;				#this block is commented out on 2011feb19
		#}
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			next;
		}
		
		if ($mut_allele =~ m/([^,]+),([\w,]+)/) {	#there could be more than two alternative alleles
			$mut_allele = $1;
			$mut_allele2 = $2;
		}
		
		if(length($ref_allele)==1 && length($mut_allele)==1) {  	### output snv
			if ($ref_allele =~ m/[^ACGTacgt]/ ) {
				print STDERR "WARNING: invalid allele record found in VCF4 file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
				$ref_allele = 0;
			}
			if ( $mut_allele =~ m/[^ACGTacgt]/) {
				print STDERR "WARNING: invalid allele record found in VCF4 file (ACGT expected): <$ref_allele> and <$mut_allele> in line <$_>\n";
				$mut_allele = 0;
			}
				
			my ($unfiltered_read_depth) = $info =~ /DP=(\d+)/;
			my ($MappingQuality) = $info =~ /MQ=([^;]+)/; 
			my ($QualityByDepth) = $info =~ /QD=([^;]+)/;		
			

			
			if ($coverage) {
				defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
				if ($maxcoverage) {
					defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
				}
			}
			
			if ($snpqual) {
				defined $QualityByDepth and $QualityByDepth >= $snpqual || next;		#the QD was used here as quality score
			}			
			
			if (defined $format) {
				my @format = split (/:/, $format);
				undef $gtpos;
				for my $i (0 .. @format-1) {
					if ($format[$i] eq 'GT') {
						$gtpos = $i;
						last;
					}
				}
				if (defined $sample and defined $gtpos) {
					my @sample = split (/:/, $sample);
					#if ($sample[$gtpos] =~ m#^0/1# or $sample[$gtpos] =~ m#^1/0#) {
					#	$zygosity = 'het';
					#	$counthet++;
					#} elsif ($sample[$gtpos] =~ m#^1/1#) {
					#	$zygosity = 'hom';
					#	$counthom++;
					#change the above lines to the following on 20120618
					if ($sample[$gtpos] =~ m#^(\d)/(\d)#) {
						if ($1 == $2) {
							$zygosity = 'hom';
							$counthom++;
						} else {
							$zygosity = 'het';
							$counthet++;
						}
					
					} else {
						$zygosity = 'unknown';
						$countunknown++;
					}
				} else {		#sometimes the input VCF file does not contain the GT field!!!
					$zygosity = 'unknown';
					$countunknown++;
				}
			} else {
				$zygosity = 'unknown';
				$countunknown++;
			}

			#the subject is called as homozygous for the first alternative allele (genotype 1/1. i.e. C/C), but since there was one read containing A, samtools still keep both alleles in the VCF file (but gives a very low probabilities for it).
			#1       11297762        .       T       C,A     98      . DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ 1/1:131,51,0,120,28,117:99			
			if ($mut_allele2 and $zygosity eq 'hom') {
				$mut_allele2 = '';
			}

			if (not $mut_allele2) {
				if ($ref_allele eq 'A' and $mut_allele eq 'G' or $ref_allele eq 'G' and $mut_allele eq 'A' or $ref_allele eq 'C' and $mut_allele eq 'T' or $ref_allele eq 'T' and $mut_allele eq 'C') {
					$countti++;
					
				} else {
					$counttv++;
				}
			}
			
			#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
			if ($includeinfo) {
				print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
			} else {
				print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', "\n";
			}
			
			if ($allallele) {
				if ($mut_allele2) {
					my @mut_allele2 = split (/,/, $mut_allele2);
					for my $i (0 .. @mut_allele2-1) {
						if ($includeinfo) {
							print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
						} else {
							print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele2[$i], "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', "\n";
						}
					}
				}
			}
			
			$countsnp++;
		} elsif (length($ref_allele) > 1 || length($mut_allele) > 1) {  ### output indel
			my ($indel_read_depth1, $indel_read_depth2) = $info =~ /\bAC=([^,;]+),([^,;]+)/;		#number of reads supporting consensus indel, any indel
			my ($unfiltered_read_depth) = $info =~ /\bDP=(\d+)/;
			my ($MappingQuality) = $info =~ /\bMQ=([^;]+)/;
			my ($QualityByDepth) = $info =~ /\bQD=([^;]+)/;	
					
			if ($coverage) {
				defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
				if ($maxcoverage) {
					defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
				}
			}
			
			if ($snpqual) {
				defined $QualityByDepth and $QualityByDepth >= $snpqual || next;		#the QD was used here as quality score
			}
			
			if (defined $indel_read_depth1 and $unfiltered_read_depth) {		#deleted "defined" before $unfiltered_read_depth on 2012may25
				$indel_read_depth1/$unfiltered_read_depth >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
				if ($indel_read_depth2) {
					$indel_read_depth1/$indel_read_depth2 >= $confraction or next;
				}
			}
			

			#example VCF4 records below:
			#20      2       .       TCG     T       .       PASS    DP=100
			#Chr1    5473    .       AT      ATT     23.5    .       INDEL;DP=16;AF1=0.5;CI95=0.5,0.5;DP4=4,2,3,1;MQ=42;PV4=1,0.41,0.042,0.24
			#Chr1    6498    .       ATTTT   ATTTTT  53.5    .       INDEL;DP=9;AF1=1;CI95=1,1;DP4=0,0,5,3;MQ=28
			
			if(length($ref_allele) > length ($mut_allele)) { 		# deletion or block substitution
				my $head = substr($ref_allele, 0, length ($mut_allele));
				if ($head eq $mut_allele) {
					print $chr,"\t";
					print $start+length($head),"\t";
					print $start+length($ref_allele)-1,"\t";
					
					my $ref_allele1 = substr ($ref_allele, length ($mut_allele));
					print $ref_allele1,"\t";
					print "-";
				} else {
					print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
				}
			} elsif(length($mut_allele) >= length ($ref_allele)) { 		# insertion or block substitution
				my $head = substr ($mut_allele, 0, length ($ref_allele));
				if ($head eq $ref_allele) {
					print $chr,"\t";	
					print $start+length($ref_allele)-1,"\t";
					print $start+length($ref_allele)-1,"\t";
					
					$mut_allele = substr ($mut_allele, length ($ref_allele));
					print "-\t";
					print $mut_allele;
				} else {
					print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele;
				}
			}
			

			if (defined $format) {
				my @format = split (/:/, $format);
				undef $gtpos;
				for my $i (0 .. @format-1) {
					if ($format[$i] eq 'GT') {
						$gtpos = $i;
						last;
					}
				}
				if (defined $sample and defined $gtpos) {
					my @sample = split (/:/, $sample);
					if ($sample[$gtpos] =~ m#^0/1# or $sample[$gtpos] =~ m#^1/0#) {
						$zygosity = 'het';
						$counthet++;
					} elsif ($sample[$gtpos] =~ m#^1/1#) {
						$zygosity = 'hom';
						$counthom++;
					} else {
						$zygosity = 'unknown';
						$countunknown++;
					}
				}
			} else {
				$zygosity = 'unknown';
				$countunknown++;
			}

			#commented out on 2011May21
			#if (defined $sample) {
			#	if ($sample =~ m#^0/1# or $sample =~ m#^1/0#) {
			#		print "\thet";
			#		$counthet++;
			#	} elsif ($sample =~ m#^1/1#) {
			#		print "\thom";
			#		$counthom++;
			#	} # BEGIN ARQ
			#	elsif ($sample =~ m#^./.#) {
			#		print "\tunknown";
			#		$countunknown++;
			#	} # END ARQ
			#}
			
			
			if ($includeinfo) {
				 print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
			} else {
				print "\t$zygosity";
				defined $quality_score and print "\t", $quality_score;
				defined $unfiltered_read_depth and print "\t", $unfiltered_read_depth;
				
				#defined $indel_read_depth1 and print "\t", $indel_read_depth1;		#commented out Nov 2011
				defined $MappingQuality and print "\t", $MappingQuality;
				defined $QualityByDepth and print "\t", $QualityByDepth;		#added in Nov 2011 to be consistent with SNP output
				#$includeinfo and print "\t", $otherinfo;	#commented Sep 2011
			}
			print "\n";
			$countindel++;


			#do the same thing again, exactly like above, except that we work on second mutation;
			#in the future, consider rewrite this paragraph to make the code more elegant	
			if ($allallele and $mut_allele2) {
				my @mut_allele2 = split (/,/, $mut_allele2);
				for my $mut_allele2 (@mut_allele2) {
					if(length($ref_allele) > length ($mut_allele2)) { 		# deletion or block substitution
						my $head = substr($ref_allele, 0, length ($mut_allele2));
						if ($head eq $mut_allele2) {
							print $chr,"\t";
							print $start+length($head),"\t";
							print $start+length($ref_allele)-1,"\t";
							
							my $ref_allele1 = substr ($ref_allele, length ($mut_allele2));
							print $ref_allele1,"\t";
							print "-";
						} else {
							print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
						}
					} elsif(length($mut_allele2) > length ($ref_allele)) { 		# insertion or block substitution
						my $head = substr ($mut_allele2, 0, length ($ref_allele));
						if ($head eq $ref_allele) {
							print $chr,"\t";	
							print $start+length($ref_allele)-1,"\t";
							print $start+length($ref_allele)-1,"\t";
							
							$mut_allele2 = substr ($mut_allele2, length ($ref_allele));
							print "-\t";
							print $mut_allele2;
						} else {
							print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
						}
					} else {		#identical length of alleles
						print $chr, "\t", $start, "\t", $start+length($ref_allele)-1, "\t", $ref_allele, "\t", $mut_allele2;
					}

					if (defined $sample) {
						if ($sample =~ m#^0/1# or $sample =~ m#^1/0#) {
							$zygosity = "het";
							$counthet++;
						} elsif ($sample =~ m#^1/1#) {
							$zygosity =  "hom";
							$counthom++;
						} # BEGIN ARQ
						elsif ($sample =~ m#^./.#) {
							$zygosity = "unknown";
							$countunknown++;
						} # END ARQ
					}
										
					if ($includeinfo) {
						print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
					} else {
						print "\t", $zygosity;
						print "\t", $quality_score;
						defined $unfiltered_read_depth and print "\t", $unfiltered_read_depth;
						
						defined $indel_read_depth1 and print "\t", $indel_read_depth1;
						defined $MappingQuality and print "\t", $MappingQuality;
						#$includeinfo and print "\t", $otherinfo;
					}
					print "\n";

				}
			}
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
}


=head1 SYNOPSIS

 convert2annovar.pl [arguments] <variantfile>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --format <string>		input format (default: pileup)
            --outfile <file>		output file name (default: STDOUT)
            --snpqual <float>		quality score threshold in pileup file (default: 20)
            --snppvalue <float>		SNP P-value threshold in GFF3-SOLiD file (default: 1)
            --coverage <int>		read coverage threshold in pileup file (default: 0)
            --maxcoverage <int>		maximum coverage threshold (default: none)
            --includeinfo		include supporting information in output
            --chr <string>		specify the chromosome (for CASAVA format)
            --chrmt <string>		chr identifier for mitochondria (default: M)
            --altcov <int>		alternative allele coverage threshold (for pileup format)
            --allelicfrac		print out allelic fraction rather than het/hom status (for pileup format)
            --fraction <float>		minimum allelic fraction to claim a mutation (for pileup/vcf4_indel format)
            --species <string>		if human, convert chr23/24/25 to X/Y/M (for gff3-solid format)
            --filter <string>		output variants with this filter (case insensitive, for vcf4 format)
            --confraction <float>	minimum consensus indel / all indel fraction (for vcf4 format)
            --allallele			print all alleles when multiple calls are present (for vcf4 format)
            --withzyg			print zygosity when -includeinfo is used (for vcf4 format)

 Function: convert variant call file generated from various software programs 
 into ANNOVAR input format
 
 Example: convert2annovar.pl -format pileup -outfile variant.query variant.pileup
          convert2annovar.pl -format cg -outfile variant.query variant.cg
          convert2annovar.pl -format gff3-solid -outfile variant.query variant.snp.gff
          convert2annovar.pl -format soap variant.snp > variant.avinput
          convert2annovar.pl -format maq variant.snp > variant.avinput
          convert2annovar.pl -format casava -chr 1 variant.snp > variant.avinput
          convert2annovar.pl -format vcf4 variantfile > variant.avinput
          convert2annovar.pl -format vcf4 -filter pass variantfile > variant.avinput

 Version: $LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--format>

the format of the input files.

=item B<--outfile>

specify the output file name. By default, output is written to STDOUT.

=item B<--snpqual>

quality score threshold in the pileup file, such that variant calls with lower 
quality scores will not be printed out in the output file. When VCF4 file is 
used, this argument works on the Quality-by-Depth measure, rather than the raw 
quality measure.

=item B<--coverage>

read coverage threshold in the pileup file, such that variants calls generated 
with lower coverage will not be printed in the output file.

=item B<--includeinfo>

specify that the output should contain additional information in the input line. 
By default, only the chr, start, end, reference allele, observed allele and 
homozygosity status are included in output files.

=item B<--chr>

specify the chromosome for CASAVA format

=item B<--chrmt>

specify the name of mitochondria chromosome (default is MT)

=item B<--altcov>

the minimum coverage of the alternative (mutated) allele to be printed out in 
output

=item B<--fraction>

specify the minimum fraction of alternative allele, to print out the mutation. 
For example, a site has 10 reads, 3 supports alternative allele. A -fraction of 
0.4 will not allow the mutation to be printed out.

=item B<--species>

specify the species from which the sequencing data is obtained. For the GFF3-
SOLiD format, when species is human, the chromosome 23, 24 and 25 will be 
converted to X, Y and M, respectively.

=item B<--filter>

for VCF4 file, only print out variant calls with this filter annotated. For 
example, if using GATK VariantFiltration walker, you will see PASS, 
GATKStandard, HARD_TO_VALIDATE, etc in the filter field. Using 'pass' as a 
filter is recommended in this case.

=item B<--confraction>

consesus indel fraction, calculated as reads supporting consensus indels divided 
by reads supporting any indels

=item B<--allallele>

print all alleles for mutations at a locus, rather than the first allele, if the 
input VCF4 file contains multiple alternative alleles for a mutation. By 
default, this option is off. When it is on, two lines will be printed out in the 
output, and both will have the same quality scores as VCF4 does not provide 
separate quality scores for individual alleles.

=back

=head1 DESCRIPTION

This program is used to convert variant call file generated from various 
software programs into ANNOVAR input format. Currently, the program can handle 
Samtools genotype-calling pileup format, Solid GFF format, Complete Genomics 
variant format, SOAP format. These formats are described below.

=over 8

=item * B<pileup format>

The pileup format can be produced by the Samtools genotyping calling subroutine. 
Note that the phrase pileup format can be used in several instances, and here I 
am only referring to the pileup files that contains the actual genotype calls. 

Using SamTools, given an alignment file in BAM format, a pileup file with 
genotype calls can be produced by the command below:

	samtools pileup -vcf ref.fa aln.bam> raw.pileup
	samtools.pl varFilter raw.pileup > final.pileup

ANNOVAR will automatically filter the pileup file so that only SNPs reaching a 
quality threshold are printed out (default is 20, use --snpqual argument to 
change this). Most likely, users may want to also apply a coverage threshold, 
such that SNPs calls from only a few reads are not considered. This can be 
achieved using the -coverage argument (default value is 0).

An example of pileup files for SNPs is shown below:

	chr1 556674 G G 54 0 60 16 a,.....,...,.... (B%A+%7B;0;%=B<:
	chr1 556675 C C 55 0 60 16 ,,..A..,...,.... CB%%5%,A/+,%....
	chr1 556676 C C 59 0 60 16 g,.....,...,.... .B%%.%.?.=/%...1
	chr1 556677 G G 75 0 60 16 ,$,.....,...,.... .B%%9%5A6?)%;?:<
	chr1 556678 G K 60 60 60 24 ,$.....,...,....^~t^~t^~t^~t^~t^~t^~t^~t^~t B%%B%<A;AA%??<=??;BA%B89
	chr1 556679 C C 61 0 60 23 .....a...a....,,,,,,,,, %%1%&?*:2%*&)(89/1A@B@@
	chr1 556680 G K 88 93 60 23 ..A..,..A,....ttttttttt %%)%7B:B0%55:7=>>A@B?B;
	chr1 556681 C C 102 0 60 25 .$....,...,....,,,,,,,,,^~,^~. %%3%.B*4.%.34.6./B=?@@>5.
	chr1 556682 A A 70 0 60 24 ...C,...,....,,,,,,,,,,. %:%(B:A4%7A?;A><<999=<<
	chr1 556683 G G 99 0 60 24 ....,...,....,,,,,,,,,,. %A%3B@%?%C?AB@BB/./-1A7?

The columns are chromosome, 1-based coordinate, reference base, consensus base, 
consensus quality, SNP quality, maximum mapping quality of the reads covering 
the sites, the number of reads covering the site, read bases and base qualities.

An example of pileup files for indels is shown below:

	seq2  156 *  +AG/+AG  71  252  99  11  +AG  *  3  8  0

ANNOVAR automatically recognizes both SNPs and indels in pileup file, and process them correctly.

=item * B<GFF3-SOLiD format>

The SOLiD provides a GFF3-compatible format for SNPs, indels and structural 
variants. A typical example file is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version 2
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files Yoruban_snp_10x.txt
	##run-path 
	chr_name        AB_SOLiD SNP caller     SNP     coord   coord   1       .       .       coverage=# cov;ref_base=ref;ref_score=score;ref_confi=confi;ref_single=Single;ref_paired=Paired;consen_base=consen;consen_score=score;consen_confi=conf;consen_single=Single;consen_paired=Paired;rs_id=rs_id,dbSNP129
	1       AB_SOLiD SNP caller     SNP     997     997     1       .       .       coverage=3;ref_base=A;ref_score=0.3284;ref_confi=0.9142;ref_single=0/0;ref_paired=1/1;consen_base=G;consen_score=0.6716;consen_confi=0.9349;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     2061    2061    1       .       .       coverage=2;ref_base=G;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=C;consen_score=1.0000;consen_confi=0.8985;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     4770    4770    1       .       .       coverage=2;ref_base=A;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=G;consen_score=1.0000;consen_confi=0.8854;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     4793    4793    1       .       .       coverage=14;ref_base=A;ref_score=0.0723;ref_confi=0.8746;ref_single=0/0;ref_paired=1/1;consen_base=G;consen_score=0.6549;consen_confi=0.8798;consen_single=0/0;consen_paired=9/9
	1       AB_SOLiD SNP caller     SNP     6241    6241    1       .       .       coverage=2;ref_base=T;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=C;consen_score=1.0000;consen_confi=0.7839;consen_single=0/0;consen_paired=2/2
	
Newer version of ABI BioScope now use diBayes caller, and the output file is given below:

	##gff-version 3
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##List of SNPs. Date Sat Dec 18 10:30:45 2010    Stringency: medium Mate Pair: 1 Read Length: 50 Polymorphism Rate: 0.003000 Bayes Coverage: 60 Bayes_Single_SNP: 1 Filter_Single_SNP: 1 Quick_P_Threshold: 0.997000 Bayes_P_Threshold: 0.040000 Minimum_Allele_Ratio: 0.150000 Minimum_Allele_Ratio_Multiple_of_Dicolor_Error: 100
	##1     chr1
	##2     chr2
	##3     chr3
	##4     chr4
	##5     chr5
	##6     chr6
	##7     chr7
	##8     chr8
	##9     chr9
	##10    chr10
	##11    chr11
	##12    chr12
	##13    chr13
	##14    chr14
	##15    chr15
	##16    chr16
	##17    chr17
	##18    chr18
	##19    chr19
	##20    chr20
	##21    chr21
	##22    chr22
	##23    chrX
	##24    chrY
	##25    chrM
	# source-version SOLiD BioScope diBayes(SNP caller)
	#Chr    Source  Type    Pos_Start       Pos_End Score   Strand  Phase   Attributes
	chr1    SOLiD_diBayes   SNP     221367  221367  0.091151        .       .       genotype=R;reference=G;coverage=3;refAlleleCounts=1;refAlleleStarts=1;refAlleleMeanQV=29;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=27;diColor1=11;diColor2=33;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     555317  555317  0.095188        .       .       genotype=Y;reference=T;coverage=13;refAlleleCounts=11;refAlleleStarts=10;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=29;diColor1=00;diColor2=22;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     555327  555327  0.037582        .       .       genotype=Y;reference=T;coverage=12;refAlleleCounts=6;refAlleleStarts=6;refAlleleMeanQV=19;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=29;diColor1=12;diColor2=30;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     559817  559817  0.094413        .       .       genotype=Y;reference=T;coverage=9;refAlleleCounts=5;refAlleleStarts=4;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=14;diColor1=11;diColor2=33;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     714068  714068  0.000000        .       .       genotype=M;reference=C;coverage=13;refAlleleCounts=7;refAlleleStarts=6;refAlleleMeanQV=25;novelAlleleCounts=6;novelAlleleStarts=4;novelAlleleMeanQV=22;diColor1=00;diColor2=11;het=1;flag= 
	The file conforms to standard GFF3 specifications, but the last column is solid-
	specific and it gives certain parameters for the SNP calls.

An example of the short indel format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version SOLiD Corona Lite v.4.0r2.0, find-small-indels.pl v 1.0.1, process-small-indels v 0.2.2, 2009-01-12 12:28:49
	##type DNA
	##date 2009-01-26
	##time 18:33:20
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files ../../mp-results/JOAN_20080104_1.pas,../../mp-results/BARB_20071114_1.pas,../../mp-results/BARB_20080227_2.pas
	##run-path /data/results2/Yoruban-frag-indel/try.01.06/mp-w2x25-2x-4x-8x-10x/2x
	##Filter-settings: max-ave-read-pos=none,min-ave-from-end-pos=9.1,max-nonreds-4filt=2,min-insertion-size=none,min-deletion-size=none,max-insertion-size=none,max-deletion-size=none,require-called-indel-size?=T
	chr1    AB_SOLiD Small Indel Tool       deletion        824501  824501  1       .       .       del_len=1;tight_chrom_pos=824501-824502;loose_chrom_pos=824501-824502;no_nonred_reads=2;no_mismatches=1,0;read_pos=4,6;from_end_pos=21,19;strands=+,-;tags=R3,F3;indel_sizes=-1,-1;read_seqs=G3021212231123203300032223,T3321132212120222323222101;dbSNP=rs34941678,chr1:824502-824502(-),EXACT,1,/GG
	chr1    AB_SOLiD Small Indel Tool       insertion_site  1118641 1118641 1       .       .       ins_len=3;tight_chrom_pos=1118641-1118642;loose_chrom_pos=1118641-1118642;no_nonred_reads=2;no_mismatches=0,1;read_pos=17,6;from_end_pos=8,19;strands=+,+;tags=F3,R3;indel_sizes=3,3;read_seqs=T0033001100022331122033112,G3233112203311220000001002

The keyword deletion or insertion_site is used in the fourth column to indicate 
that file format.

An example of the medium CNV format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version SOLiD Corona Lite v.4.0r2.0, find-small-indels.pl v 1.0.1, process-small-indels v 0.2.2, 2009-01-12 12:28:49
	##type DNA
	##date 2009-01-27
	##time 15:54:36
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files big_d20e5-del12n_up-ConsGrp-2nonred.pas.sum
	##run-path /data/results2/Yoruban-frag-indel/try.01.06/mp-results-lmp-e5/big_d20e5-indel_950_2050
	chr1    AB_SOLiD Small Indel Tool       deletion        3087770 3087831 1       .       .       del_len=62;tight_chrom_pos=none;loose_chrom_pos=3087768-3087773;no_nonred_reads=2;no_mismatches=2,2;read_pos=27,24;from_end_pos=23,26;strands=-,+;tags=F3,F3;indel_sizes=-62,-62;read_seqs=T11113022103331111130221213201111302212132011113022,T02203111102312122031111023121220311111333012203111
	chr1    AB_SOLiD Small Indel Tool       deletion        4104535 4104584 1       .       .       del_len=50;tight_chrom_pos=4104534-4104537;loose_chrom_pos=4104528-4104545;no_nonred_reads=3;no_mismatches=0,4,4;read_pos=19,19,27;from_end_pos=31,31,23;strands=+,+,-;tags=F3,R3,R3;indel_sizes=-50,-50,-50;read_seqs=T31011011013211110130332130332132110110132020312332,G21031011013211112130332130332132110132132020312332,G20321302023001101123123303103303101113231011011011
	chr1    AB_SOLiD Small Indel Tool       insertion_site  2044888 2044888 1       .       .       ins_len=18;tight_chrom_pos=2044887-2044888;loose_chrom_pos=2044887-2044889;no_nonred_reads=2;bead_ids=1217_1811_209,1316_908_1346;no_mismatches=0,2;read_pos=13,15;from_end_pos=37,35;strands=-,-;tags=F3,F3;indel_sizes=18,18;read_seqs=T31002301231011013121000101233323031121002301231011,T11121002301231011013121000101233323031121000101231;non_indel_no_mismatches=3,1;non_indel_seqs=NIL,NIL
	chr1    AB_SOLiD Small Indel Tool       insertion_site  74832565        74832565        1       .       .       ins_len=16;tight_chrom_pos=74832545-74832565;loose_chrom_pos=74832545-74832565;no_nonred_reads=2;bead_ids=1795_181_514,1651_740_519;no_mismatches=0,2;read_pos=13,13;from_end_pos=37,37;strands=-,-;tags=F3,R3;indel_sizes=16,16;read_seqs=T33311111111111111111111111111111111111111111111111,G23311111111111111111111111111111111111111311011111;non_indel_no_mismatches=1,0;non_indel_seqs=NIL,NIL

An example of the large indel format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version ???
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files /data/results5/yoruban_strikes_back_large_indels/LMP/five_mm_unique_hits_no_rescue/5_point_6x_del_lib_1/results/NA18507_inter_read_indels_5_point_6x.dat
	##run-path 
	chr1    AB_SOLiD Large Indel Tool       insertion_site  1307279 1307791 1       .       .       deviation=-742;stddev=7.18;ref_clones=-;dev_clones=4
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2042742 2042861 1       .       .       deviation=-933;stddev=8.14;ref_clones=-;dev_clones=3
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2443482 2444342 1       .       .       deviation=-547;stddev=11.36;ref_clones=-;dev_clones=17
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2932046 2932984 1       .       .       deviation=-329;stddev=6.07;ref_clones=-;dev_clones=14
	chr1    AB_SOLiD Large Indel Tool       insertion_site  3166925 3167584 1       .       .       deviation=-752;stddev=13.81;ref_clones=-;dev_clones=14

An example of the CNV format by GFF3-SOLiD if given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version ???
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files Yoruban_cnv.coords
	##run-path 
	chr1    AB_CNV_PIPELINE repeat_region   1062939 1066829 .       .       .       fraction_mappable=51.400002;logratio=-1.039300;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   1073630 1078667 .       .       .       fraction_mappable=81.000000;logratio=-1.409500;copynum=1;numwindows=2
	chr1    AB_CNV_PIPELINE repeat_region   2148325 2150352 .       .       .       fraction_mappable=98.699997;logratio=-1.055000;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   2245558 2248109 .       .       .       fraction_mappable=78.400002;logratio=-1.042900;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   3489252 3492632 .       .       .       fraction_mappable=59.200001;logratio=-1.119900;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   5654415 5657276 .       .       .       fraction_mappable=69.900002;logratio=1.114500;copynum=4;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   9516165 9522726 .       .       .       fraction_mappable=65.850006;logratio=-1.316700;numwindows=2
	chr1    AB_CNV_PIPELINE repeat_region   16795117        16841025        .       .       .       fraction_mappable=44.600002;logratio=1.880778;copynum=7;numwindows=9

The keyword repeat_region is used here, although it actually refers to CNVs.

An example of the inversion format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.2
	##generated by SOLiD inversion tool
	chr10   AB_SOLiD        inversion       46443107        46479585        268.9   .       .       left=chr10:46443107-46443146;right=chr10:46479583-46479585;leftscore=295.0;rightscore=247.0;count_AAA_further_left=117;count_AAA_left=3;count_AAA_right=3;count_AAA_further_right=97;left_min_count_AAA=chr10:46443107-46443112;count_AAA_min_left=0;count_AAA_max_left=3;right_min_count_AAA=chr10:46479585-46479585;count_AAA_min_right=1;count_AAA_max_right=3;homozygous=UNKNOWN
	chr4    AB_SOLiD        inversion       190822813       190850112       214.7   .       .       left=chr4:190822813-190822922;right=chr4:190850110-190850112;leftscore=140.0;rightscore=460.0;count_AAA_further_left=110;count_AAA_left=78;count_AAA_right=74;count_AAA_further_right=77;left_min_count_AAA=chr4:190822813-190822814;count_AAA_min_left=69;count_AAA_max_left=77;right_min_count_AAA=chr4:190850110-190850112;count_AAA_min_right=74;count_AAA_max_right=74;homozygous=NO
	chr6    AB_SOLiD        inversion       168834969       168837154       175.3   .       .       left=chr6:168834969-168835496;right=chr6:168836643-168837154;leftscore=185.4;rightscore=166.2;count_AAA_further_left=67;count_AAA_left=43;count_AAA_right=40;count_AAA_further_right=59;left_min_count_AAA=chr6:168835058-168835124,chr6:168835143-168835161,chr6:168835176-168835181,chr6:168835231-168835262;count_AAA_min_left=23;count_AAA_max_left=29;right_min_count_AAA=chr6:168836643-168836652;count_AAA_min_right=23;count_AAA_max_right=31;homozygous=NO

The program should be able to recognize all the above GFF3-SOLiD format 
automatically, and handle them accordingly.

=item * B<Complete Genomics format>

This format is provided by the Complete Genomics company to their customers. The 
file var-[ASM-ID].tsv.bz2 includes a description of all loci where the assembled 
genome differs from the reference genome.

An example of the Complete Genomics format is shown below:

	#BUILD  1.5.0.5
	#GENERATED_AT   2009-Nov-03 19:52:21.722927
	#GENERATED_BY   dbsnptool
	#TYPE   VAR-ANNOTATION
	#VAR_ANN_SET    /Proj/Pipeline/Production_Data/REF/HUMAN-F_06-REF/dbSNP.csv
	#VAR_ANN_TYPE   dbSNP
	#VERSION        0.3
	
	>locus  ploidy  haplotype       chromosome      begin   end     varType reference       alleleSeq       totalScore      hapLink xRef
	1       2       all     chr1    0       959     no-call =       ?                       
	2       2       all     chr1    959     972     =       =       =                       
	3       2       all     chr1    972     1001    no-call =       ?                       
	4       2       all     chr1    1001    1008    =       =       =                       
	5       2       all     chr1    1008    1114    no-call =       ?                       
	6       2       all     chr1    1114    1125    =       =       =                       
	7       2       all     chr1    1125    1191    no-call =       ?                       
	8       2       all     chr1    1191    1225    =       =       =                       
	9       2       all     chr1    1225    1258    no-call =       ?                       
	10      2       all     chr1    1258    1267    =       =       =                       
	12      2       all     chr1    1267    1275    no-call =       ?                       
	13      2       all     chr1    1275    1316    =       =       =                       
	14      2       all     chr1    1316    1346    no-call =       ?                       
	15      2       all     chr1    1346    1367    =       =       =                       
	16      2       all     chr1    1367    1374    no-call =       ?                       
	17      2       all     chr1    1374    1388    =       =       =                       
	18      2       all     chr1    1388    1431    no-call =       ?                       
	19      2       all     chr1    1431    1447    =       =       =                       
	20      2       all     chr1    1447    1454    no-call =       ?                       

The following information is provided in documentation from Complete Genomics, that describes the var-ASM format.

	1. locus. Identifier of a particular genomic locus
	2. ploidy. The ploidy of the reference genome at the locus (= 2 for autosomes, 2 for pseudoautosomal regions on the sex chromosomes, 1 for males on the non-pseudoautosomal parts of the sex chromosomes, 1 for mitochondrion, '?' if varType is 'no-ref' or 'PAR-called-in-X'). The reported ploidy is fully determined by gender, chromosome and location, and is not inferred from the sequence data.
	3. haplotype. Identifier for each haplotype at the variation locus. For diploid genomes, 1 or 2. Shorthand of 'all' is allowed where the varType field is one of 'ref', 'no-call', 'no-ref', or 'PAR-called-in-X'. Haplotype numbering does not imply phasing; haplotype 1 in locus 1 is not necessarily in phase with haplotype 1 in locus 2. See hapLink, below, for phasing information.
	4. chromosome. Chromosome name in text: 'chr1','chr2', ... ,'chr22','chrX','chrY'. The mitochondrion is represented as 'chrM'. The pseudoautosomal regions within the sex chromosomes X and Y are reported at their coordinates on chromosome X.
	5. begin. Reference coordinate specifying the start of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	6. end. Reference coordinate specifying the end of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	7. varType. Type of variation, currently one of:
		snp: single-nucleotide polymorphism
		ins: insertion
		del: deletion
		sub: Substitution of one or more reference bases with the bases in the allele column
		'ref' : no variation; the sequence is identical to the reference sequence on the indicated haplotype
		no-call-rc: 'no-call reference consistent 'one or more bases are ambiguous, but the allele is potentially consistent with the reference
		no-call-ri: 'no-call reference inconsistent' one or more bases are ambiguous, but the allele is definitely inconsistent with the reference
		no-call: an allele is completely indeterminate in length and composition, i.e. alleleSeq = '?'
		no-ref: the reference sequence is unspecified at this locus.
		PAR-called-in-X: this locus overlaps one of the pseudoautosomal regions on the sex chromosomes. The called sequence is reported as diploid sequence on Chromosome X; on chromosome Y the sequence is reported as varType = 'PAR-called-in-X'.
	8. reference. The reference sequence for the locus of variation. Empty when varType is ins. A value of '=' indicates that the user must consult the reference for the sequence; this shorthand is only used in regions where no haplotype deviates from the reference sequence.
	9. alleleSeq. The observed sequence at the locus of variation. Empty when varType is del. '?' isused to indicate 0 or more unknown bases within the sequence; 'N' is used to indicate exactly one unknown base within the sequence.'=' is used as shorthand to indicate identity to the reference sequence for non-variant sequence, i.e. when varType is 'ref'.
	10. totalScore. A score corresponding to a single variation and haplotype, representing the confidence in the call.
	11. hapLink. Identifier that links a haplotype at one locus to haplotypes at other loci. Currently only populated for very proximate variations that were assembled together. Two calls that share a hapLink identifier are expected to be on the same haplotype,
	12. xRef. Field containing external variation identifiers, currently only populated for variations corroborated directly by dbSNP. Format: dbsnp:[rsID], with multiple entries separated by the semicolon (;).

In older versions of the format specification, the sub keyword used to be insdel 
keyword. ANNOVAR takes care of this.

=item * B<SOAPsnp format>

An example of the SOAP SNP caller format is shown below:

	chr8  35782  A  R  1  A  27  1  2  G  26  1  2  5   0.500000  2.00000  1  5   
	chr8  35787  G  R  0  G  25  4  6  A  17  2  4  10  0.266667  1.60000  0  5   

The following information is provided in documentation from BGI who developed 
SOAP suite. It differs slightly from the description at the SOAPsnp website, and 
presumably the website is outdated.

	Format description:(left to right)
	1. Chromosome name
	2. Position of locus
	3. Nucleotide at corresponding locus of reference sequence
	4. Genotype of sequencing sample
	5. Quality value
	6. nucleotide with the highest probability(first nucleotide)
	7. Quality value of the nucleotide with the highest probability
	8. Number of supported reads that can only be aligned to this locus 
	9. Number of all supported reads that can be aligned to this locus
	10. Nucleotide with higher probability 
	11. Quality value of nucleotide with higher probability 
	12. Number of supported reads that can only be aligned to this locus 
	13. Number of all supported reads that can be aligned to this locus 
	14. Total number of reads that can be aligned to this locus 
	15. Order and quality value
	16. Estimated copy number for this locus 
	17. Presence of this locus in the dbSNP database. 1 refers to presence and 0 refers to inexistence
	18. The distance between this locus and another closest SNP

=item * B<SOAPindel format>

The current version of ANNOVAR handles SoapSNP and SoapIndel automatically via a 
single argument '--format soap'. An example of SOAP indel caller format is shown 
below:

	chr11   44061282        -       +2      CT      Hete
	chr11   45901572        +       +1      C       Hete
	chr11   48242562        *       -3      TTC     Homo
	chr11   57228723        *       +4      CTTT    Homo
	chr11   57228734        *       +4      CTTT    Homo
	chr11   57555685        *       -1      C       Hete
	chr11   61482191        -       +3      TCC     Hete
	chr11   64608031        *       -1      T       Homo
	chr11   64654936        *       +1      C       Homo
	chr11   71188303        +       -1      T       Hete
	chr11   75741034        +       +1      T       Hete
	chr11   76632438        *       +1      A       Hete
	chr11   89578266        *       -2      AG      Homo
	chr11   104383261       *       +1      T       Hete
	chr11   124125940       +       +4      CCCC    Hete
	chr12   7760052 *       +1      T       Homo
	chr12   8266049 *       +3      ACG     Homo

I do not see a documentation describing this format yet as of September 2010.

=item B<--SOAPsv format>

An example is given below:

	Chr2 Deletion 42894 43832 43167 43555 388 0-0-0 FR 41

An explanation of the structural variation format is given below:

	Format description (from left to right)
	1. Chromosome name
	2. Type of structure variation
	3. Minimal value of start position in cluster
	4. Maximal value of end position in cluster
	5. Estimated start position of this structure variation
	6. Estimated end position of this structure variation
	7. Length of SV
	8. Breakpoint of SV (only for insertion)
	9. Unusual matching mode (F refers to align with forward sequence, R refers
	to align with reverse
	sequence)
	10. number of paired-end read which support this structure variation

=item * B<MAQ format>

MAQ can perform alignment and generate genotype calls, including SNP calls and 
indel calls. The format is described below:

For indel header: The output is TAB delimited with each line consisting of chromosome, start 
position, type of the indel, number of reads across the indel, size of the indel 
and inserted/deleted nucleotides (separated by colon), number of indels on the 
reverse strand, number of indels on the forward strand, 5' sequence ahead of the 
indel, 3' sequence following the indel, number of reads aligned without indels 
and three additional columns for filters.

An example is below:

	chr10   110583  -       2       -2:AG   0       1       GCGAGACTCAGTATCAAAAAAAAAAAAAAAAA        AGAAAGAAAGAAAAAGAAAAAAATAGAAAGAA        1       @2,     @72,   @0,
	chr10   120134  -       8       -2:CA   0       1       CTCTTGCCCGCTCACACATGTACACACACGCG        CACACACACACACACACATCAGCTACCTACCT        7       @65,62,61,61,45,22,7,   @9,12,13,13,29,52,67,   @0,0,0,0,0,0,0,
	chr10   129630  -       1       -1:T    1       0       ATGTTGTGACTCTTAATGGATAAGTTCAGTCA        TTTTTTTTTAGCTTTTAACCGGACAAAAAAAG        0       @       @      @
	chr10   150209  -       1       4:TTCC  1       0       GCATATAGGGATGGGCACTTTACCTTTCTTTT        TTCCTTCCTTCCTTCCTTCCCTTTCCTTTCCT        0       @       @      @
	chr10   150244  -       2       -4:TTCT 0       1       CTTCCTTCCTTCCTTCCCTTTCCTTTCCTTTC        TTCTTTCTTTCTTTCTTTCTTTTTTTTTTTTT        0       @       @      @
	chr10   159622  -       1       3:AGG   0       1       GAAGGAGGAAGGACGGAAGGAGGAAGGAAGGA        AGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGA        0       @       @      @
	chr10   206372  -       2       2:GT    1       0       ATAATAGTAACTGTGTATTTGATTATGTGTGC        GTGTGTGTGTGTGTGTGTGTGTGTGCGTGCTT        1       @37,    @37,   @8,
	chr10   245751  -       11      -1:C    0       1       CTCATAAATACAAGTCATAATGAAAGAAATTA        CCACCATTTTCTTATTTTCATTCATTTTTAGT        10      @69,64,53,41,30,25,22,14,5,4,   @5,10,21,33,44,49,52,60,69,70,  @0,0,0,0,0,0,0,0,0,0,
	chr10   253066  -       1       2:TT    0       1       TATTGATGAGGGTGGATTATACTTTAGAACAC        TATTCAAACAGTTCTTCCACATATCTCCCTTT        0       @       @      @
	chr10   253455  -       2       -3:AAA  1       0       GTTGCACTCCAGCCTGGCGAGATTCTGTCTCC        AAAAAAAAAAAAAAAAATTGTTGTGAAATACA        1       @55,    @19,   @4,

For snp output file: Each line consists of chromosome, position, reference base, 
consensus base, Phred-like consensus quality, read depth, the average number of 
hits of reads covering this position, the highest mapping quality of the reads 
covering the position, the minimum consensus quality in the 3bp flanking regions 
at each side of the site (6bp in total), the second best call, log likelihood 
ratio of the second best and the third best call, and the third best call.

An example is below:

	chr10   83603   C       T       28      12      2.81    63      34      Y       26      C
	chr10   83945   G       R       59      61      4.75    63      62      A       47      G
	chr10   83978   G       R       47      40      3.31    63      62      A       21      G
	chr10   84026   G       R       89      22      2.44    63      62      G       49      A
	chr10   84545   C       T       54      9       1.69    63      30      N       135     N
	chr10   85074   G       A       42      5       1.19    63      38      N       108     N
	chr10   85226   A       T       42      5       1.00    63      42      N       107     N
	chr10   85229   C       T       42      5       1.00    63      42      N       112     N
	chr10   87518   A       G       39      4       3.25    63      38      N       9       N
	chr10   116402  T       C       39      4       1.00    63      38      N       76      N


=item * B<CASAVA format>

An example of Illumina CASAVA format is given below:

	#position       A       C       G       T       modified_call   total   used    score           reference       type
	14930   3       0       8       0       GA      11      11      29.10:11.10             A       SNP_het2
	14933   4       0       7       0       GA      11      11      23.50:13.70             G       SNP_het1
	14976   3       0       8       0       GA      11      11      24.09:9.10              G       SNP_het1
	15118   2       1       4       0       GA      8       7       10.84:6.30              A       SNP_het2

An example of the indels is given below:

	# ** CASAVA depth-filtered indel calls **
	#$ CMDLINE /illumina/pipeline/install/CASAVA_v1.7.0/libexec/CASAVA-1.7.0/filterIndelCalls.pl--meanReadDepth=2.60395068970547 --indelsCovCutoff=-1 --chrom=chr1.fa /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0000.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0001.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0002.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0003.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0004.txt
	#$ CHROMOSOME chr1.fa
	#$ MAX_DEPTH undefined
	#
	#$ COLUMNS pos CIGAR ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) max2_gtype bp1_reads ref_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count
	948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
	978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
	1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
	1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0

=item * B<VCF4 format>

VCF4 can be used to describe both population-level variation information, or for 
reads derived from a single individual.

One example of the indel format for one individual is given below:

	##fileformat=VCFv4.0
	##IGv2_bam_file_used=MIAPACA2.alnReAln.bam
	##INFO=<ID=AC,Number=2,Type=Integer,Description="# of reads supporting consensus indel/any indel at the site">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="total coverage at the site">
	##INFO=<ID=MM,Number=2,Type=Float,Description="average # of mismatches per consensus indel-supporting read/per reference-supporting read">
	##INFO=<ID=MQ,Number=2,Type=Float,Description="average mapping quality of consensus indel-supporting reads/reference-supporting reads">
	##INFO=<ID=NQSBQ,Number=2,Type=Float,Description="Within NQS window: average quality of bases from consensus indel-supporting reads/from reference-supporting reads">
	##INFO=<ID=NQSMM,Number=2,Type=Float,Description="Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
	##INFO=<ID=SC,Number=4,Type=Integer,Description="strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
	##IndelGenotyperV2=""
	##reference=hg18.fa
	##source=IndelGenotyperV2
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Miapaca_trimmed_sorted.bam      
	chr1    439     .       AC      A       .       PASS    AC=5,5;DP=7;MM=7.0,3.0;MQ=23.4,1.0;NQSBQ=23.98,25.5;NQSMM=0.04,0.0;SC=2,3,0,2   GT      1/0
	chr1    714048  .       T       TCAAC   .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.266666,21.932203;NQSMM=0.0,0.15254237;SC=3,0,3,3 GT      0/1
	chr1    714049  .       G       GC      .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.233334,21.83051;NQSMM=0.0,0.15254237;SC=3,0,3,3  GT      0/1
	chr1    813675  .       A       AATAG   .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=25.74,25.166666;NQSMM=0.0,0.033333335;SC=4,1,1,2       GT      0/1
	chr1    813687  .       AGAGAGAGAGAAG   A       .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=24.54,25.2;NQSMM=0.02,0.06666667;SC=4,1,1,2    GT      1/0


=back

The code was written by Dr. Kai Wang and modified by Dr. Germán Gastón Leparc. 
Various users have provided sample input files for many SNP callin software, for 
the development of conversion subroutines. We thank these users for their 
continued support to improve the functionality of the script.

For questions or comments, please contact kai@openbioinformatics.org.

=cut
