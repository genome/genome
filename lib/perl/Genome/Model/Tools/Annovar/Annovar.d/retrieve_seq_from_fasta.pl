#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

our $VERSION = 			'$Revision: 499 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-02-23 17:47:27 -0800 (Thu, 23 Feb 2012) $';

our ($verbose, $help, $man);
our ($regionfile);
our ($outfile, $format, $seqdir, $seqfile, $tabout, $altchr);
our $localtime = localtime;

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'seqdir=s'=>\$seqdir, 'seqfile=s'=>\$seqfile,
	'tabout'=>\$tabout, 'altchr'=>\$altchr) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

my $path = $0;
$path =~ s/[^\\\/]+$//;
$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in
$verbose and print STDERR "NOTICE: Adding $path into system executable path\n";

($regionfile) = @ARGV;
if (not defined $format) {
	$format = 'tab';
	print STDERR "NOTICE: The default regionfile format is set as 'tab'\n";
}
if (not defined $outfile) {
	$outfile = "$regionfile.fa";
	print STDERR "NOTICE: The output file is written to $outfile (use --outfile to change this)\n";
}
open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";

$format =~ m/^(refGene|knownGene|ensGene|genericGene|tab|simple)$/ or pod2usage ("Error in argument: the --format argument can be only 'refGene', 'knownGene', 'ensGene', 'genericGene', 'tab' or 'simple'");

my ($allregion, $allchr, $sorted_region) = readAllRegion ($regionfile, $format);

my %seqhash;				#database sequence for each chromosome
my %name_seq;				#sequence for each region
my (%seqlen, %discordlen, %badorf);	#store the length of each sequence, and the ID of sequences with discordant length, ORF contains stop codon
my ($count_success, @failure) = (0);
for my $curchr (sort keys %$allchr) {
	my ($seqid, $curseq) = ('', '');
	my $fastafile;
	
	if ($seqfile) {			#if the user specify the actual sequence file, the program will read directly from this file
		$fastafile = $seqfile;
	} else {			#otherwise, the program will read from chr1.fa, chr2.fa, chr3.fa, etc.
		%seqhash = ();		#clear the seqhash storage
		
		if ($curchr =~ m/^chr/) {
			$fastafile = File::Spec->catfile($seqdir, "$curchr.fa");		#by default, all FASTA files should be saved at fastadir, with the same name
		} else {
			$fastafile = File::Spec->catfile($seqdir, "chr$curchr.fa");
		}
		
		if (not -e $fastafile) {
			$fastafile = File::Spec->catfile($seqdir, "softMask", "$curchr.fa");	#rheMac2 has this directory configuration
		}
		
		if (not -e $fastafile) {
			if ($curchr =~ m/^chr/) {
				my $curchr1 = $curchr;
				$curchr1 =~ s/^chr//;
				$curchr1 =~ s/_random$//;		#the chr6_random.fa file will be saved in the 6/ directory
				$fastafile = File::Spec->catfile($seqdir, $curchr1, "$curchr.fa");	#panTro2_seq has this directory configuration
			}
		}
		if (not -e $fastafile) {			#to handle cases where no "chr" prefix is given
			print STDERR "WARNING: the FASTA file $curchr.fa cannot be retrieved from the specified directory $seqdir. Sequences in this chromosome will not be processed\n";
			next;
		}
	}
	
	if (not %seqhash) {
		open (FASTA, $fastafile) or print STDERR "WARNING: cannot read from FASTA file $fastafile so sequences in $curchr will not be processed: $!\n" and next;
		while (<FASTA>) {
			if (m/^>(\S+)/) {
				$seqid and $seqhash{$seqid} = $curseq;			#finish reading the sequence for seqid and save it
				$seqid = $1;
				$curseq = '';
			} else {
				s/[\r\n]+$//;
				$curseq .= uc $_;					#only use upper case characters
			}
		}
		$seqhash{$seqid} = $curseq;
	}
	
	if (not $seqhash{$curchr}) {			#this chromosome just do not have FASTA sequences (maybe users used a wrong seqdir
		print STDERR "WARNING: Unable to retrieve regions at $curchr due to lack of sequence information\n";
		next;
	}
	
	print STDERR "NOTICE: Finished reading ", scalar keys %seqhash, " sequences from $fastafile\n";
	
	for my $i (0 .. @{$allregion->{$curchr}}-1) {
		my ($name, $start, $end, $strand, $exonpos, $mrnastart, $mrnaend) = @{$allregion->{$curchr}[$i]};
		my @start = split (/,/, $start);
		my @end = split (/,/, $end);
		my $seq;
		for my $i (0..@start-1){
			if ($start[$i] >= length ($seqhash{$curchr})) {		#here there must be an annotation error in user-specified gene/region definition file
				print STDERR "WARNING: Ignoring the start position start=$start[$i] since it is longer than the $curchr sequence (length=" , length($seqhash{$curchr}), ")\n";
				undef $seq;
				last;
			}
			$seq .= substr ($seqhash{$curchr}, $start[$i], $end[$i]-$start[$i]);
		}
		if ($strand eq '-' and defined $seq) {
			$seq = reverse $seq;
			$seq =~ tr/acgtACGT/tgcaTGCA/;
		}
		if (defined $seq) {
			if (defined $seqlen{$name}) {
				$seqlen{$name} != length ($seq) and $verbose and warn "WARNING: the sequence $name was found more than once with different sequence lengths\n";
				$seqlen{$name} != length ($seq) and $discordlen{$name}++;
			} else {
				$seqlen{$name} = length ($seq);
			}			
			#if ($tabout) {
			#	print $name, "\t", $seq, "\n";
			#} else {
			#	if ($seqdir) {
			#		print ">$name ", $discordlen{$name}?"Warning: $name occur more than once in file with discordant length. ":"", "Comment: this sequence (leftmost exon at $curchr:$start[0]) is generated by ANNOVAR on $localtime, based on regions speficied in $regionfile and sequence files stored at $seqdir.\n$seq\n";
			#	} else {
			#		print ">$name ", $discordlen{$name}?"Warning: $name occur more than once in file with discordant length. ":"", "Comment: this sequence (leftmost exon at $curchr:$start[0]) is generated by ANNOVAR on $localtime, based on regions speficied in $regionfile and the sequence file $seqfile.\n$seq\n";
			#	}
			#}
			$name_seq{$name, $exonpos} = $seq;		#added 2011feb19
			$count_success++;
			
			if (defined $mrnastart and defined $mrnaend) {
				my $dna = substr ($seq, $mrnastart-1, $mrnaend-$mrnastart+1);
				my $protein = translateDNA ($dna);
				if ($protein =~ m/\*\w/) {
					$verbose and print STDERR "WARNING: Potential FASTA sequence error for $name (ORF with premature stop codon)!!!\n";
					$badorf{$name, $exonpos}++;
				}
			}
		} else {
			print STDERR "WARNING: DNA sequence for $name cannot be inferred\n";
			push @failure, $name;
		}
	}
}

for my $i (0 .. @$sorted_region-1) {
	my ($name, $exonpos) = @{$sorted_region->[$i]};
	if (not $name_seq{$name, $exonpos}) {
		print STDERR "WARNING: Cannot identify sequence for $name (starting from $exonpos)\n";
		next;
	}
	
	if ($tabout) {
		print $name, "\t", $name_seq{$name, $exonpos}, "\n";
	} else {
		print ">$name ", $discordlen{$name}?"Warning: $name occur more than once in file with discordant length. ":" ", $badorf{$name, $exonpos}?"Warning: this coding transcript does not have correct ORF annotation. ":" ", 
			"Comment: this sequence (leftmost exon at $exonpos) is generated by ANNOVAR on $localtime, based on regions speficied in $regionfile and ", $seqdir?"sequence files stored at $seqdir.\n":"the sequence file $seqfile.\n",
			$name_seq{$name, $exonpos}, "\n";
	}
}
	
	

print STDERR "NOTICE: Finished writting FASTA for $count_success genomic regions to $outfile\n";
if (%discordlen) {
	my @discordlen = keys %discordlen;
	@discordlen = splice (@discordlen, 0, 5);
	print STDERR "WARNING: ${\(scalar keys %discordlen)} regions occur more than once with discordant sequence length (for example, ", join(", ", @discordlen), ")\n";
}
if (%badorf) {
	my @badorf = keys %badorf;
	@badorf = splice (@badorf, 0, 5);
	print STDERR "WARNING: ${\(scalar keys %badorf)} gene regions do not have complete ORF (for example, ", join(", ", @badorf), ")\n";
}
if (@failure) {
	print STDERR "WARNING: DNA sequences for ${\(scalar @failure)} regions cannot be inferred (see $outfile.failure for identifiers).\n";
	open (FAILURE, ">$outfile.failure") or die "Error: cannot write to failure file $outfile.failure: $!\n";
	print FAILURE join ("\n", @failure), "\n";
}

sub translateDNA {
	my ($seq) = @_;
	my ($nt3, $protein);
	my %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");

	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		if (defined $codon1{$1}) {
			$protein .= $codon1{$1};
		} else {
			$verbose and print STDERR "WARNING: invalid triplets found in DNA sequence to be translated: <$1> within <$seq>\n";
			$protein .= "X";
		}
	}
	return $protein;
}

sub readAllRegion {
	my ($regionfile, $format) = @_;
	my (@field, $name, $curchr, $start, $end, $strand);
	my (%allregion, %allchr, @sorted_region);		#allregion is a hash with chr as key; allchr has chr as key; sorted_region is an array with name and linecount as elements
	my ($countregion, %counthap, $exonpos) = (0);
	open (REGION, $regionfile) or die "Error: cannot read from regionfile $regionfile: $!\n";
	print STDERR "NOTICE: Reading region file $regionfile ... ";
	while (<REGION>) {
		s/[\r\n]+$//;
		$strand = '+';
		my ($mrnastart, $mrnaend);
		if ($format eq 'knownGene' or $format eq 'refGene' or $format eq 'ensGene' or $format eq 'genericGene') {
			#example: 1475    NM_021648       chr6    -       116677823       116681954       116680619       116681864       1       116677823,      116681954,      0       TSPYL4 cmpl     cmpl    0,
			@field = split (/\t/, $_);
			if ($format eq 'knownGene') {
				@field >= 11 or die "Error: invalid record found in '$format' file $regionfile (>=11 fields expected but found ${\(scalar @field)} fields): <$_>\n";
			} elsif ($format eq 'refGene') {
				@field == 16 or die "Error: invalid record found in '$format' file $regionfile (16 fields expected): <$_>\n";
				shift @field;		#shifting out the first field
			} elsif ($format eq 'ensGene') {
				if (@field == 16) {		#many organisms (yeast, etc) do not have ENS-prefixed gene names in ENSEMBL
					shift @field;
				} elsif (@field == 10) {
					1;
				} else {
					die "Error: invalid record found in '$format' file $regionfile (16 or 10 fields expected but found ${\(scalar @field)}): <$_>\n";
				}
			} elsif ($format eq 'genericGene') {
				@field >= 11 or die "Error: invalid record found in '$format' file $regionfile (>=11 fields expected but found ${\(scalar @field)} fields): <$_>\n";
				shift @field;
			}
						
			$name = $field[0];
			($curchr, $start, $end) = @field[1,8,9];
			$strand = $field[2];
			$strand eq '+' or $strand eq '-' or die "Error: invalid strand designations in the region file $regionfile: <$_>\n";
			
			if (not $altchr and $curchr =~ m/hap\d+$/) {		#for example, chr5_h2_hap1, chr6_cox_hap1
				$counthap{$curchr}++;
				next;
			}
			$exonpos = "$curchr:$field[3]";
			
			#the following paragraph was added 2011Nov12 to identify transcripts without complete coding information
			my ($txstart, $txend, $cdsstart, $cdsend) = @field[3, 4, 5, 6];
			if ($cdsstart != $cdsend) {			#this is a real protein coding gene by annotation
				my @exonstart = split (/,/, $start);
				my @exonend = split (/,/, $end);
				
				$txstart++;
				$cdsstart++;
				@exonstart = map {$_+1} @exonstart;
				if ($strand eq '+') {
					#<---->-----<--->----<------>----<----->-----<--->
					#             **********
					my $intron = 0;
					for my $i (0 .. @exonstart-1) {
						$i and $intron += ($exonstart[$i]-$exonend[$i-1]-1);
						if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
							$mrnastart = $cdsstart-$txstart+1-$intron;
						}
						if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
							$mrnaend = $cdsend-$txstart+1-$intron;
						}
						
					}
				} elsif ($strand eq '-') {
					#<---->-----<--->----<------>----<----->-----<--->
					#             **********
					my $intron = 0;
					for (my $i=@exonstart-1; $i>=0; $i--) {
						$i<@exonstart-1 and $intron += ($exonstart[$i+1]-$exonend[$i]-1);
						if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
							$mrnastart = $txend-$cdsend+1-$intron;
						}
						if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
							$mrnaend = $txend-$cdsstart+1-$intron;
						}
						
					}
				}
			}
			
			
		} elsif ($format eq 'tab') {
			#example: 1	1000	2000
			#example: chr1	1000	2000
			@field = split (/\t/, $_);
			@field >= 3 or die "Error: invalid record found in input file (at least 3 tab-delimited fields expected for '-format tab'): <$_>\n";
			$field[1] =~ m/^\d+$/ and $field[2] =~ m/^\d+$/ or die "Error: invalid record found in input file (2nd and 3rd column should be positive integer for '-format tab': <$_> <$field[1]><$field[2]\n";
			$name = "$field[0]:$field[1]-$field[2]";
			($curchr, $start, $end) = @field[0..2];
			if (not $curchr =~ m/^chr/) {
				$curchr = "chr$curchr";
			}
			
			$start--;		#make zero-start coordinate, to be consistent with UCSC
			$exonpos = "$curchr:$start";
		} elsif ($format eq 'simple') {
			#example: 1:1000-2000
			#example: chr1:1000-2000
			m/^((chr)?\w+):(\d+)\-(\d+)/ or die "Error: invalid record found in input line (chr:start-end expected): <$_>\n";
			($curchr, $start, $end) = ($1, $3, $4);
			if (not $curchr =~ m/^chr/) {
				$curchr = "chr$curchr";
			}
			$name = "$curchr:$start-$end";
			
			$start--;		#make zero-start coordinate, to be consistent with UCSC
			$exonpos = "$curchr:$start";
		}
		
		$allchr{$curchr}++;
		if (defined $mrnastart and defined $mrnaend) {
			push @{$allregion{$curchr}}, [$name, $start, $end, $strand, $exonpos, $mrnastart, $mrnaend];
		} else {
			push @{$allregion{$curchr}}, [$name, $start, $end, $strand, $exonpos];
		}
		push @sorted_region, [$name, $exonpos];		#these two elements uniquely identifies a sequence ID
		$countregion++;
	}
	print STDERR "Done with $countregion regions from ", scalar (keys %allchr), " chromosomes";
	if (%counthap) {
		my @counthap = keys %counthap;
		if (@counthap > 5) {
			@counthap = splice (@counthap, 0, 3);
			print STDERR " (${\(scalar keys %counthap)} 'alternative haplotype' chromosomes are skipped, such as ", join (", ", @counthap), ", etc)\n";
		} else {
			print STDERR " (${\(scalar keys %counthap)} 'alternative haplotype' chromosomes are skipped, including ", join (", ", @counthap), ")\n";
		}
	} else {
		print STDERR "\n";
	}
	return (\%allregion, \%allchr, \@sorted_region);
}


=head1 SYNOPSIS

 retrieve_seq_from_fasta.pl [arguments] <regionfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --format <string>		format of the regionfile
 	    --outfile <file>		output file name (default: input.fa)
 	    --seqdir <directory>	directory containing sequence files for genomes
 	    --seqfile <file>		sequence file name
 	    --tabout			use tab-delimited output (default is FASTA file)
 	    --altchr			process alternative haplotype chromosome (default: OFF)

 Function: reformat sequences at specific genomic positions from whole-genome FASTA files

 Example: retrieve_seq_from_fasta.pl -format tab -seqdir humandb/ region.tabinput
 	  retrieve_seq_from_fasta.pl -format tab -seqfile chrall.fa region.tabinput
          retrieve_seq_from_fasta.pl -format simple -seqdir humandb/ region.input
          retrieve_seq_from_fasta.pl -format refGene -seqdir humandb hg18_refGenet.txt
          retrieve_seq_from_fasta.pl -format genericGene -seqdir humandb/ hg18_gencodeGene.txt
          

 Version: $LastChangedDate: 2012-02-23 17:47:27 -0800 (Thu, 23 Feb 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--format>

specify the format of the input file (default: tab). Other supported formats 
include simple, refGene, knownGene, ensGene, genericGene and many other format 
support will be added in the future.

=item B<--outfile>

specify the output file name (default: <input>.fa)

=item B<--seqdir>

specify the directory of the sequence files (usually these sequence files are 
organized by chromosomes, and that the file name represent the chromosome).

=item B<--seqfile>

specify the sequence file name (assuming that all sequences are stored within a 
single file).

=item B<--tabout>

use tab-delimited output. By default, FASTA file will be printed out.

=item B<--altchr>

process alternative haplotype chromosome, that look like chr6_cox_hap1. By 
default this will not be processed.


=back

=head1 DESCRIPTION

This program is designed for retrieving sequences at specific genomic regions 
from whole-genome FASTA files. These files should be saved at a specific 
directory, with one chromosome per file.

The FASTA files can be easily downloaded automatically by the ANNOVAR software:

	annotate_variation.pl -buildver hg18 -downdb seq humandb/

The FASTA files (one for each chromosome) will be saved at the humandb/ 
directory. Users can then use the retrieve_seq_from_fasta.pl program to retrieve 
sequences from specific genomic regions in these chromosomes.

Several region input files can be supported, via the -format argument:

=over 8

=item * B<simple>

The simple format is something like chr1:200-3000, one region per line.

=item * B<tab>

the tab format is tab-delimited format, with first 3 columns as chr, start, end, 
at each line.

=item * B<refgene>

the refgene format is provided by the UCSC Genome Browser, in the refGene.txt 
database annotation files.

=item * B<knowngene>

the knowngene format is provided by the UCSC Genome Browser, in the knownGene.txt 
database annotation files.

=item * B<ensgene>

the ensgene format is provided by the UCSC Genome Browser, in the ensGene.txt 
database annotation files.

=back

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut