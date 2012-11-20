#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 502 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-03-08 08:43:06 -0800 (Thu, 08 Mar 2012) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $step, $skip, $remove, $model, $checkfile, $dispensable, $ver1000g, $mceway, $verdbsnp, $genetype, $maf_threshold);
our ($file1000g);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'step=s'=>\$step, 'skip'=>\$skip, 'remove'=>\$remove, 'model=s'=>\$model,
	'checkfile!'=>\$checkfile, 'dispensable=s'=>\$dispensable, 'ver1000g=s'=>\$ver1000g, 'mceway=i'=>\$mceway, 'verdbsnp=i'=>\$verdbsnp, 'genetype=s'=>\$genetype,
	'maf_threshold=f'=>\$maf_threshold) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

my $path = $0;
$path =~ s/[^\\\/]+$//;
$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in

($queryfile, $dbloc) = @ARGV;

$outfile ||= $queryfile;
$genetype ||= 'refgene';
$genetype =~ m/^refgene|knowngene|ensgene$/i or pod2usage ("Error in argument: the --genetype can be 'refgene', 'knowngene' or 'ensgene' only");

if (not defined $buildver) {
	$buildver = 'hg18';
	print STDERR "NOTICE: the --buildver argument is set as 'hg18' by default\n";
}


if (defined $model) {
	$model =~ m/^(recessive|dominant)$/ or pod2usage ("Error in argument: the --model argument can be only 'recessive' or 'dominant' (for autosomes or X-linked recessive disease in males");
} else {
	$model = 'recessive';
	print STDERR "NOTCE: the --model argument is set as 'recessive' by default\n";
}

if (not defined $ver1000g) {
	if ($buildver eq 'hg18') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g' by default\n";
		$ver1000g = '1000g';
	} elsif ($buildver eq 'hg19') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g2010nov' by default\n";
		$ver1000g = '1000g2010nov';
	}
}

if (defined $maf_threshold) {
	$maf_threshold >= 0 and $maf_threshold <= 1 or pod2usage ("Error: the --maf_threshold argument must be between 0 and 1 inclusive");
	$maf_threshold = " --maf_threshold $maf_threshold";
} else {
	$maf_threshold = '';
}

if ($ver1000g eq '1000g') {
	$file1000g = '2009_04';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010') {
	$file1000g = '2010_03';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010jul') {
	$file1000g = '2010_07';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010nov') {
	$file1000g = '2010_11';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g eq '1000g2011may') {
	$file1000g = '2011_05';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g eq '1000g2012feb') {
	$file1000g = '2012_02';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	$file1000g = $1 . '_' . $monthhash{$2};
} else {
	pod2usage ("Error in argument: the --ver1000g $ver1000g is not yet supported by this program");
}

if ($buildver eq 'hg18') {
	defined $mceway or $mceway = 44;
} elsif ($buildver eq 'hg19') {
	defined $mceway or $mceway = 46;
} else {
	defined $mceway or pod2usage ("Error in argument: please specify the --mceway argument (the default values work on hg18/hg19 builds but not $buildver)");
}

defined $verdbsnp or $verdbsnp = 130;

not defined $checkfile and $checkfile = 1;

my %valistep;
if ($step) {
	my @step = split (/,/, $step);
	for my $i (0 .. @step-1) {
		if ($step[$i] =~ m/^(\d+)-(\d+)$/) {
			for my $nextstep ($1 .. $2) {
				$valistep{$nextstep}++;
			}
		} elsif ($step[$i] =~ m/^(\d+)$/) {
			$valistep{$1}++;
		} else {
			pod2usage ("Error: invalid -step argument ($step) is specified. Please use comma-separated number only (dash line such as 1-5 is accepted)");
		}
	}
} else {
	for my $nextstep (1..9) {
		$valistep{$nextstep}++;
	}
}

if (not $valistep{1}) {
	if (not $skip) {
		pod2usage ("Error in argument: the step 1 in the procedure cannot be skipped");
	}
}

if ($valistep{10}) {
	$dispensable or pod2usage ("Error in argument: the --dispensable file (one gene name per line) must be supplied for executing step 10");
}
$checkfile and checkFileExistence ();


my $sc;
my $linecount;

#run step 1
if ($valistep{1}) {
	$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $genetype -outfile $outfile.step1 $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 1 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	system ("fgrep -v -w synonymous $outfile.step1.exonic_variant_function | fgrep -v -w nonframeshift | cut -f 4- > $outfile.step2.varlist");
	system ("fgrep -w splicing $outfile.step1.variant_function | cut -f 3- >> $outfile.step2.varlist");
	system ("sort $outfile.step2.varlist | uniq > $outfile.step2.varlist.temp; mv $outfile.step2.varlist.temp $outfile.step2.varlist");
	

} else {
	$skip or system ("cp $queryfile $outfile.step2.varlist") and die "Error cp";
}
$remove and unlink ("$outfile.step1.varlist");
$linecount = qx/cat $outfile.step2.varlist | wc -l/; chomp $linecount;
$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
	
#run step2
if ($valistep{2}) {
	$sc = "annotate_variation.pl -regionanno -dbtype mce${mceway}way -buildver $buildver -outfile $outfile.step2 $outfile.step2.varlist $dbloc";
	print STDERR "\nNOTICE: Running step 2 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	system ("cut -f 3- $outfile.step2.${buildver}_phastConsElements${mceway}way > $outfile.step3.varlist");
} else {
	$skip or system ("cp $outfile.step2.varlist $outfile.step3.varlist") and die "Error cp";
}
$remove and unlink ("$outfile.step2.varlist", "$outfile.step2.${buildver}_phastConsElements${mceway}way");
$linecount = qx/cat $outfile.step3.varlist | wc -l/; chomp $linecount;
$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;


#run step3
if ($valistep{3}) {

	$sc = "annotate_variation.pl -regionanno -dbtype segdup -buildver $buildver -outfile $outfile.step3 $outfile.step3.varlist $dbloc";
	print STDERR "\nNOTICE: Running step 3 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	system ("cut -f 4- $outfile.step3.${buildver}_genomicSuperDups > $outfile.step4.varlist1; fgrep -v -f $outfile.step4.varlist1 $outfile.step3.varlist > $outfile.step4.varlist");
	unlink ("$outfile.step4.varlist1");
} else {
	$skip or system ("cp $outfile.step3.varlist $outfile.step4.varlist") and die "Error cp";
}
$remove and unlink ("$outfile.step3.varlist", "$outfile.step4.varlist1", "$outfile.step3.${buildver}_segdup");
$linecount = qx/cat $outfile.step4.varlist | wc -l/; chomp $linecount;
$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;

if ($buildver eq 'hg19') {
	if ($valistep{4} or $valistep{5} or $valistep{6}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_all $maf_threshold -buildver $buildver -outfile $outfile.step4 $outfile.step4.varlist $dbloc";
		print STDERR "\nNOTICE: Running step 4/5/6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
		system ("cp $outfile.step4.${buildver}_ALL.sites.${file1000g}_filtered $outfile.step7.varlist");
	} else {
		$skip or system ("cp $outfile.step4.varlist $outfile.step7.varlist") and die "Error cp";
	}
	$remove and unlink ("$outfile.step4.varlist", "$outfile.step4.${buildver}_${ver1000g}_ceu_filtered", "$outfile.step4.${buildver}_${ver1000g}_ceu_reported");
	$linecount = qx/cat $outfile.step7.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
} else {

	#run step4
	if ($valistep{4}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_ceu $maf_threshold -buildver $buildver -outfile $outfile.step4 $outfile.step4.varlist $dbloc";
		print STDERR "\nNOTICE: Running step 4 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
		system ("cp $outfile.step4.${buildver}_CEU.sites.${file1000g}_filtered $outfile.step5.varlist");
	} else {
		$skip or system ("cp $outfile.step4.varlist $outfile.step5.varlist") and die "Error cp";
	}
	$remove and unlink ("$outfile.step4.varlist", "$outfile.step4.${buildver}_${ver1000g}_ceu_filtered", "$outfile.step4.${buildver}_${ver1000g}_ceu_reported");
	$linecount = qx/cat $outfile.step5.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
	
	
	#run step5
	if ($valistep{5}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_yri $maf_threshold -buildver $buildver -outfile $outfile.step5 $outfile.step5.varlist $dbloc";
		print STDERR "\nNOTICE: Running step 5 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
		system ("cp $outfile.step5.${buildver}_YRI.sites.${file1000g}_filtered $outfile.step6.varlist");
	} else {
		$skip or system ("cp $outfile.step5.varlist $outfile.step6.varlist") and die "Error cp";
	}
	$remove and unlink ("$outfile.step5.varlist", "$outfile.step5.${buildver}_${ver1000g}_yri_filtered", "$outfile.step5.${buildver}_${ver1000g}_yri_reported");
	$linecount = qx/cat $outfile.step6.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
	
	#run step6
	if ($valistep{6}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_jptchb $maf_threshold -buildver $buildver -outfile $outfile.step6 $outfile.step6.varlist $dbloc";
		print STDERR "\nNOTICE: Running step 6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
		system ("cp $outfile.step6.${buildver}_JPTCHB.sites.${file1000g}_filtered $outfile.step7.varlist");
	} else {
		$skip or system ("cp $outfile.step6.varlist $outfile.step7.varlist") and die "Error cp";
	}
	$remove and unlink ("$outfile.step6.varlist", "$outfile.step6.${buildver}_${ver1000g}_jptchb_filtered", "$outfile.step6.${buildver}_${ver1000g}_jptchb_reported");
	$linecount = qx/cat $outfile.step7.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
}

#run step7
if ($valistep{7}) {
	$sc = "annotate_variation.pl -filter -dbtype snp${verdbsnp} -buildver $buildver -outfile $outfile.step7 $outfile.step7.varlist $dbloc";
	print STDERR "\nNOTICE: Running step 7 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	system ("cp $outfile.step7.${buildver}_snp${verdbsnp}_filtered $outfile.step8.varlist");
} else {
	$skip or system ("cp $outfile.step7.varlist $outfile.step8.varlist") and die "Error cp";
}
$remove and unlink ("$outfile.step7.varlist", "$outfile.step7.${buildver}_snp${verdbsnp}_filtered", "$outfile.step7.${buildver}_${verdbsnp}_reported");
$linecount = qx/cat $outfile.step8.varlist | wc -l/; chomp $linecount;
$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;

#run step8 (this cannot be skipped any more as of Mar 2011)
if (1) {
	$sc = "fgrep -f $outfile.step8.varlist $outfile.step1.exonic_variant_function | cut -f 2- > $outfile.step8;";		#function, gene name, plus original input
	$sc .= "cut -f 3- $outfile.step8 > $outfile.step8.temp;";		#avinput
	$sc .= "fgrep -v -f $outfile.step8.temp $outfile.step8.varlist > $outfile.step8.temp1;";	#splicing variants
	$sc .= "fgrep -f $outfile.step8.temp1 $outfile.step1.variant_function >> $outfile.step8;";	#splicing functions
	print STDERR "\nNOTICE: Running step 8 with system command <$sc>\n";
	system ($sc);			#this command may generate error, because the $outfile.step8.temp1 file may be empty
	system ("cp $outfile.step8 $outfile.step9.varlist");
} else {
	$skip or system ("cp $outfile.step8.varlist $outfile.step9.varlist") and die "Error cp";
}
$remove and unlink ("$outfile.step8.varlist", "$outfile.step8.temp", "$outfile.step8.temp1", "$outfile.step1.exonic_variant_function");
$linecount = qx/cat $outfile.step9.varlist | wc -l/; chomp $linecount;
$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;

#run step9, which cannot be skipped
if (1) {
	my %found;
	open (VAR, "$outfile.step9.varlist") or die "Error: cannot read from varlist file $outfile.step9.varlist: $!\n";
	while (<VAR>) {
		my @field = split (/\t/, $_);
		$field[1] =~ s/,$//;
		my @gene = split (/,/, $field[1]);
		my %gene;
		for my $i (0 .. @gene-1) {
			$gene[$i] =~ m/^(\w+)/ or die "Error: invalid record in input file $outfile.step9.varlist (gene name expected): <$_>\n";
			$gene{$1}++;
		}
		for my $key (keys %gene) {
			$found{$key}++;
			if (m/\bhom\b/) {		#if the word "hom" is printed in the input line, add the count by 1.
				$found{$key}++;
			}
		}
	}
	
	print STDERR "\nNOTICE: a list of potentially important genes and the number of variants in them are written to $outfile.genelist\n";
	open (OUT, ">$outfile.genelist") or die "Error: cannot write to output file $outfile.genelist: $!\n";
	for my $key (keys %found) {
		if ($model eq 'recessive') {
			if ($found{$key} >= 2) {
				print OUT "$key\t$found{$key}\n";
			}
		} elsif ($model eq 'dominant') {
			if ($found{$key} >= 1) {
				print OUT "$key\t$found{$key}\n";
			}
		}
	}
	print STDERR "NOTICE: Consider filter out the list of dispensable genes from the $outfile.genelist file to identify the final candidate gene list.\n";
}

#run step10
if ($valistep{10}) {
	if ($dispensable) {
		$sc = "fgrep -v -w -f $dispensable $outfile.genelist > $outfile.genelist.nondispensible";
		print STDERR "\nNOTICE: Running step 10 with system command <$sc>\n";
		system ($sc);
	} else {
		print STDERR "ERROR: Step 10 requires a dispensible file as input to filter out non-important genes from $outfile.genelist\n";
	}
}


sub checkFileExistence {
	my @file = ("${buildver}_refGene.txt", "${buildver}_refLink.txt", "${buildver}_refGeneMrna.fa", "${buildver}_genomicSuperDups.txt", 
		"${buildver}_snp$verdbsnp.txt", );
	if ($buildver eq 'hg18') {
		push @file, "${buildver}_phastConsElements44way.txt";
		push @file, "${buildver}_CEU.sites.${file1000g}.txt", "${buildver}_YRI.sites.${file1000g}.txt", "${buildver}_JPTCHB.sites.${file1000g}.txt";
	} elsif ($buildver eq 'hg19') {
		push @file, "${buildver}_phastConsElements46way.txt";
		push @file, "${buildver}_ALL.sites.${file1000g}.txt";
	}
	for my $i (0 .. @file-1) {
		my $dbfile = File::Spec->catfile ($dbloc, $file[$i]);
		-f $dbfile or die "Error: the required database file $dbfile does not exist. Please download it via -downdb argument by annotate_variation.pl.\n";
	}
}


=head1 SYNOPSIS

 auto_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --outfile <string>		output file name prefix
            --buildver <string>		genome build version (default: hg18)
            --step <string>		comma-delimited steps (by default, run all steps)
            --skip			skip the first a few steps before execution of the pipeline
            --remove			remove all temporary files
            --model <string>		dominant or recessive model
            --<no>checkfile		check to see whether DB file exist before execution
            --dispensable <file>	specify a file containing list of dispensable genes (one per line)
            --ver1000g <string>		version of 1000G data (1000g, 1000g2010, 1000g2010jul, 1000g2010nov, 1000g2011may)
            --mceway <int>		number of species in MCE alignment (default: 44 for hg18, 46 for hg19)
            --verdbsnp <int>		dbSNP version (default: 130)
            --genetype <string>		gene definition (default: refgene)
            --maf_threshold <float>	MAF threshold in 1000G data

 Function: automatically run a pipeline on a list of variants (potentially 
 whole-genome SNPs from a patient with Mendelian disease) and identify a small 
 subset that are most likely causal for Mendelian diseases
 
 Example: auto_annovar.pl ex2.human humandb/
 
 Version: $LastChangedDate: 2012-03-08 08:43:06 -0800 (Thu, 08 Mar 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--step>

specify a few specific steps to perform, rather than running all steps from 
scratch. A comma-delimited list should be supplied, and dash can be correctly 
recognized. For example, 1-5,7 means step 1 through step 5 plus step 7.

=item B<--skip>

skip file preparation step. By default, if a step is not specified, a file 
should still be prepared for that step. Using --skip argument omit this process, 
so that one can resume a previously interrupted auto_annovar session.

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--model>

specify whether dominant or recessive model is used. The input file MUST contain 
het or hom status for this argument to work properly.

=item B<--checkfile>

check to make sure that database files exist, before executing the current 
program.

=item B<--dispensable>

specify a list of dispensable genes for use in step 10. The file contains one 
gene name per line.

=item B<--ver1000g>

specify the version of the 1000 Genomes Project data set to be used in 
annotation. For -buildver hg18, default is 1000g, but values such as 1000g2010 
and 1000g2010jul can be alterntaively specified. For -buildver hg19, default is 
1000g2010nov.

=item B<--mceway>

specify the multi-way alignment to be used in conserved region annotation.

=item B<--verdbsnp>

specify the dbSNP version to be used in annotation. Default is 130.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

Thea aut_annovar script automate a "variants reduction" procedure that tries to 
identify a subset of most likely causal variants from a subject with exome-
sequencing or whole-genome sequencing data, assuming that the subject has a rare 
Mendelian disease caused by non-synonymous mutations. These steps are described 
below:

=item B<step1>

identify splicing and exonic variants, but remove synonymous mutations (these 
subsets of variants are more likely to be functional.

=item B<step2>

identify the variants located in genomic regions that are annotated as Most 
Conserved Elements (more likely to be functional)

=item B<step3>

identify variants that are not located in segmental duplications regions (less 
likely to be affected by genotype calling issues)

=item B<step4/step5/step6>

remove variants observed in the 1000 Genomes Project (CEU, YRI, JPTCHB, 
respectively) because these variants are less likely to be causing Mendelian 
diseases. Howeve, I do see MANY known mendelian disease variants in 1000G, so 
some MAF threshold should really be applied here!

=item B<step7>

remove variants observed in the dbSNP 130 because these variants are 
less likely to be causing Mendelian diseases. Howeve, I do see MANY known 
mendelian disease variants in dbSNP130!

=item B<step8>

collect the remaining variants, map them to genes, and prepare the next step

=item B<step9>

identify a list of genes that are likely to be affected by deleterious mutations

=item B<step10; by default OFF>

remove some genes that are regarded as "dispensable" 

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut