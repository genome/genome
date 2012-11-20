#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 511 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-11-01 10:38:05 -0700 (Thu, 01 Nov 2012) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $step, $checkfile, $remove, $verdbsnp, $ver1000g, $genetype, $veresp, $alltranscript);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'step=s'=>\$step, 
	'checkfile!'=>\$checkfile, 'remove'=>\$remove, 'verdbsnp=i'=>\$verdbsnp, 'ver1000g=s'=>\$ver1000g, 'genetype=s'=>\$genetype,
	'veresp=i'=>\$veresp, 'alltranscript!'=>\$alltranscript) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

my $path = $0;
$path =~ s/[^\\\/]+$//;
$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in

($queryfile, $dbloc) = @ARGV;

$outfile ||= $queryfile;
$buildver ||= 'hg18';
$buildver eq 'hg18' or $buildver eq 'hg19' or pod2usage ("Error in argument: the --buildver argument can be 'hg18' and 'hg19' only");
not defined $checkfile and $checkfile = 1;

defined $verdbsnp or $verdbsnp = 130;
$genetype ||= 'refgene';
$genetype =~ m/^refgene|knowngene|ensgene|gencodegene$/i or $genetype =~ m/wgEncodeGencode\w+$/ or pod2usage ("Error in argument: the --genetype can be 'refgene', 'knowngene', 'ensgene' or 'gencodegene' only");

my ($file1000g);
if (not defined $ver1000g) {
	if ($buildver eq 'hg18') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g' by default\n";
		$ver1000g = '1000g';
	} elsif ($buildver eq 'hg19') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g2010nov' by default\n";
		$ver1000g = '1000g2010nov';
	}
}

if (not defined $veresp) {
	$veresp = 5400;
}

if ($genetype eq 'gencodegene') {
	if ($buildver eq 'hg18') {
		$genetype = 'wgEncodeGencodeManualV3';
	} elsif ($buildver eq 'hg19') {
		$genetype = 'wgEncodeGencodeManualV4';
	}
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
} elsif ($ver1000g =~ m/^1000g(20\d\d)([a-z]{3})$/) {
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	$file1000g = $1 . '_' . $monthhash{$2};
} else {
	pod2usage ("Error in argument: the --ver1000g $ver1000g is not yet supported by this program");
}

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
	for my $nextstep (1..20) {
		$valistep{$nextstep}++;
	}
}

$checkfile and checkFileExistence ();


my $sc;

#run step 1
if ($valistep{1}) {
	$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $genetype -outfile $outfile -exonsort $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 1 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step2
if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		$sc = "annotate_variation.pl -regionanno -dbtype mce44way -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 2 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	} elsif ($buildver eq 'hg19') {
		$sc = "annotate_variation.pl -regionanno -dbtype mce46way -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
}

#run step3
if ($valistep{3}) {
	$sc = "annotate_variation.pl -regionanno -dbtype segdup -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 3 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}


if ($buildver eq 'hg19') {
	if ($valistep{4} or $valistep{5} or $valistep{6}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 4/5/6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
} 
else {
	#run step4
	if ($valistep{4}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_ceu -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 4 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step5
	if ($valistep{5}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_yri -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 5 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step6
	if ($valistep{6}) {
		$sc = "annotate_variation.pl -filter -dbtype ${ver1000g}_jptchb -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
}

#run step7
if ($valistep{7}) {
	$sc = "annotate_variation.pl -filter -dbtype snp$verdbsnp -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 7 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step8
if ($valistep{8}) {
	$sc = "annotate_variation.pl -filter -dbtype avsift -buildver $buildver -sift 0 -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 8 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step9 to step13
if ($valistep{9} or $valistep{10} or $valistep{11} or $valistep{12} or $valistep{13}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_all -buildver $buildver -outfile $outfile $queryfile $dbloc -otherinfo";
	print STDERR "\nNOTICE: Running step 9-13 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step14
if ($valistep{14}) {
	$sc = "annotate_variation.pl -filter -dbtype esp${veresp}_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 14 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

=head1
#run step9
if ($valistep{9}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_pp2 -score_threshold 0 -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 9 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step10
if ($valistep{10}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_phylop -score_threshold 0 -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 10 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step11
if ($valistep{11}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_mt -score_threshold 0 -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 11 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step12
if ($valistep{12}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_lrt -score_threshold 0 -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 12 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step13
if ($valistep{13}) {
	$sc = "annotate_variation.pl -filter -dbtype ljb_gerp++ -score_threshold 0 -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 12 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}
=cut


open (FUNCTION, "$outfile.variant_function") or die "Error: cannot read from variant function file: $!\n";

if ($valistep{1}) {
	open (STEP1, "$outfile.exonic_variant_function") or die "Error: cannot read from exonic variant function file: $!\n";
}

if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		open (STEP2, "$outfile.hg18_phastConsElements44way") or die "Error: cannot read from mce file: $!\n";
	} 
	elsif ($buildver eq 'hg19') {
		open (STEP2, "$outfile.hg19_phastConsElements46way") or die "Error: cannot read from mce file: $!\n";
	}
}

if ($valistep{3}) {
	open (STEP3, "$outfile.${buildver}_genomicSuperDups") or die "Error: cannot read from segdup file: $!\n";
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg19_ALL.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.${file1000g}_dropped: $!\n";
	}
} 
else {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg18_CEU.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_CEU.sites.${file1000g}_dropped: $!\n";
	}
	if ($valistep{5}) {
		open (STEP5, "$outfile.hg18_YRI.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_YRI.sites.${file1000g}_dropped: $!\n";
	}
	if ($valistep{6}) {
		open (STEP6, "$outfile.hg18_JPTCHB.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_JPTCHB.sites.${file1000g}_dropped: $!\n";
	}
}

if ($valistep{7}) {
	open (STEP7, "$outfile.${buildver}_snp${verdbsnp}_dropped") or die "Error: cannot read from snp$verdbsnp drop file: $!\n";
}

if ($valistep{8}) {
	open (STEP8, "$outfile.${buildver}_avsift_dropped") or die "Error: cannot read from avsift drop file: $!\n";
}

if ($valistep{9} or $valistep{10} or $valistep{11} or $valistep{12} or $valistep{13}) {
	open (STEP9, "$outfile.${buildver}_ljb_all_dropped") or die "Error: cannot read from ljb_all drop file: $!\n";
}

if ($valistep{14}) {
	open (STEP14, "$outfile.${buildver}_esp${veresp}_all_dropped") or die "Error: cannot read from esp${veresp}_all drop file: $!\n";
}



my (@allstep);
for my $i (1 .. 3) {
	$allstep[$i] = [];
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^1000g\w+_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
} 
else {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^1000g\w*_ceu\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}

	if ($valistep{5}) {
		while (<STEP5>) {
			m/^1000g\w*_yri\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 5 : <$_>\n";
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	
	if ($valistep{6}) {
		while (<STEP6>) {
			m/^1000g\w*_jptchb\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 6 : <$_>\n";
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
}

if ($valistep{7}) {
	while (<STEP7>) {
		m/^snp\d+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 7 : <$_>\n";
		$allstep[7]->{$2} = $1;
	}
}

if ($valistep{8}) {
	while (<STEP8>) {
		m/^avsift\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 8 : <$_>\n";
		$allstep[8]->{$2} = $1;
	}
}

if ($valistep{9} or $valistep{10} or $valistep{11} or $valistep{12} or $valistep{13}) {
	while (<STEP9>) {
		m/^ljb_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 9 : <$_>\n";
		$allstep[9]->{$2} = $1;
	}
}

if ($valistep{14}) {
	while (<STEP14>) {
		m/^esp\d+_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 8 : <$_>\n";
		$allstep[14]->{$2} = $1;
	}
}



print STDERR "NOTICE: Finished loading filterstep database file\n";

open (OUT, ">$outfile.genome_summary.csv") or die "Error: cannot write to output file: $!\n";
open (OUTE, ">$outfile.exome_summary.csv") or die "Error: cannot write to output file: $!\n";

if ($buildver eq 'hg19') {
	print OUT join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP${veresp}_ALL", "${ver1000g}_ALL", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP${veresp}_ALL", "${ver1000g}_ALL", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
} else {
	print OUT join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP${veresp}_ALL", "${ver1000g}_CEU,${ver1000g}_YRI,${ver1000g}_JPTCHB", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join (',', qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP${veresp}_ALL", "${ver1000g}_CEU,${ver1000g}_YRI,${ver1000g}_JPTCHB", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
}

while (<FUNCTION>) {
	s/[\r\n]+$//;
	m/^(\S+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)(.*)/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
	my ($function, $gene, $varstring, $otherinfo) = ($1, $2, $3, $4||'');
	my $exonic;
	if ($function =~ m/\bsplicing\b/ or $function =~ m/\bexonic\b/) {
		$exonic = 1;
	}
	print OUT qq/"$function","$gene"/;
	$exonic and print OUTE qq/"$function","$gene"/;
	
	if (not @{$allstep[1]}) {
		if ($valistep{1}) {
			if (defined ($_ = <STEP1>)) {
				m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				my ($efun, $aachange, $varstring) = ($1, $2, $3);
				my @aachange = split (/:|,/, $aachange);
				if ($alltranscript) {
					push @{$allstep[1]}, $varstring, $efun, $aachange;
				} else {
					if (@aachange >= 5) {
						push @{$allstep[1]}, $varstring, $efun, "$aachange[1]:$aachange[3]:$aachange[4]";
					} else {
						push @{$allstep[1]}, $varstring, $efun, $aachange;		#aachange could be "UNKNOWN"
					}
				}
			}
		}
	}
	
	if (not @{$allstep[2]}) {
		if ($valistep{2}) {
			if (defined ($_ = <STEP2>)) {
				m/^mce\d+way\tScore=(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				push @{$allstep[2]}, $2, $1;
			}
		}
	}
	
	if (not @{$allstep[3]}) {
		if ($valistep{3}) {
			if (defined ($_ = <STEP3>)) {
				m/^segdup\tScore=(\S+);\S+\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 3 : <$_>\n";
				push @{$allstep[3]}, $2, ($1>0.01)?sprintf("%.2f", $1):$1;
			}
		}
	}
	
	for my $i (1 .. 3) {			#starting from step 1 to step 7
		my $curstep = $allstep[$i];
		if (@$curstep and $curstep->[0] eq $varstring) {
			if ($i == 1) {
				print OUT qq/,"$curstep->[1]"/, qq/,"$curstep->[2]"/;
				$exonic and print OUTE qq/,"$curstep->[1]"/, qq/,"$curstep->[2]"/;
			} else {
				print OUT qq/,"$curstep->[1]"/;
				$exonic and print OUTE qq/,"$curstep->[1]"/;
			}
			@$curstep = ();
		}
		else {
			if ($i == 1) {
				print OUT ",,";
				$exonic and print OUTE ",,";
			} else {
				print OUT ",";
				$exonic and print OUTE ",";
			}
		}
	}
	
	if (defined $allstep[14]->{$varstring}) {
		print OUT qq/,$allstep[14]->{$varstring}/;
		$exonic and print OUTE qq/,$allstep[14]->{$varstring}/;
	} else {
		print OUT ",";
		$exonic and print OUTE ",";
	}
	
	for my $i (4 .. 8) {		#step9 already includes step 9 through step 13
		if ($buildver eq 'hg19') {
			$i == 5 and next;
			$i == 6 and next;
		}
		if (defined $allstep[$i]->{$varstring}) {
			print OUT qq/,$allstep[$i]->{$varstring}/;
			$exonic and print OUTE qq/,$allstep[$i]->{$varstring}/;
		} else {
			print OUT ",";
			$exonic and print OUTE ",";
		}
	}
	if (defined $allstep[9]->{$varstring}) {
		print OUT qq/,$allstep[9]->{$varstring}/;
		$exonic and print OUTE qq/,$allstep[9]->{$varstring}/;
	} else {
		print OUT ",,,,,,,,,,,";
		$exonic and print OUTE ",,,,,,,,,,,";
	}
	
	
	
	my @varstring = split (/\s+/, $varstring);
	$otherinfo =~ s/^\s+//;
	my @otherinfo = split (/\t/, $otherinfo);
	for my $i (0 .. @otherinfo-1) {
		$otherinfo[$i] = qq/"$otherinfo[$i]"/;
	}

	print OUT ',', join (',', @varstring), ",", join (',', @otherinfo), "\n";
	$exonic and print OUTE ',', join (',', @varstring), ",", join (',', @otherinfo), "\n";
}

print STDERR "NOTICE: Final whole-genome summary was written to $outfile.genome_summary.csv file\n";
print STDERR "NOTICE: Final whole-exome summary was written to $outfile.exome_summary.csv file\n";

if ($remove) {
	unlink ("$outfile.variant_function", "$outfile.exonic_variant_function", "$outfile.hg18_phastConsElements44way", "$outfile.hg19_phastConsElements46way", 
		"$outfile.${buildver}_genomicSuperDups", 
		"$outfile.${buildver}_esp${veresp}_all_dropped",
		"$outfile.hg19_ALL.sites.${file1000g}_dropped", 
		"$outfile.hg18_CEU.sites.${file1000g}_dropped", "$outfile.hg18_YRI.sites.${file1000g}_dropped", "$outfile.hg18_JPTCHB.sites.${file1000g}_dropped", 
		"$outfile.${buildver}_snp${verdbsnp}_dropped", "$outfile.${buildver}_avsift_dropped", "$outfile.${buildver}_ljb_all_dropped");
}

sub checkFileExistence {
	my @file = ("${buildver}_refGene.txt", "${buildver}_refLink.txt", "${buildver}_refGeneMrna.fa", "${buildver}_genomicSuperDups.txt", 
		"${buildver}_snp$verdbsnp.txt", "${buildver}_avsift.txt", "${buildver}_ljb_all.txt",  "${buildver}_esp${veresp}_all.txt");
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

 summarize_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
		-h, --help			print help message
		-m, --man			print complete documentation
		-v, --verbose			use verbose output
		    --outfile <string>		output file name prefix
		    --buildver <string>		genome build version (default: hg18)
		    --remove			remove all temporary files
		    --verdbsnp <int>		dbSNP version to use (default: 130)
		    --ver1000g <string>		1000G version (default: 1000g2010nov)
		    --veresp <int>		ESP version (default: 5400)
		    --genetype <string>		gene definition can be refgene (default), knowngene, ensgene
		    --checkfile			check existence of database file (default: ON)
		    --(no)alltranscript		output all transcript names for exonic variants (default: OFF)

 Function: automatically run a pipeline on a list of variants and summarize 
 their functional effects in a comma-delimited file, to be opened by Excel for 
 manual filtering
 
 Example: summarize_annovar.pl ex2.human humandb/
 
 Version: $LastChangedDate: 2012-11-01 10:38:05 -0700 (Thu, 01 Nov 2012) $

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

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--verdbsnp>

version of dbSNP to use in annotation. By default, 130 is used.

=item B<--ver1000g>

version of 1000 Genomes Project dataset to use in annotation. By default, 1000g2010nov is used.

=item B<--veresp>

version of ESP dataset to use in annotation. By default, 5400 is used, but users can also choose 6500.

=item B<--genetype>

gene definition systems to use, including refgene (default), knowngene, ensgene

=item B<--checkfile>

check to make sure that database files exist, before executing the current 
program.

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

summarize_annovar is a script that automate some routines in ANNOVAR and 
generates an Excel-compatible file for users to manually browse and filter.

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut
