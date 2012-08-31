#!/usr/bin/env genome-perl

use strict;
use English;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Carp;

my ( $locus_name, $config_template, $help );
GetOptions(
	"locus-name=s" => \$locus_name,
	"template=s"   => \$config_template,
	"help"         => \$help,
) or die "Incorrect usage!\n";

pod2usage( -verbose => 2 ) if ($help);
pod2usage( -verbose => 2 ) unless ( defined($locus_name) );

## We will set the $todays_date
my ( $day, $month, $year ) = (localtime)[ 3, 4, 5 ];
my $todays_date = sprintf( "%02d%02d%02d", $year + 1900, $month + 1, $day, ) . getppid;

$config_template =
  "/gscuser/ssurulir/test-genome/BIFBRE2011TST/BIFBRE2011TST_template_config_date"
  unless ($config_template);

$locus_name =~ m/(.*)TST/;
my $locus_id = $1;

my ( $file_name, $dir ) = fileparse($config_template);
my $new_file_name = $locus_name . "_config_" . $todays_date;

## We will write to $new_file_name
open(WRITE, "> $dir$new_file_name")
  || die "Unable to open $dir.$new_file_name to write: $!\n";

## We will parse the file into hash
open CONFIG,
  "<" . $config_template || die "Error opening file: " . $OS_ERROR;

# Read the file
my @param_list = <CONFIG>;
chomp @param_list;
close CONFIG;

my %config;
foreach my $record (@param_list) {
	## We will remove single quote and/or space character from the param values.
	$record =~ s/'|\s//g;
	my @params = split( /:/, $record );
	$config{ $params[0] } = $params[1] if ( $params[1] );
	
	if ($params[0] eq "assembly_name") {
		$params[1] =~ s/<locus-name>/$locus_name/;
	}
	
	if ($params[0] eq "locus_tag") {
		$params[1] =~ s/<locus-name>/$locus_name/;
	}
	
	if ($params[0] eq "locus_id") {
		$params[1] =~ s/<locus-id>/$locus_id/;
	}
	
	print WRITE $params[0]."\n" unless ($params[1]);
	print WRITE $params[0].": ".$params[1]."\n" if ($params[1]);
}

close WRITE;

=head1 NAME

PERL script to generate config file based on the config template for "gmt hgmi hap" run 

Usage:
generate_config.pl --locus-name <locus-name>
	or
perl generate_config.pl --locus-name ENTFACTST --template ~/test-genome/Enterococcus_faecalis/EntFaecalis_template
