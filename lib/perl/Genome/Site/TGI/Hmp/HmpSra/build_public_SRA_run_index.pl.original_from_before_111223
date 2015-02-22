#!/usr/bin/env perl

=head1  NAME 

build_public_SRA_run_index.pl - Builds an index of all the SRA runs (SRR accessions) available on the public SRA site.

=head1 SYNOPSIS

build_public_SRA_run_index.pl
        [--output_dir=.
         --ftp_subdirs=sra0,sra1,sra2
         --reuse_files
    ]

=head1 OPTIONS

<--output_dir,-o>
    Where to put the intermediate and final output.

<--ftp_subdirs,-s>
    Comma-delimited list of top-level directories to examine on the FTP site.  Default = 'sra0,sra1,sra2'

<--reuse_files,-r>
    Reuse downloaded files (no questions asked) if they already exist on disk.  (Default = no)

=head1 DESCRIPTION

This script uses wget to build an index of all the SRR accessions available for fasp/Aspera download from
the NCBI's FTP site at ftp-private.ncbi.nlm.nih.gov.  It assumes that the directory structure (above the
SRR runs themselves) is two levels deep: it first lists every subdirectory of each directory named in 
--ftp_subdirs, and then it recursively lists each of those subdirectories, expecting to find [directories
corresponding to] SRR accessions therein.

=head1 INPUT

A comma-delimited list of the top-level directories to traverse on the NCBI FTP site.

=head1 OUTPUT

A tab-delimited list mapping each available SRA SRR run accession to an Aspera fasp URL.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use File::Spec;
use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

## globals
#my $SRA_FTP_HOST = 'ftp-private.ncbi.nlm.nih.gov';
my $SRA_FTP_HOST = 'ftp.ncbi.nlm.nih.gov';
my $DEFAULT_OUTPUT_DIR = '.';
# as of 9/14/10 the 'sra0' ftp directory was not visible on NCBI FTP site and thus this script does not run...but hopefully this is just a temporary bug on their end ... jmartin 100914
my $DEFAULT_FTP_SUBDIRS = 'sra0,sra1,sra2';

## input
my $options = {};
&GetOptions($options,
            "output_dir|o=s",
            "ftp_subdirs|s=s",
            "reuse_files|r!",
            "help|h",
            ) || pod2usage();

## display documentation
if( $options->{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters($options);

## main program
my($output_dir, $ftp_subdirs, $reuse_files) = map { $options->{$_} } ('output_dir', 'ftp_subdirs', 'reuse_files');
die "output directory ($output_dir) does not exist or is not writable" if ((!-e $output_dir) || (!-w $output_dir));

my @ftp_sd = split(/\s*,\s*/, $ftp_subdirs);
die "no --ftp_subdirs specified to search" if (scalar(@ftp_sd) == 0);

# build index 1 level deep
my $index1 = [];
foreach my $sd (@ftp_sd) {
    my $dir = "/" . $sd . "/SRR/";
    my $list = &get_ftp_dir_listing($SRA_FTP_HOST, $dir, $output_dir, $reuse_files);
    map { $_->{'d1'} = $sd; } @$list;
    push(@$index1, @$list);
    my $ll = scalar(@$list);
    print STDERR "read $ll level 1 directories from $SRA_FTP_HOST${dir}\n";
}
my $n1 = scalar(@$index1);
print STDERR "read $n1 level 1 directories from $SRA_FTP_HOST\n";

# build index 2 levels deep
my $index2 = [];
foreach my $sd (@$index1) {
    my($date, $url, $name, $d1) = map {$sd->{$_}} ('date', 'url', 'name', 'd1');
    my $dir = "/" . $d1 . "/SRR/" . $name . "/";
    my $list = &get_ftp_dir_listing($SRA_FTP_HOST, $dir, $output_dir, $reuse_files);
    map { $_->{'d2'} = $name; } @$list;
    push(@$index2, @$list);
    my $ll = scalar(@$list);
    print STDERR "read $ll level 2 directories from $SRA_FTP_HOST${dir}\n";
}
my $n2 = scalar(@$index2);
print STDERR "read $n2 level 2 directories from $SRA_FTP_HOST\n";

foreach my $i2 (@$index2) {
    my($date, $url, $name, $d1, $d2) = map {$i2->{$_}} ('date', 'url', 'name', 'd1', 'd2');
    print join("\t", $name, $url, $date) . "\n";
}

exit(0);

## subroutines
sub check_parameters {
    my $options = shift;
    
    ## defaults
    $options->{'output_dir'} = $DEFAULT_OUTPUT_DIR if (!defined($options->{'output_dir'}));
    $options->{'ftp_subdirs'} = $DEFAULT_FTP_SUBDIRS if (!defined($options->{'ftp_subdirs'}));
}

# retrieve the listing of a directory on the remote FTP site
sub get_ftp_dir_listing {
    my($host, $dir, $output_dir, $reuse_files) = @_;
    my $url = "ftp://" . $host . $dir;

    my $ld = $dir;
    $ld =~ s/^\/+//;
    $ld =~ s/\/+$//;
    $ld =~ s/\//\-/g;

    my $log_file = File::Spec->catfile($output_dir, $ld . ".log");
    my $html_file = File::Spec->catfile($output_dir, $ld . ".html");
    my $list = [];

    if ($reuse_files && (-e $html_file) && (-e $log_file)) {
        print STDERR "$html_file already downloaded and --reuse_files specified\n";
    } else {
        my $cmd = "wget '$url' -o $log_file -O $html_file";

        print STDERR "cmd=$cmd\n";
        system($cmd);

        # make sure that the exit value is propagated
        my $exitval = undef;
        
        if ($? == -1) {
            print STDERR "$0 failed executing '$cmd': $!\n";
            $exitval = 1;
        }
        elsif ($? & 127) {
            printf STDERR "$0 failed executing '$cmd' - child died with signal %d, %s coredump\n", ($? & 127),  ($? & 128) ? 'with' : 'without';
            $exitval = 1;
        }
        else {
            $exitval = $? >> 8;
            if ($exitval != 0) {
                print STDERR "$0 failed executing '$cmd' - child exited with value $exitval\n";
            }
        }
        exit($exitval) if ($exitval != 0);
    }        

    # parse the results
    my $fh = FileHandle->new();
    $fh->open($html_file, 'r') || die "unable to read from $html_file";
    my $lnum = 0;
    while (my $line = <$fh>) {
        chomp($line);
        ++$lnum;

        if ($line =~ /\s+(\d{4} \S{3} \d+ [\d:\s]+)\s+Directory\s+\<a href=\"([^\"]+)\"\>([^\/]+)\//) {
            my($date, $url, $name) = ($1, $2, $3);
            $date =~ s/\s+$//;
            push(@$list, { 'date' => $date, 'url' => $url, 'name' => $name });
            print STDERR "date=$date url=$url name=$name\n";
        }
    }

    $fh->close();
    return $list;
}
