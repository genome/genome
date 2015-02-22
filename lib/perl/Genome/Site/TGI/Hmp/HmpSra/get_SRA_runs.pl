#!/usr/bin/env genome-perl

# Use Aspera command-line client to download a specified list of SRA runs (from a list of SRR accession numbers)

use strict;
use FileHandle;

## globals
my $USAGE = "Usage: $0 ascp_path index_file SRR_accession_list";

## input
my $ascp_path = shift || die $USAGE;
my $ifile = shift || die $USAGE;
my $srr_accs = shift || die $USAGE;
# TODO - allow specification of download dir
# TODO - clean up script to same format as the rest

## main program
my $key_file = $ascp_path;
print "ASCP path: $ascp_path\n";
####$key_file =~ s/\/bin\/ascp/\/etc\/asperaweb_id_dsa.putty/;
#modified for WUGC use ... jmartin 100806
$key_file = "/opt/aspera/etc/asperaweb_id_dsa.putty";

# read index file
my $index = {};
my $ifh = FileHandle->new();
$ifh->open($ifile, 'r') || die "unable to read from $ifile";
while (my $line = <$ifh>) {
    chomp($line);
    my($srr_acc, $ftp_url, $date) = split(/\t/, $line);
    if (defined($index->{$srr_acc})) {
        print "WARNING - ignoring duplicate entry for '$srr_acc'\n";
        next;
    }
    $index->{$srr_acc} = { 'accession' => $srr_acc, 'ftp_url' => $ftp_url, 'date' => $date };
}
$ifh->close();
print "read " . scalar(keys %$index) . " entries from $ifile\n";

# read SRR accessions
my $srr_list = [];
my $afh = FileHandle->new();
$afh->open($srr_accs, 'r') || die "unable to read from $srr_accs";
my $lnum = 0;
my $na = 0;
while (my $line = <$afh>) {
    chomp($line);
    ++$lnum;
    if ($line =~ /^(SRR\d+)/) {
        ++$na;
        # TODO - add option to halt instead of printing a warning when an accession is missing
        my $acc = $1;
        my $srr = $index->{$acc};
        if (!defined($srr)) {
            print "no URL found for $acc in $ifile\n";
        } else {
            push(@$srr_list, $srr);

            # convert FTP URL to FASP location like anonftp@ftp-private.ncbi.nlm.nih.gov:22/sra0/SRR/000004/SRR004103/
            my $ftp_url = $srr->{'ftp_url'};
            my($host, $port, $path) = ($ftp_url =~ /ftp:\/\/([^:]+):(\d+)(\/.*)/);
            die "couldn't parse FTP URL $ftp_url" if (!defined($path));
            # HACK
            $host =~ s/ftp\.ncbi/ftp-private\.ncbi/;
            $srr->{'fasp_location'} = "anonftp\@${host}:${path}";
        }
    } else {
        die "unrecognized SRR accession at line $lnum of $srr_accs";
    }
}
$afh->close();
print "found URLs for " . scalar(@$srr_list) . "/$na SRR accessions\n";
my $nret = 0;
my $nfailed = 0;

# retrieve SRR entries with ascp
foreach my $srr (@$srr_list) {
    my $fasp_locn = $srr->{'fasp_location'};
    my $accn = $srr->{'accession'};
    my $date = `date`;
    chomp($date);

    # -Q = fair transfer policy
    # -T = no encryption
    # -i = pointer to key file

#modified for WUGC use ... jmartin 100812
    ####my $cmd = "$ascp_path -QT -i $key_file '$fasp_locn' ./$accn &>$accn.log";
    my $cmd = "$ascp_path -T -l 100M -i $key_file '$fasp_locn' ./$accn";
    #my $cmd = "$ascp_path -T -l 20M -i $key_file '$fasp_locn' ./$accn";
    #my $cmd = "$ascp_path -T -l 50M -i $key_file '$fasp_locn' ./$accn";

    #print STDERR "cmd=$cmd\n";

    # TODO - make this configurable
    if ((-e $accn) && (-d $accn)) {
        print "$date: $fasp_locn (already retrieved, skipping)\n";
        next;
    }

    print "$date: $fasp_locn\n";


    ####system($cmd);
    print "RUNNING: $cmd\n";
    print "OUTPUT LOCATION: " . `pwd`;
#modified for WUGC use ... jmartin 100806
   system("
$cmd <<EOF
y
EOF
");



    my $failed = 0;

    if ($? == -1) {
        print "ascp command failed to execute: $!\n";
        $failed = 1;
    }
    elsif ($? & 127) {
        printf "ascp command died with signal %d, %s coredump\n", ($? & 127),  ($? & 128) ? 'with' : 'without';
        $failed = 1;
    }
    else {
        my $exitval = $? >> 8;
        if ($exitval != 0) {
            printf "ascp command exited with value $exitval\n";
            $failed = 1;
        }
    }

    if ($failed) {
        print "retrieval of $accn failed, removing target dir\n";
        system("/bin/rm -rf $accn/");
        ++$nfailed;
    } else {
        ++$nret;
    }
}
print "transferred $nret SRA runs with ascp, $nfailed failure(s) detected\n";
exit(0);
