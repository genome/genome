package PAP::Command::KEGGScan::RunKeggScan;

# bdericks: This module was originally written by Todd Wylie and named
# KEGGscan_KO.v50.pl. The module can be found in version control at 
# svn+ssh://svn/srv/svn/gscpan/hgmi_annotation/proteinAnnot. I've made
# some small changes so this script can be included as a tool in the PAP
# namespace and included in all further snapshots.

use strict;
use warnings;

use PAP;

use Carp qw(confess);
use Cwd;
use File::Copy;
use File::Basename;
use App::DBI;
use File::Path qw(make_path);

# FIXME Really need to remove this dependency, but porting this script
# into the PAP namespace may not be straightforward... Used to generate
# the top and full reports after blast is run
use lib "/gsc/scripts/gsc/compbio/lib";
use BPdeluxe;

class PAP::Command::KEGGScan::RunKeggScan {
    is => 'PAP::Command',
    has => [
        species_name => {
            is => 'Text',
            doc => 'Name of input species',
        },
        query_fasta_path => {
            is => 'Path',
            doc => 'Path to query fasta file',
        },
        subject_fasta_path => {
            is => 'Path',
            doc => 'Path to subject fasta file',
        },
        output_directory => {
            is => 'Path',
            doc => 'All output is placed in this directory',
        },
    ],
    has_optional => [
        query_sequence_type => {
            is => 'Text',
            default => 'CONTIG',
            doc => 'Type of query sequence',
        },
        kegg_release_version => {
            is => 'Text',
            default => 'RELEASE-52',
            doc => 'Version of Kegg to use',
        },
        blast_lsf_queue => {
            is => 'Text',
            default => 'short',
            doc => 'Queue in which LSF jobs should be scheduled',
        },
        blast_lsf_resource => {
            is => 'Text',
            default => 'select[mem>1024] rusage[mem=1024]',
            doc => 'Resources needed to run LSF jobs',
        },
        blast_lsf_max_memory => {
            is => 'Number',
            default => '1024000',
            doc => 'Maximum memory limit of LSF jobs',
        },
        blast_lsf_job_limit => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of concurrent LSF jobs allowed',
        },
        fasta_chunk_size => {
            is => 'Number',
            default => 50,
            doc => 'Maximum number of sequence allowed in a fasta chunk',
        },
        db_format => {
            is => 'Text',
            calculate_from => 'kegg_release_version',
            calculate => q|
                $kegg_release_version =~ /\D*(\d+)\D*/;
                my $version_num = $1;
                return "new" if $version_num >= 41;
                return "old"; |,
            doc => 'The format of KO entries in genes db changed in version 41, ' .
                   'older versions require special handling',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->status_message("Creating kegg output directory at " . $self->output_directory);

    if (-d $self->output_directory) {
        $self->warning_message("Existing output directory found at " . $self->output_directory . ", removing!");
        my $rv = system("rm -rf " . $self->output_directory);
        confess "Could not remove " . $self->output_directory unless defined $rv and $rv == 0;
    }

    make_path($self->output_directory);
    chmod(0775, $self->output_directory);

    $self->status_message("Done creating output directory, now creating EC indices.");

    # We will be running BLAST (blastp) on the set of all sequences in the
    # query (species) FASTA file against ONLY the sequences in the KEGG
    # genes FASTA file that have associated EC numbers. Place the EC-only
    # FASTA file in the session output directory. We will be splitting up
    # the BLAST jobs across multiple CPUs via the blade center and
    # BladeBlastBatcher.pm.
    my $revised_subject = $self->revise_subject_for_ECs;
    my ($EC_index,$KO_index) = $self->build_EC_index($revised_subject);

    $self->status_message("Using BladeBlastBatcher to schedule blast jobs.");

    my $blast_batcher = PAP::Command::Blast::BladeBlastBatcher->create(
	    query_fasta_path => $self->query_fasta_path,
		subject_fasta_path => $revised_subject, 
        lsf_queue => $self->blast_lsf_queue,   
		lsf_job_limit => $self->blast_lsf_job_limit, 
		output_directory => $self->output_directory, 
		blast_name => "blastp",
		blast_params => "hitdist=40 wordmask=seg postsw B=10000 topcomboN=1",
        lsf_resource => $self->blast_lsf_resource,
        lsf_max_memory => $self->blast_lsf_max_memory,
        fasta_chunk_size => $self->fasta_chunk_size,
	);
    my $batcher_rv = $blast_batcher->execute;
    confess "Trouble executing blast batcher!" unless defined $batcher_rv and $batcher_rv == 1;
    my @reports = @{$blast_batcher->reports};

    $self->status_message("BladeBlastBatcher completed, aggregating reports.");
    $self->status_message("There are " . scalar @reports . " that need to be aggregated!");

    my $master_report_path = $self->output_directory . "/MASTER_BLAST.report";
    my $aggregate = IO::File->new($master_report_path, "a");
    REPORT: for my $report (@reports) {
        unless (-e $report) {
            $self->warning_message("Report at $report doesn't exist!");
            next REPORT;
        }

        unless (-s $report) {
            $self->warning_message("Report at $report has no size!");
            unlink $report;
            next REPORT;
        }

        my $fh = new IO::File->new($report, "r");
        while (my $line = $fh->getline) { 
            $aggregate->print($line);
        }
        $fh->close;
        unlink($report);   # Saves a little space.
    }
    $aggregate->close;

    unless (-s $master_report_path) {
        confess "Master blast report at $master_report_path has no size!";
    }

    $self->status_message("All reports aggregated, now parsing.");

    # Walk through all of the blast report files and gather the needed info. We
    # will be using BPdeluxe.pm to parse the BLAST reports. Our goal is to get
    # one best (top) match per query with regard to E-value & bit score. If
    # there are multiple/identical E-value matches under a single query, then
    # best bitscore will determine the top hit. If both E-value and bitscore
    # are identical we will simply use the first match in the list by default.
    my $full_report_path    = $self->output_directory . "/" . "REPORT-full.ks";
    my $top_hit_report_path = $self->output_directory . "/" . "REPORT-top.ks";
    $self->generate_blast_reports($master_report_path, $top_hit_report_path, $full_report_path, $EC_index, $KO_index);
    $self->status_message("Report parsing completed, starting clean up.");

    # Move files to subdirectories to reduce clutter
    my $ancillary_dir = $self->create_directory($self->output_directory . "/ancillary/");
    my $log_dir = $self->create_directory($self->output_directory . "/logs/");

    opendir(DIR, $self->output_directory) or confess "Couldn't open " . $self->output_directory;
    while (my $file = readdir(DIR)) {
        my $destination;
        if ($file =~ /^CHUNK/) {
            unlink ($self->output_directory . "/$file");
        }
        elsif (($file =~ /OU$/) || ($file =~ /EC_only/)) {
            $destination = $ancillary_dir
        }
        elsif (($file =~ /.out$/) || ($file =~ /.err$/) || ($file =~ /.log$/)) {
            $destination = $log_dir;
        }
        next unless defined $destination;

        my $full_path = $self->output_directory . "/$file";
        my $mv_rv = system("mv $full_path $destination");
        $self->warning_message("Could not move $full_path to $destination!") unless $mv_rv == 0;
    }

    $self->status_message("Keggscan finished, reports can be found at:\n$full_report_path\n$top_hit_report_path");
    return 1;
}

# ===========================================================================
# S U B R O U T I N E S
# ===========================================================================

# Simply creates a path and changes the permissions
sub create_directory {
    my ($self, $dir) = @_;
    make_path($dir);
    chmod(0755, $dir);
    return $dir;
}

# REVISE SUBJECT FOR ECS
# This routine will take a KEGG genes FASTA file and revise it,
# creating a version with only entries that have associated Enzyme
# Commission (EC) numbers. The file will be written to the session
# output directory. We will return the path to the revised file.
sub revise_subject_for_ECs {
    my $self = shift;
    my $subject_fasta = $self->subject_fasta_path;
    my $revised_subject = $self->output_directory . "/genes.EC_only";

	$ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    local $/ = "\n>";
    my $fh_in  = IO::File->new($self->subject_fasta_path, "r"); 
    my $fh_out = IO::File->new($revised_subject, "w");

    while (my $line = $fh_in->getline) {
        chomp $line;

        if ($line =~ /EC\:[0-9]+\.[0-9-]+\.[0-9-]+\.[0-9-]+/) {
            if ($line =~ /^>/) {
                $fh_out->print($line);
            } 
            else {
                $fh_out->print("\n>$line");
            }
        }
    }

    $fh_in->close;
    $fh_out->close;

    $self->status_message("Running xdformat on $revised_subject.");
    my $rv = system("xdformat", "-p", $revised_subject);
    $self->status_message("Return value from system command: $rv");
	confess "Failed during xdformat: $!\n" unless defined $rv and $rv == 0;

	$ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 0;

    return $revised_subject;
}

# BUILD EC INDEX
# This routine will take the revised subject file and create an index of
# subject name to EC number association.
sub build_EC_index {
    my $self = shift;
    my $revised_subject = shift;

    my $fh = IO::File->new($revised_subject, "r");
    my $errors = $self->output_directory . "/error.log";
    my $error_log = IO::File->new($errors, "a");

    my $inter = 0;
    my %ECs;
    my %KOs;

    while (my $line = $fh->getline) {
        $inter++;
        # NOTE: Looks like KEGG genes files always have subject_name in 1st
        # col & EC in # last. If this ever changes, we'll have to update this routine....
        # Genesdb started to append KO numbers after EC numbers and this routine was modified to capture the KO numbers
        if ( $line =~ /^>/ ) {
            my $gi;
            my $meta;
            my $ECs;
            my $KOs;
            if ($self->db_format eq "old") { # This branch handles pre-release-41 KO values
                if ($line =~ /\[KO\:.+\]/) {
                    ($gi, $meta, $ECs, $KOs) = $line =~ /^\>\s*([a-zA-Z0-9-_.]+\:[a-zA-Z0-9-_.]+)\s*(.+.)\s*(\[EC\:.+\]\s*).+(\[KO\:.+\]\s*)/;
                    my @ko_pull = split(/\]/, $KOs);
                    $KOs = shift(@ko_pull);
                    $KOs =~ s/\[KO\://;
                    my @kids = split(/\s+/, $KOs);
                    my $kids = join(":", @kids);

                    if ( ($gi) && ($meta) && ($kids) ) {
                        $KOs{$gi}   = [ $meta, $kids ];	   
                    }
                    else {
                        $error_log->print("$inter FALL-THROUGH EVENT: $line\n\n");
                    }
                }
                else {
                    ($gi, $meta, $ECs) = $line =~ /^\>\s*([a-zA-Z0-9-_.]+\:[a-zA-Z0-9-_.]+)\s*(.+.)\s*(\[EC\:.+\]\s*)/;	
                }
            } 
            elsif ($self->db_format eq "new") {# This branch handles release-41 and above KO values
                if ($line =~ /; K/){
                    ($gi, $meta, $ECs, $KOs) = $line =~ /^\>\s*([a-zA-Z0-9-_.]+\:[a-zA-Z0-9-_.]+)\s*(.+.)\s*(\[EC:.+\]|\(EC:.+\)\s*)\;\s+(K.+)$/;

                    if (defined $KOs) {
                        my @ko_pull = split(/\;/, $KOs);
                        my @kids;

                        for my $ko_entry (@ko_pull) {
                            $ko_entry =~ s/^\s+//;
                            my @ko_string = split(/\s+/,$ko_entry);
                            push(@kids,$ko_string[0]);
                        }

                        my $kids = join(":", @kids);
                        if ( ($gi) && ($meta) && ($kids) ) {
                            $KOs{$gi}   = [ $meta, $kids ];	   
                        }
                        else {
                            $error_log->print("$inter FALL-THROUGH EVENT: $line\n\n");
                        }
                    }
                    else {
                        ( $gi, $meta, $ECs ) = $line =~/^\>\s*([a-zA-Z0-9-_.]+\:[a-zA-Z0-9-_.]+)\s*(.+.)\s*(\[EC:.+\]|\(EC:.+\)\s*)/;
                    }
                }
                else {
                    ( $gi, $meta, $ECs ) = $line =~/^\>\s*([a-zA-Z0-9-_.]+\:[a-zA-Z0-9-_.]+)\s*(.+.)\s*(\[EC:.+\]|\(EC:.+\)\s*)/;	
                }
            } 
            else {
                die "\n\nUnable to determine Kegg release format of genes db.\n\n";
            }

            ## We will substitute () with []
            $ECs =~ s/\(/\[/;
            $ECs =~ s/\)/\]/;

            my @ec_pull = split(/\]/, $ECs);
            $ECs        = shift(@ec_pull);
            $ECs        =~ s/\[EC\://;
            my @ids     = split(/\s+/, $ECs);
            my $ids     = join(":", @ids);

            if ( ($gi) && ($meta) && ($ids) ) {
                $ECs{$gi}   = [ $meta, $ids ];	    
            }
            else {
                $error_log->print("$inter FALL-THROUGH EVENT: $line\n\n");
            }
        }
    }

    $error_log->close;
    return(\%ECs,\%KOs);
}

# REPORT TO FILE
# This routine will take a given reference to a hash created by calling
# BlastTopHitLogic.pm and then print the appropriate output as needed for
# KEGGscan. Output is tab-delimited, the columns corresponding to fields in
# the KEGG MySQL database.
sub create_report_line {
    my ($self, $blast_parse, $EC_index, $KO_index)  = @_;
    my $errors = $self->output_directory . "/error.log";
    my $error_log = IO::File->new($errors, "a");
    my %check;

    # Species:
    if ($self->species_name) { 
        $check{1} = [ "species", "yes", $self->species_name];
    } 
    else { 
        $check{1} = [ "species", "no" ];
    }
    # ECs:
    if ( $$EC_index{$blast_parse->[1]}[1] ) {
        $check{2} = [ "ECs", "yes", $$EC_index{$blast_parse->[1]}[1] ];
    }
    else {
        $check{2} = [ "ECs", "no" ];
    }
    # KOs:
    if ( $$KO_index{$blast_parse->[1]}[1] ) {
        $check{9} = [ "KOs", "yes", $$KO_index{$blast_parse->[1]}[1] ];
    }
    else {
        $check{9} = [ "KOs", "yes", "none" ];
    }
    # QUERY GI:
    if ( $blast_parse->[0] ) {
        $check{3} = [ "query GI", "yes", $blast_parse->[0] ];
    }
    else {
        $check{3} = [ "query GI", "no" ];
    }
    # SUBJECT GI:
    if ( $blast_parse->[1] ) {
        $check{4} = [ "subject GI", "yes", $blast_parse->[1] ];
    }
    else {
        $check{4} = [ "subjct GI", "no" ];
    }
    # P VALUE:
    if ( $blast_parse->[2]->{p_value} ) {
        $check{5} = [ "p_value", "yes", $blast_parse->[2]->{p_value} ];
    }
    else {
        $check{5} = [ "p_value", "no" ];
    }
    # BIT SCORE:
    if ( $blast_parse->[2]->{bit_score} ) {
        $check{6} = [ "bit score", "yes", $blast_parse->[2]->{bit_score} ];
    }
    else {
        $check{6} = [ "bit score", "yes" ];
    }
    # SUBJECT META:
    if ( ${$$EC_index{$blast_parse->[1]}}[0] ) {
        $check{7} = [ "subject meta", "yes", ${$$EC_index{$blast_parse->[1]}}[0] ];
    }
    else {
        $check{7} = [ "subject meta", "no" ];
    }
    # QUERY SEQUENCE TYPE:
    if ($self->query_sequence_type) {
        $check{8} = [ "query sequence type", "yes", $self->query_sequence_type ];
    }
    else {
        $check{8} = [ "query sequence type", "no" ];
    }

    # Failure check & output:
    my @failures;
    for my $field (keys %check) {
        if ($check{$field}[1] eq "no") {
            push (@failures, $check{$field}[0]);
        }
    }

    my $report_line;
    if (@failures > 0) {
        for my $field (keys %check) {
            if ($check{$field}[1] eq "yes") {
                $error_log->print("VIABLE: ", $check{$field}[2], "\n");
            }
        }
        $error_log->print("ERROR(S): ", join("--", @failures), "\n||\n");
    }
    else {
        my @report_line;
        foreach my $field (sort keys %check) {
            push(@report_line, $check{$field}[2]);
        }
        $report_line = join("\t", @report_line);
    }

    $error_log->close;
    return $report_line;
}

# This used to be the BlastTopHitLogic module, but has been incorporated
# into this module to remove the need for a giant hash that was leading
# to some memory issues
sub generate_blast_reports {
    $DB::single = 1;
    my ($self, $master_report, $top_report_path, $full_report_path, $EC_index, $KO_index) = @_;

    my $full_report_fh = IO::File->new($full_report_path, "w");
    my $top_report_fh = IO::File->new($top_report_path, "w");

    my $BP = new BPdeluxe::Multi(IO::File->new($master_report));

    while (my $multi = $BP->nextReport) {
        my $query = $multi->query;
        my %grouped;
        while(my $sbjct = $multi->nextSbjct) {	
            my $group_hashref = $sbjct->group_list("$query","return_hashref");
            $grouped{$group_hashref->{subject_gi}} = $group_hashref;
        }

        my %final_vals;
        my $top_p_value;
        foreach my $val (sort {$grouped{$a}->{p_value} <=> $grouped{$b}->{p_value}} keys %grouped) {
            my @full_report_values = ($query, $val, $grouped{$val});
            my $full_report_line = $self->create_report_line(\@full_report_values, $EC_index, $KO_index);
            $full_report_fh->print($full_report_line . "\n") if defined $full_report_line and $full_report_line ne '';;

            # Keep track of top value for top report
            $top_p_value = $grouped{$val}->{p_value} if not defined $top_p_value;
            $final_vals{$val} = $grouped{$val} if ($grouped{$val}->{p_value} eq $top_p_value);
        }

        # Write top p value entry to top report
        my @p_vals = sort {$final_vals{$b}->{bit_score} <=> $final_vals{$a}->{bit_score}} keys %final_vals;
        my $top_p = shift @p_vals;
        next unless defined $top_p;
        my @top_report_values = ($query, $top_p, $final_vals{$top_p});
        my $top_report_line = $self->create_report_line(\@top_report_values, $EC_index, $KO_index);
        $top_report_fh->print($top_report_line . "\n") if defined $top_report_line and $top_report_line ne '';
    }

    $full_report_fh->close;
    $top_report_fh->close;
    return 1;
}

# *|**********************************************************************|*
# *| POD: Section below is reserved for documentation.                    |*
# *| Last update: Thu Oct 20 10:08:34 CDT 2005                            |*
# *|**********************************************************************|*

=head1 NAME

RunKeggScan.pm

=head1 USAGE

 Your query sequence file should be of nucleotide FASTA format. The subject database should be a FASTA of protein sequences as retrieved from the KEGG FTP site (ftp://ftp.genome.ad.jp/pub/kegg/). As of this writing, the genes file can be found by downloading and unpacking sequences.tar.gz from KEGG; the file needed is simply called "genes".

 NOTE: The 'KeggRelease' line entry must ONLY contain a single numeric portion.  This code uses this number to determine the expected formatting of the genes db fasta file with regards to the KO numbers.

=head1 DESCRIPTION

KEGGscan was written to compare nucleotide sequence (full length cDNAs, genes, clustered ESTs) to the KEGG gene database of /currated protein sequences. Sequences in the KEGG database have known, annotated Enzyme Commission (EC) ids associated with them. By aligning query sequence against annotated sequence we can assign putative function by EC number association. Obviously, the better the resolution of the query sequence, the more believable our results will be.  The BLAST used to make our alignments uses the topcomboN=1 argument to return only the "best" result. A routine takes the KEGG genes FASTA file and revises it, creating a version with only entries that have associated Enzyme Commission (EC) numbers.

The following logic filters the match list. Top hits are determined on a one-to-one basis (1 group hit per query). The rules for hits are:

=over 5

=item 1

Lowest p_value determines top match.                                   

=item 2

If other HSPs tie the lowest p_value, highest bit score prevails.      

=item 3

If p_value, bit scores are identical, then choosing any 1 of the hits is equally representative.              

=back

=head1 AUTHOR

 Todd Wylie
 CompBio
 Genome Sequencing Center
 Washington University 
 School of Medicine
 4444 Forest Park Boulevard
 St. Louis, MO 63108

CONTACT: twylie@.wustl.edu

 Modifications made by Brian Derickson
 Analysis Pipeline Group
 The Genome Center at Washington University in St. Louis

CONTACT: bdericks@genome.wustl.edu

=head1 LIMITATION/BUGS

Requires:

=over 5

=item 1

BladeBlastBatcher.pm

=item 2

BlastTopHitLogic.pm  

=back

NOTE: Looks like KEGG genes files always have subject_name in 1st col & EC in last. If this ever changes, we'll have to update this routine....

=head1 COPYRIGHT

Copyright (C) 2005 by Todd Wylie and Washington University School of Medicine 
Genome Sequencing Center.

=head1 NOTE:

This software was written using the latest version of GNU Emacs, the extensible, 
real-time text editor. Please see http://www.gnu.org/software/emacs/ for more 
information and download sources.

=cut
