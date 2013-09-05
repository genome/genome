package Genome::Model::GenePrediction::Command::Pap::KEGGScan::KeggScanCommon;

# Subroutines used by Genome::Model::GenePrediction::Command::Pap::KEGGScan::RunKeggScan
# and Genome::Model::Tools::Predictor::Keggscan

use strict;
use warnings;

use Carp;

use base 'Exporter';
our @EXPORT_OK = qw(create_directory revise_subject_for_ECs build_EC_index create_report_line generate_blast_reports);

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

    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

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
    my ($self, $master_report, $top_report_path, $full_report_path, $EC_index, $KO_index) = @_;

    my $full_report_fh = IO::File->new($full_report_path, "w");
    my $top_report_fh = IO::File->new($top_report_path, "w");

    my $BP = $self->bpdelux_module_name->new(IO::File->new($master_report));

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

1;
