package Genome::Model::Tools::Assembly::ReScaffoldMsiAce;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Assembly::ReScaffoldMsiAce {
    is => 'Command',
    has => [
	acefile => {
	    is => 'Text',
	    doc => 'Assembly ace file',
	},
	scaffold_file => {
	    is => 'Text',
	    doc => 'Assembly scaffold file',
	    is_optional => 1,
	},
	auto_report => {
	    is => 'Boolean',
	    doc => 'Run consed autoreport to get scaffold info',
	    is_optional => 1,
	},
	assembly_directory => {
	    is => 'Text',
	    doc => 'Main assembly directory',
	},
        min_contig_length => {
            is => 'Integer',
            doc => 'Minimum contig length to export to new ace file',
            default => 0,
            is_optional => 1,
        },
        min_scaffold_length => {
            is => 'Integer',
            doc => 'Minimum scaffold length to export to new ace file',
            default => 0,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to re-scaffold msi assemblies';
}

sub help_synopsis {
    return <<EOS
gmt assembly re-scaffold-msi-ace --acefile /gscmnt/111/assembly/edit_dir/newbler.ace --assembly-directory /gscmnt/111/assembly
gmt assembly re-scaffold-msi-ace --acefile /gscmnt/111/assembly/edit_dir/newbler.ace --scaffold-file /gscmnt/111/assembly/edit_dir/scaffolds --assembly-directory /gscmnt/111/assembly
gmt assembly re-scaffold-msi-ace --acefile /gscmnt/111/assembly/edit_dir/newbler.ace --assembly-directory /gscmnt/111/assembly --auto-report
EOS
}

sub help_detail {
    return <<EOS
Tool to re-scaffold manually edited ace files.  It will read in
a text file of scaffolding info and re-name contigs to form
pcap style scaffolds, eg, line 8.9-7,15-1.2-1.3 will create the
following scaffolds:

Contig0.1-Contig0.2           (from Contigs8.9 and 7)
Contig1.1-Contig1.2-Contig1.3 (from Contigs15, 1.2 and 1.3)

Tool will also create msi.gap.txt file which is later used
to determine gap sizes between scaffold contigs.  Each gap
will be assigned default, unknown value of 100 bp.

This tool will work with any acefiles from any assembler.
EOS
}

sub execute {
    my $self = shift;

    #make sure inputs are correct
    if ( not $self->_validate_inputs ) {
        $self->error_message( "Failed to validate tool inputs" );
        return;
    }

    #get contig lengths from input ace file
    $self->debug_message( "Determining old scaffolds" );
    my $contig_lengths;
    if ( not $contig_lengths = $self->_get_unpadded_contig_lengths( $self->acefile ) ) {
        $self->error_message( "Failed to get contig lengths from input ace file" );
        return;
    }
 
    #determine new scaffolds
    $self->debug_message( "Determining new scaffolds" );
    my $new_scaffolds;
    if ( not $new_scaffolds = $self->_determine_new_scaffolds( $contig_lengths ) ) {
        $self->error_message( "Failed to determine new scaffolds" );
        return;
    }
    
    #make partial ace of new scaffolds
    $self->debug_message( "Writing new scaffolds" );
    my ($partial_ace, $contigs, $reads ) = $self->_write_new_scaffolds( $self->acefile, $new_scaffolds );

    #update DS line .. needed for older newbler ace files .. and update read and contig counts
    $self->debug_message( "Updating DS line and writing new ace file: ace.msi" );
    my $final_ace = $self->_update_ds_line_write_wa_tags( $partial_ace, $contigs, $reads );

    #TODO sort this ace file numerically by contig number

    $self->debug_message( "Done" );
    return 1;
}

sub _determine_new_scaffolds {
    my ( $self, $contig_lengths ) = @_;
    
    my $scaffolds_file;
    if ( not $scaffolds_file = $self->_set_new_scaffolds_file ) {
        $self->error_message( "Failed to set scaffolds file from either supplied or autoreport" );
        return;
    }

    #create new scaffolds
    my $new_scaffolds;
    if ( $scaffolds_file ) {
	my $scaffolds;
        unless ( $scaffolds = $self->_parse_scaffolds_file( $scaffolds_file ) ) {
            $self->error_message( "Failed to parse scaffolds file" );
            return;
        }
	my $valid_scaffolds;
        unless( $valid_scaffolds = $self->_check_for_contigs_to_complement($scaffolds) ) {
            $self->error_message( "Failed to check contigs to complement" );
            return;
        }
	unless( $new_scaffolds = $self->_create_new_scaffolds($contig_lengths, $valid_scaffolds) ) {
            $self->error_message( "Failed to create new scaffolds" );
            return;
        }
    }
    else {
        unless ( $new_scaffolds = $self->_create_new_scaffolds($contig_lengths) ) {
            $self->error_message( "FAILED to create new scaffolds" );
        }
    }

    return $new_scaffolds;
}


sub _validate_inputs {
    my $self = shift;

    unless ( $self->acefile and -s $self->acefile ) {
        $self->error_message( "Can't find ace file or file is zero size or file was not supplied" );
        return;
    }
    if ($self->scaffold_file and $self->auto_report) {
	$self->error_message("You can't select to run auto report and supply scaffold file");
	return;
    }
    return 1;
}

sub _set_new_scaffolds_file {
    my $self = shift;

    my $report_file;
    if ($self->auto_report) {
	unless ($report_file = $self->_run_auto_report()) {
	    $self->error_message("Failed to get report file by running auto report");
	    return;
	}
    }

    if ($self->scaffold_file) {
	unless (-s $self->scaffold_file) {
	    $self->error_message("Can't find scaffold file: ".$self->scaffold_file);
	    return;
	}
	$report_file = $self->scaffold_file;
    }

    return $report_file;
}

sub _run_auto_report {
    my $self = shift;
    $self->debug_message("Running consed auto report");
    my $acefile = $self->acefile;
    if (system("consed -ace $acefile -autoreport")) {
	$self->error_message("Failed to run consed auto report on ace file: $acefile");
	return;
    }
    my $dir = $self->assembly_directory.'/edit_dir';
    my @out_files = `ls -t $dir/*[0-9]\.out`; #grap autoreport output files
    unless (@out_files) {
	$self->error_message("Failed to find any auto report output files with file format ####.out");
	return;
    }
    chomp @out_files;

    return shift @out_files;
}

sub _parse_scaffolds_file {
    my ($self, $file) = @_;
    my @scaffolds;
    my $fh = Genome::Sys->open_file_for_reading($file);
    foreach my $line ($fh->getlines) {
	next unless ($line =~ /^\d+/ or $line =~ /^E-\d+/);
	chomp $line;
	if ($line =~ /,/) {
	    my @tmp = split (',', $line);
	    foreach (@tmp) {
                $_ =~ s/^\s+//;
		push @scaffolds, $_;
	    }
	} else {
	    push @scaffolds, $line;
	}
    }

    $fh->close;

    return \@scaffolds;
}

sub _check_for_contigs_to_complement {
    my ($self, $scaffolds) = @_;
    my $contigs_to_complement;
    foreach (@$scaffolds) {
	my @tmp = split ('-', $_);
	foreach (@tmp) {
	    $contigs_to_complement .= $_.', ' if $_ =~ /c/;
	}
    }
    if ($contigs_to_complement) {
	if ($self->auto_report) {
	    $self->debug_message("\n\nConsed autoreport suggests that the following contigs must be complemented:\n".
				  "\t$contigs_to_complement .. please complement these contigs in the ace file and run the program again\nExiting");
	}
	if ($self->scaffold_file) {
	    $self->debug_message("\nPlease complement the following contigs in the ace file: $contigs_to_complement\n".
				  "Then remove c from contig numbers then run the program again to reflect correct compelementation in post assembly files");
	}
	return;
    }
    return $scaffolds;
}

#need to get unpadded contig lengths w/o using
#ace object .. ace files can get too big to load
sub _get_unpadded_contig_lengths {
    my ( $self, $ace ) = @_;
    my %contig_lengths;
    my $fh = Genome::Sys->open_file_for_reading( $ace );
    my $contig_name;
    my $is_sequence = 0;
    while ( my $line = $fh->getline ) {
        if ( $line =~ /^CO\s+/ ) {
            ( $contig_name ) = $line =~ /^CO\s+(\S+)/;
            $contig_name =~ s/contig//i;
            $is_sequence = 1;
            next;
        }
        if ( $line =~ /^BQ\s+/ ) {
            $is_sequence = 0;
            next;
        }
        if ( $is_sequence == 1 ) {
            next if $line =~ /^\s+$/;
            chomp $line;
            #$line =~ s/[nx*]//ig;
            $line =~ s/n//ig;
            $line =~ s/x//ig;
            $line =~ s/\*//g;
            $contig_lengths{$contig_name} += length $line;
        }
    }
    $fh->close;

    return \%contig_lengths;
}

sub _create_new_scaffolds {
    my ($self, $old_contigs, $scaffolds) = @_;

    my %scaffold_lengths;
    if( $scaffolds ) {
        foreach my $scaf ( @$scaffolds ) {
            $scaf =~ s/\s+//;
            my @tmp = split (/-/, $scaf);
            foreach my $scaf_ctg (@tmp) {
                next if $scaf_ctg eq 'E'; #eg E-12.1-E
                my $scaf_num;
                if( $scaf_ctg =~ /^\d+\.\d+$/ ) {
                    ($scaf_num) = $scaf_ctg =~ /^(\d+)\./;
                }
                elsif ( $scaf_ctg =~ /^\d+$/ ) {
                    $scaf_num = $scaf_ctg;
                }
                else {
                    die "Could not get scaffold number from contig: $scaf_ctg\n";
                }
                die "Contig $scaf_ctg exists in scaffolds file but not in ace file\n" if
                    not exists $old_contigs->{$scaf_ctg};
                if( $old_contigs->{$scaf_ctg} >= $self->min_contig_length ) {
                    $scaffold_lengths{$scaf_num} += $old_contigs->{$scaf_ctg};
                }
            }
        }
    }

    #TODO - this is pretty bad .. sorry will clean up
    my $new_scafs = {};
    #hash of scaffolds with array of contigs in scaffold as value
    #$new_scafs->{scaffold?}->{scaffold_contigs} = [
    #                                               contig??
    #                                               contig??
    #                                              ]
    my $scaf_lengths = {};
    #my %valid_contigs_to_export;
    #hash of scaffold name and scaffold size
    if ($scaffolds) {
	foreach my $scaf (@$scaffolds) {
	    $scaf =~ s/\s+//;
	    #TODO - don't differentiate between scaf with - and w/o .. no need
	    my @tmp = split (/-/, $scaf);
	    my $scaf_ctg_1;
	    foreach my $scaf_ctg (@tmp) {
		next if $scaf_ctg eq 'E'; #eg E-12.1-E
                my $scaf_num;
                if( $scaf_ctg =~ /^\d+\.\d+$/ ) {
                    ($scaf_num) = $scaf_ctg =~ /^(\d+)\./;
                }
                elsif ( $scaf_ctg =~ /^\d+$/ ) {
                    $scaf_num = $scaf_ctg;
                }
                else {
                    die "Could not get scaffold number from contig: $scaf_ctg\n";
                }
                if( not $scaffold_lengths{$scaf_num} or $scaffold_lengths{$scaf_num} < $self->min_scaffold_length ) {
                    #print "removing $scaf_num .. in scaffold with length ".$scaffold_lengths{$scaf_num}."\n";
                    delete $old_contigs->{$scaf_ctg} and next;
                }
                if ( $old_contigs->{$scaf_ctg} < $self->min_contig_length ) {
                    print "removing $scaf_ctg with length ".$old_contigs->{$scaf_ctg}."\n";
                    delete $old_contigs->{$scaf_ctg} and next;
                }
		$scaf_ctg_1 = $scaf_ctg unless $scaf_ctg_1; #sets first scaf contig
		push @{$new_scafs->{$scaf_ctg_1}->{scaffold_contigs}}, $scaf_ctg;
		$scaf_lengths->{$scaf_ctg_1} += $old_contigs->{$scaf_ctg};
		delete $old_contigs->{$scaf_ctg};
	    }
	}
    }


    #rename the remaining, non-scaffold contigs
    foreach my $contig (keys %$old_contigs) {
        if ( $old_contigs->{$contig} <= $self->min_contig_length ) {
            print "Excluding contig"."$contig with length: ".$old_contigs->{$contig}."\n";
            delete $old_contigs->{$contig};
            next;
        }
	push @{$new_scafs->{$contig}->{scaffold_contigs}}, $contig;
	$scaf_lengths->{$contig} = $old_contigs->{$contig};
	delete $old_contigs->{$contig};
    }

    #new scaffold numbers start with 0; new contig numbers start with 1
    #Contig0.1 is first scaffold, first contig

    my $new_scaf_num = 0;
    my $new_ctg_num = 1;
    my $new_scaf_names = {};

    #write a new gap file
    my $gap_file = $self->assembly_directory.'/edit_dir/msi.gap.txt';
    unlink $gap_file;
    my $gap_fh = Genome::Sys->open_file_for_writing( $gap_file );
    
    foreach my $scaf (sort {$scaf_lengths->{$b} <=> $scaf_lengths->{$a}} keys %{$scaf_lengths}) {
	foreach my $scaf_ctg ( @{$new_scafs->{$scaf}->{scaffold_contigs}} ) {
	    my $new_ctg_name = 'Contig'.$new_scaf_num.'.'.$new_ctg_num;
	    $new_scaf_names->{$scaf_ctg} = $new_ctg_name;
	    #print gap info to gap file only if part of multi contigs scaffold
	    if (scalar @{$new_scafs->{$scaf}->{scaffold_contigs}} > 1) {
		#dont' print gap size if last contig in scaffold
		next if $new_ctg_num == scalar @{$new_scafs->{$scaf}->{scaffold_contigs}};
		$gap_fh->print("$new_ctg_name 100\n");
	    }
	    $new_ctg_num++; #increment for next contig
	}
	$new_ctg_num = 1; #reset for next scaffold
	$new_scaf_num++;
    }

    $gap_fh->close;

    unless ( $new_scaf_names ) {
        $self->debug_message( "Could not get any new scaffolds .. check scaffolds file to make sure contigs match or set min-contig-length lower because all contigs may have been filtered out" );
        return;
    }

    return $new_scaf_names;
}

sub _write_new_scaffolds {
    my ( $self, $ace, $scaffold ) = @_;
    my $ace_out = $self->assembly_directory.'/edit_dir/ace.msi.scaffolds';
    unlink $ace_out;
    my $fh_out = Genome::Sys->open_file_for_writing( $ace_out );
    my $fh_in = Genome::Sys->open_file_for_reading( $ace );
    my $skip_this_contig = 0;

    my ( $contig_count, $read_count ) = 0;

    while (my $line = $fh_in->getline) {
        #stop writing when contig or wa tags are reached
        #only writing contigs here
        last if $line =~ /^CT\{$/ or $line =~ /^WA\{$/;
	if ($line =~ /^CO\s+/) {
	    chomp $line;
	    my ($contig_name) = $line =~ /^CO\s+(\S+)/;
	    my $rest_of_line = "$'";
	    my ($contig_number) = $contig_name =~ /contig(\S+)/i;
	    if (exists $scaffold->{$contig_number}) {
		$fh_out->print('CO '.$scaffold->{$contig_number}." $rest_of_line\n");
                $contig_count++;
                $skip_this_contig = 0;
	    }
	    else {
                #skip any contigs that are not specified in $scaffold
                $skip_this_contig = 1;
	    }
	}
        else {
	    $fh_out->print($line) unless $skip_this_contig == 1;
            $read_count++ if $line =~ /^RD\s+/;
	}
    }

    $fh_in->close;
    $fh_out->close;

    return $ace_out, $contig_count, $read_count;
}

#Transfer tags not used as of 7/1/11
sub _write_contig_tags {
    my ( $self, $ace, $scaffolds ) = @_;
    my $tags_out = $self->assembly_directory.'/edit_dir/ace.msi.tags';
    unlink $tags_out;
    my $fh_out = Genome::Sys->open_file_for_writing( $tags_out );
    my $fh = Genome::Sys->open_file_for_reading( $ace );
    my $in_tag_lines = 0;
    my $in_tag_comment_lines = 0;
    my $print_tag_lines = 0;
    while ( my $line = $fh->getline ) {
        if ( $line =~ /^CT\{/ ) {
            $in_tag_lines = 1;
            $print_tag_lines = 0;
            $in_tag_comment_lines = 0;
        }
        elsif ( $line =~ /^WA\{/ ) { #don't print WA tags
            $in_tag_lines = 0;
            $in_tag_comment_lines = 0;
            $print_tag_lines = 0;
        }
        elsif ( $in_tag_lines == 1 and $print_tag_lines == 1 and $line =~ /^COMMENT{/ ) {
            $in_tag_comment_lines = 1;
            $fh_out->print( $line );
        }
        elsif ( $in_tag_lines == 1 and $line =~ /^C{/ ) {
            $in_tag_comment_lines = 0;
            $fh_out->print( $line );
        }
        elsif ( $in_tag_comment_lines == 1 and $print_tag_lines == 1) {
            $fh_out->print( $line );
        }
        elsif ( $in_tag_lines == 1 and $line =~ /^contig(\S+)\s+/i and $in_tag_comment_lines == 0 ) {
            chomp $line;
	    my ($contig_name) = $line =~ /^(\S+)/;
	    my $rest_of_line = "$'";
	    my ($contig_number) = $contig_name =~ /contig(\S+)/i;
	    if (exists $scaffolds->{$contig_number}) {
                $fh_out->print( "CT{\n" );
		$fh_out->print( $scaffolds->{$contig_number}."$rest_of_line\n" );
                $print_tag_lines = 1;
	    }
	    else{
                $print_tag_lines = 0;
	    }
        }
        elsif ( $in_tag_lines == 1 and $print_tag_lines == 1 ) { #not printing anything
            $fh_out->print( $line );
        }
        #else #do nothing
    }
    
    $fh_out->close;
    $fh->close;

    return $tags_out;
}

#not used
sub _merge_files_get_contig_read_counts {
    my ( $self, $scaf_file, $tags_file ) = @_;
    my $ace_out = $self->assembly_directory.'/edit_dir/ace.msi.int1';
    unlink $ace_out;
    my $ace_out_fh = Genome::Sys->open_file_for_writing( $ace_out );
    my $s_fh = Genome::Sys->open_file_for_reading( $scaf_file );
    my $read_count = 0;
    my $contig_count = 0;
    while ( my $line = $s_fh->getline ) {
        next if $line =~ /^AS\s+/;
        $read_count++ if $line =~ /^RD\s+/;
        $contig_count++ if $line =~ /^CO\s+/;
        $ace_out_fh->print( $line );
    }
    $s_fh->close;
    my $t_fh = Genome::Sys->open_file_for_reading( $tags_file );
    while ( my $line = $t_fh->getline ) {
        $ace_out_fh->print( $line );
    }
    $t_fh->close;
    $ace_out_fh->close;

    unlink $scaf_file, $tags_file;

    return $ace_out, $contig_count, $read_count;
}

sub _update_ds_line_write_wa_tags {
    my ($self, $ace, $contig_count, $read_count ) = @_;
    my $fh = Genome::Sys->open_file_for_reading( $ace );
    my $ace_out = $self->assembly_directory.'/edit_dir/ace.msi';
    unlink $ace_out;
    my $out_fh = Genome::Sys->open_file_for_writing( $ace_out );
    #write as line
    $out_fh->print( "AS $contig_count $read_count\n\n" );
    while (my $line = $fh->getline) {
	if ($line =~ /^DS\s+/) {
	    if ($line =~ /PHD_FILE/) {
		$out_fh->print($line);
	    }
	    else {
		$line =~ s/DS /DS VERSION: 1 /;
		$out_fh->print($line);
	    }
	}
	else {
	    $out_fh->print($line);
	}
    }
    $fh->close;

    my $ball_dir = $self->assembly_directory.'/phdball_dir';
    if (-d $ball_dir) {
	my @phd_ball_file = glob ("$ball_dir/*");
	if (scalar @phd_ball_file > 0) {
	    foreach (@phd_ball_file) {
		$out_fh->print("\nWA{\n"."phdBall newbler 080416:144002\n".$_."\n}\n\n");
	    }
	}
    }

    $out_fh->close;
    unlink $ace;

    return $ace_out;
}

1
