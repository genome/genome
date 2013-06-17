package Genome::Model::Tools::BioSamtools::ErrorRate;

use strict;
use warnings;

use Genome;
use Statistics::R;
use Cwd;
use Bio::DB::Sam::Constants;

class Genome::Model::Tools::BioSamtools::ErrorRate {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A BAM format file of alignment data'
        },
        output_file => {
            is => 'Text',
            doc => 'A file path to store tab separated value output.  The file extension must be .tsv in pileup mode',
        },
        reference_fasta => {
            is_optional => 1,
        },
        pileup => {
            is => 'Boolean',
            doc => 'The pileup method takes longer, but provides an R graph as output.',
            default_value => 1,
            is_optional => 1,
        },
        use_c => {
            is => 'Boolean',
            doc => 'The C version of pileup is faster than the perl version.',
            default_value => 1,
            is_optional => 1,
        },
    ],
};

my @CIGAR_OPS = qw/M I D N S H P/;

sub execute {
    my $self = shift;

    if ($self->pileup) {
        if ($self->use_c) {
            return $self->c_pileup_error_rate;
        } else {
            return $self->perl_pileup_error_rate;
        }
    } else {
        return $self->non_pileup_error_rate;
    }
}

sub non_pileup_error_rate {
    my $self = shift;
    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    my $bam  = $refcov_bam->bio_db_bam;
    my $index = $refcov_bam->bio_db_index;
    my $header = $bam->header();
    my $text = $header->text;
    my @lines = split("\n",$text);
    my @rg_lines = grep {$_ =~ /^\@RG/} @lines;
    my %rg_libraries;
    for my $rg_line (@rg_lines) {
        unless ($rg_line =~ /ID\:(\S+)/) {
            die;
        }
        my $id = $1;
        my $lib;
        if ($rg_line =~ /LB\:(\S+)/) {
            $lib = $1;
        } elsif ($rg_line =~ /PU\:(\S+)/) {
            $lib = $1;
        } else {
            die;
        }
        $rg_libraries{$id} = $lib;
    }
    my $default_rg_id = 0;
    my @rg_ids = keys %rg_libraries;
    if (scalar(@rg_ids) == 1) {
        $default_rg_id = $rg_ids[0];
    }

    my $targets = $header->n_targets();
    my %read_groups;
    while (my $align = $bam->read1) {
        my $flag = $align->flag;
        my $read_group = $align->aux_get('RG');
        unless ($read_group) {
            if (defined($default_rg_id)) {
                $read_group = $default_rg_id;
            }
        }
        my $type;
        if ($flag & 1) {
            if ($flag & 64)  {
                $type = 'read_1';
            } else {
                $type = 'read_2';
            }
        } else {
            $type = 'fragment';
        }
        unless ($flag & 4) {
            my $match_descriptor = $align->aux_get('MD');
            my $x_descriptor = $align->aux_get('XD');
            unless ($match_descriptor) {
                $match_descriptor = $align->aux_get('XD');
            }
            if ($match_descriptor) {
                my @match_fields = split(/[^0-9]+/,$match_descriptor);
                #sum
                for my $matches (@match_fields) {
                    if ($matches =~ /^\s*$/) { next; }
                    $read_groups{$read_group}{$type}{matches} += $matches;
                }
                #character count minus ^
                my @mismatch_fields = split(/[0-9]+/,$match_descriptor);
                for my $mismatch_field (@mismatch_fields) {
                    $mismatch_field =~ s/\^//g;
                    $read_groups{$read_group}{$type}{mismatches} += length($mismatch_field);
                }
                #print $match_descriptor ."\t". $match_sum ."\t". $mismatch_sum ."\n";
            } else {
                print $align->qname ."\n";
                my $cigar_str = $align->cigar_str;
                print $cigar_str ."\n";
                my $c = $align->cigar;
                #die(Data::Dumper::Dumper($cigar));
                for (my $i=0; $i <scalar(@{$c}); $i++) {
                    my $op  = $c->[$i] & BAM_CIGAR_MASK;
                    my $len = $c->[$i] >> BAM_CIGAR_SHIFT;
                    print $CIGAR_OPS[$op] ."\t". $len ."\n";
                }
                die('Please implement cigar string parsing or use: gmt bio-samtools error-rate --pileup!');
            }
        } else {
            $read_groups{$read_group}{$type}{unaligned} += $align->l_qseq;
        }
    }
    print $output_fh "LIBRARY\tREAD_GROUP\tREAD_TYPE\tUNALIGNED\tALIGNED\tPC_ALIGNED\tMISMATCHES\tMATCHES\tERROR\n";
    for my $read_group (sort keys %read_groups) {
        for my $read_group_type (sort keys %{$read_groups{$read_group}}) {
            my $read_group_type_matches = $read_groups{$read_group}{$read_group_type}{matches};
            my $read_group_type_mismatches = $read_groups{$read_group}{$read_group_type}{mismatches};
            my $read_group_type_unaligned = $read_groups{$read_group}{$read_group_type}{unaligned};
            my $read_group_type_aligned = $read_group_type_matches + $read_group_type_mismatches;
            my $read_group_type_total = $read_group_type_aligned + $read_group_type_unaligned;
            my $read_group_type_align_rate = sprintf("%.02f",(($read_group_type_aligned / $read_group_type_total ) * 100));
            my $read_group_type_error_rate = sprintf("%.02f",(($read_group_type_mismatches / $read_group_type_aligned ) * 100));
            my $library = $rg_libraries{$read_group} || 'NA';
            print $output_fh $library ."\t". $read_group ."\t". $read_group_type ."\t".$read_group_type_unaligned ."\t". $read_group_type_aligned ."\t". $read_group_type_align_rate
                ."\t". $read_group_type_mismatches ."\t". $read_group_type_matches ."\t". $read_group_type_error_rate ."%\n";
        }
    }
    $output_fh->close;
    return 1;
}

sub c_pileup_error_rate {
    my $self = shift;

    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->output_file,qw/.tsv/);
    unless ($suffix eq '.tsv') {
        die('Unable to parse output *.tsv file name: '. $self->output_file);
    }

    my $bam_file = $self->bam_file;
    my $output_file = $self->output_file;
    # C util is running out of my home directory until it has been tested and can be deployed to the blades
    my $cmd = "samtools view $bam_file | " . Genome::Sys->swpath("bam-errorrate","0.6") . " > $output_file";
    my $temporary_errorrate = '/gscuser/iferguso/bin/bam-errorrate0.7';
    if (-e $temporary_errorrate) {
        $cmd = "samtools view $bam_file | $temporary_errorrate > $output_file";
    }
    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => [ $bam_file ],
        output_files => [ $output_file ],
        skip_if_output_is_present => 0,
    );

    $self->run_r_script;
    return 1;
}

sub run_r_script {
    my $self = shift;

    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->output_file,qw/.tsv/);

    my %read_ends;
    open my $tsv_fh, '<', $self->output_file;

    my $confirmed_header = 0;

    while (my $line = <$tsv_fh>) {
        chomp $line;
        my @fields = split /\s+/, $line;

        if (not $confirmed_header) {
            if ($fields[0] ne 'read_end') {
                die sprintf "Unable to parse header in '%s'. Expected first field to be read_end, but it was '%s'.",
                    $self->output_file,
                    $fields[0];
            }
            $confirmed_header = 1;
            next;
        }

        if ($fields[0] !~ /^[012]$/) {
            die sprintf "Read end was '%s', which was not 0, 1, or 2.", $fields[0];
        }
        $read_ends{$fields[0]} = 1;
    }
    close $tsv_fh;

    my $r_library = $self->__meta__->module_path;
    $r_library =~ s/\.pm/\.R/;
    $self->status_message('R LIBRARY: '. $r_library);
    my $tempdir = Genome::Sys->create_temp_directory();

    my $cwd = getcwd();

    for my $read_end (keys %read_ends) {
        my $input_file = $self->output_file;

        my $count_plot_file = $dirname . $basename .'_counts_'. $read_end .'.png';
        my $rate_plot_file = $dirname . $basename .'_rates_'. $read_end .'.png';
        my $rate_dist_file = $dirname . $basename .'_rate_distribution_'. $read_end .'.png';

        # The rate plot is all that is necessary at this time
        #my $r_cmd = "generatePlots('$input_file','$read_end','$count_plot_file','$rate_plot_file','$rate_dist_file')";

        my $r_script = $tempdir .'/'. $basename .'_error_rate_'. $read_end .'.R';
        my $r_script_fh = Genome::Sys->open_file_for_writing($r_script);

        print $r_script_fh "source('$r_library')\n";
        print $r_script_fh "fullFile <- readTable('$input_file')\n";
        print $r_script_fh "readEndErrorRate <- getReadEnd(fullFile,'$read_end')\n";
        print $r_script_fh "positionErrorRate <- getPositionErrorRate(readEndErrorRate)\n";
        print $r_script_fh "makeRatePlot(positionErrorRate,'$rate_plot_file','$read_end')\n";
        #print $r_script_fh "$r_cmd\n";
        $r_script_fh->close;

        my $cmd = 'Rscript '. $r_script;
        Genome::Sys->shellcmd(
            cmd => $cmd,
        );
        unlink($r_script);
    }
    chdir $cwd;
    return 1;
}

sub perl_pileup_error_rate {
    my $self = shift;

    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->output_file,qw/.tsv/);
    unless ($suffix eq '.tsv') {
        die('Unable to parse output *.tsv file name: '. $self->output_file);
    }

    my $fai = Bio::DB::Sam::Fai->load($self->reference_fasta);
    unless ($fai) {
        die('Failed to load fai for: '. $self->reference_fasta);
    }

    my $bam = Bio::DB::Bam->open($self->bam_file);
    unless ($bam) {
        die('Failed to load BAM file: '. $self->bam_file);
    }

    my $index  = Bio::DB::Bam->index($self->bam_file);
    unless ($index) {
        die('Failed to load BAM index for: '. $self->bam_file);
    }

    my $header = $bam->header;
    my $tid_names = $header->target_name;
    my $tid_lengths = $header->target_len;

    my %read_counts;
    my %del_pileups;
    my $callback = sub {
        my ($tid,$pos,$pileups,$callback_data) = @_;
        my $chr = $tid_names->[$tid];
        my $ref_base = $fai->fetch($chr .':'. ($pos+1) .'-'. ($pos+1) );
        for my $pileup (@$pileups) {
            my $b = $pileup->alignment;
            my $qname = $b->qname;
            my $flag = $b->flag;
            # Unmapped
            if ($flag & 4) {
                next;
            }
            if ($pileup->is_refskip) {
                # There is a gap in the alignment, ignore gaps.
                next;
            }
            my $read_end = 0;
            if ($flag & 1) {
                if ($flag & 64) {
                    $read_end = 1;
                } elsif ($flag & 128) {
                    $read_end = 2;
                }
            }
            $qname .= '/'. $read_end;

            my $qpos = $pileup->qpos;
            my $qbase = substr($b->qseq,$qpos,1);
            unless ($qbase) {
                warn('Failed to get query base for '. $b->qname .' at position '. $qpos .' from string: '. $b->qseq);
                next;
            }
            my $offset = 0;
            my $c = $b->cigar;
            if ($flag & 16) {
                # Hard-masking really jacks up using the cigar string
                # my $len = $b->cigar2qlen
                # TODO: Can we simply flip the position for the negative strand alignments
                my $length = $b->l_qseq;
                $qpos = ($length - $qpos - 1);

                my $op  = $c->[-1] & BAM_CIGAR_MASK;
                if ($CIGAR_OPS[$op] eq 'H') {
                    $offset = $c->[-1] >> BAM_CIGAR_SHIFT;
                }
            } else {
                my $op  = $c->[0] & BAM_CIGAR_MASK;
                if ($CIGAR_OPS[$op] eq 'H') {
                    $offset = $c->[0] >> BAM_CIGAR_SHIFT;
                }
            }
            $qpos += $offset;
            $read_counts{$read_end}{total}->[$qpos]++;

            my $indel = $pileup->indel;
            if ($indel) {
                my $indel_size = abs($indel);
                if ($indel > 0) {
                    for (1 .. $indel_size) {
                        # TODO: Verify if we need to add or subtract based on strand of alignment
                        my $ipos;
                        if ($flag & 16) {
                            $ipos = ($qpos - $_)
                        } else {
                            $ipos = ($qpos + $_)
                        }
                        $read_counts{$read_end}{insertion}->[$ipos]++;
                        $read_counts{$read_end}{total}->[$ipos]++;
                    }
                } else {
                    $read_counts{$read_end}{deletion}->[$qpos] += $indel_size;
                    for (1 .. $indel_size) {
                        # This is relative to the reference position
                        my $dpos = ($pos + $_);
                        $del_pileups{$qname}{$dpos} = 1;
                    }
                }
            } elsif ($del_pileups{$qname}{$pos}) {
                my $del_pileup = delete($del_pileups{$qname}{$pos});
                $read_counts{$read_end}{total}->[$qpos]--;
                #print $pos ."\t". $ref_base ."\t". $qname ."\t". $del_pileup ."\t". $del_pileup ."\t". $pileup->qpos ."\t". $qbase ."\t". $pileup->indel ."\n";
                next;
            }
            if ($qbase =~ /[nN]/) {
                $read_counts{$read_end}{ambiguous}->[$qpos]++;
            } elsif (uc($ref_base) ne uc($qbase)) {
                $read_counts{$read_end}{mismatch}->[$qpos]++;
            } else {
                $read_counts{$read_end}{match}->[$qpos]++;
            }
        }
    };
    for (my $i = 0; $i < scalar(@{$tid_lengths}); $i++) {
        my $end = $tid_lengths->[$i];
        $index->pileup($bam,$i,'0',$end,$callback);
    }
    my @headers = qw/read_end position total match error error_rate mismatch mismatch_rate ambiguous ambiguous_rate insertion insertion_rate deletion deletion_rate/;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        output => $self->output_file,
    );
    unless ($writer) {
        die('Failed to create output writer!');
    }
    for my $read_end (sort keys %read_counts) {
        my @positions = @{$read_counts{$read_end}{total}};
        my $sum_total;
        my $sum_match;
        my $sum_mismatch;
        my $sum_ambiguous;
        my $sum_insertion;
        my $sum_deletion;
        my $sum_error;
        for (my $i = 0; $i < scalar(@positions); $i++) {
            my $position_count = $positions[$i];
            if (!$position_count) {
                my %data = (
                    read_end => $read_end,
                    position => $i,
                    total => 0,
                    match => 0,
                    error => 0,
                    error_rate => 0,
                    mismatch => 0,
                    mismatch_rate => 0,
                    ambiguous => 0,
                    ambiguous_rate => 0,
                    insertion => 0,
                    insertion_rate => 0,
                    deletion => 0,
                    deletion_rate => 0,
                );
                $writer->write_one(\%data);
                next;
            }
            my $match = $read_counts{$read_end}{match}->[$i] || 0;
            my $mismatch = $read_counts{$read_end}{mismatch}->[$i] || 0;
            my $ambiguous = $read_counts{$read_end}{ambiguous}->[$i] || 0;
            my $insertion = $read_counts{$read_end}{insertion}->[$i] || 0;
            my $deletion = $read_counts{$read_end}{deletion}->[$i] || 0;
            
            my $total = $match + $mismatch + $ambiguous + $insertion;
            my $error = $mismatch + $ambiguous + $insertion + $deletion;
            my $error_rate = $error / $total;
            my $mismatch_rate = $mismatch / $total;
            my $ambiguous_rate = $ambiguous / $total;
            my $insertion_rate = $insertion / $total;
            my $deletion_rate = $deletion / $total;
            my %data = (
                read_end => $read_end,
                position => $i,
                total => $total,
                match => $match,
                error => $error,
                error_rate => $error_rate,
                mismatch => $mismatch,
                mismatch_rate => $mismatch_rate,
                ambiguous => $ambiguous,
                ambiguous_rate => $ambiguous_rate,
                insertion => $insertion,
                insertion_rate => $insertion_rate,
                deletion => $deletion,
                deletion_rate => $deletion_rate,
            );
            $writer->write_one(\%data);
            $sum_total += $total;
            $sum_match += $match;
            $sum_mismatch += $mismatch;
            $sum_ambiguous += $ambiguous;
            $sum_insertion += $insertion;
            $sum_deletion += $deletion;
            $sum_error += $error;
        }
        my $error_rate = $sum_error / $sum_total;
        my $mismatch_rate = $sum_mismatch / $sum_total;
        my $ambiguous_rate = $sum_ambiguous / $sum_total;
        my $insertion_rate = $sum_insertion / $sum_total;
        my $deletion_rate = $sum_deletion / $sum_total;
        my %data = (
            read_end => $read_end,
            position => 'SUM',
            total => $sum_total,
            match => $sum_match,
            error => $sum_error,
            error_rate => $error_rate,
            mismatch => $sum_mismatch,
            mismatch_rate => $mismatch_rate,
            ambiguous => $sum_ambiguous,
            ambiguous_rate => $ambiguous_rate,
            insertion => $sum_insertion,
            insertion_rate => $insertion_rate,
            deletion => $sum_deletion,
            deletion_rate => $deletion_rate,
        );
        $writer->write_one(\%data);
    }
    # Close the fh
    $writer->output->close;
    $self->run_r_script;
    return 1;
}

1;
