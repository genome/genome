package Genome::Model::Tools::SnpArray::MergeSegments;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeSegments - merges adjoining segments of similar copy number; distinguishes amplifications and deletions
#
#       AUTHOR:         Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#       CREATED:        04/01/2009 by D.K.
#       MODIFIED:       04/01/2009 by D.K.
#
#       NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();

class Genome::Model::Tools::SnpArray::MergeSegments {
    is => 'Command',
    has => [ 
        segments_file       => { 
            is => 'Text', 
            doc => "Segments with p-values from running CBS on data", 
            is_optional => 0, 
            is_input => 1 },

        amp_threshold       => { is => 'Text', doc => "Minimum seg_mean threshold for amplification", is_optional => 1, is_input => 1, default => 0.25},
        del_threshold       => { is => 'Text', doc => "Maximum seg_mean threshold for deletion", is_optional => 1, is_input => 1, default => -0.25},
        size_threshold      => { is => 'Text', doc => "Fraction of chromosome arm length above which an event is considered large-scale", is_input => 1, default => 0.25},
        min_num_mark        => { is => 'Text', doc => "Minimum number of markers to include a segment (higher = more stringent); 30 is a good start", is_optional => 0, is_input => 1, default => 1},
        output_basename     => { is => 'Text', doc => "Base name for output", is_optional => 1, is_input => 1},
        ref_arm_sizes       => { is => 'Text', doc => "Two column file of reference name and size in bp for calling by chromosome arm", example_values => ["/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/chromosome_arm_coordinates.tsv"]},
        verbose     => { is => 'Text', doc => "If set to 1, use for verbose output", is_optional => 1, is_input => 1},
        ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges adjoining segments of similar copy number"
}

sub help_synopsis {
    return <<EOS
        This command merges merges adjoining CBS segments of similar copy number and distinguishes amplifications and deletions
EXAMPLE:        gmt snp-array merge-adjoining-segments --segments-file snpArray.cbs.segments.tsv --output-basename snpArray.cbs.segments-merged
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute
{                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $segments_file = $self->segments_file;
    my $ref_sizes_file = $self->ref_arm_sizes;
    my $output_basename = $self->output_basename;

    ## Get thresholds ##

    my $amp_threshold = $self->amp_threshold;
    my $del_threshold = $self->del_threshold;
    my $size_threshold = $self->size_threshold;

    ## Load the sizes of chromosome arms ##

    my %ref_sizes = parse_ref_sizes($ref_sizes_file);
    my %stats = ();

    ## Reset various statistic counters ##
    $stats{'num_amp_del_segments'} = $stats{'num_merged_events'} = 0;
    $stats{'amplification'} = $stats{'deletion'} = 0;
    $stats{'focal amplification'} = $stats{'focal deletion'} = 0;
    $stats{'large-scale amplification'} = $stats{'large-scale deletion'} = 0;

    ## Parse the segments file ##

    my $input = new FileHandle ($segments_file);
    my $lineCounter = 0;

    my @merged_events = ();


    my $current_chrom = my $current_chr_start = my $current_chr_stop = my $current_event_type = "";
    my $current_list = "";
    my $current_segments = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;
        $line =~ s/"//g;
        my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\s+/, $line);

        if($id eq "ID")
        {
            ## Skip header ##
        }
        else
        {
            $chrom = "X" if($chrom eq "23");
            $chrom = "Y" if($chrom eq "24");
            $chrom =~ s/\"//g;

            ## Determine event type ##
            my $event_type = "neutral";

            if($seg_mean >= $amp_threshold)
            {
                $event_type = "amplification";
                $stats{'num_amp_del_segments'}++;
            }
            elsif($seg_mean <= $del_threshold)
            {
                $event_type = "deletion";
                $stats{'num_amp_del_segments'}++;
            }
            else
            {
                $stats{'num_neutral_segments'}++;
            }

            ## PROCEED ON JOINING DECISION ##

            ## Case 1: No Events Yet Parsed ##

            if(!$current_chrom)
            {
                ## Start a new event ##
                $current_chrom = $chrom;
                $current_chr_start = $chr_start;
                $current_chr_stop = $chr_stop;
                $current_event_type = $event_type;
                $current_segments = 1;
                $current_list = $line;
            }

            ## Case 2: Chromosome or event type changes ##

            elsif ($chrom ne $current_chrom || $event_type ne $current_event_type)
            {
                if($current_segments)
                {
                    $merged_events[$stats{'num_merged_events'}] = process_event($current_chrom, $current_chr_start, $current_chr_stop, $current_event_type, $current_segments, $current_list);
                    $stats{'num_merged_events'}++;
                }

#                               if($event_type ne "neutral")
#                               {
                ## Start new region ##

                $current_chrom = $chrom;
                $current_chr_start = $chr_start;
                $current_chr_stop = $chr_stop;
                $current_event_type = $event_type;
                $current_segments = 1;
                $current_list = $line;
#                               }
#                               else
#                               {
#                                       $current_chrom = $current_chr_start = $current_chr_stop = $current_event_type = $current_list = "";
#                                       $current_segments = 0;
#                               }
            }

            ## Case 3: Same type of event ##

            elsif($event_type eq $current_event_type)
            {
                $current_chr_stop = $chr_stop;
                $current_list .= "\n" . $line;
                $current_segments++;
            }
        }
    }

    close($input);

    if($current_segments)
    {
        $merged_events[$stats{'num_merged_events'}] = process_event($current_chrom, $current_chr_start, $current_chr_stop, $current_event_type, $current_segments, $current_list);
        $stats{'num_merged_events'}++;
    }


#       print "$stats{'num_amp_del_segments'} copy-change segments\n";
#       print "$stats{'num_merged_events'} merged events\n";

    my $large_scale_amp_arms = my $large_scale_del_arms = "";

    ## Open outfile for amps and dels ##

    open(EVENTS, ">$output_basename.events.tsv") or die "Can't open outfile: $!\n";
    print EVENTS "chrom\tchr_start\tchr_stop\tseg_mean\tnum_segments\tnum_markers\tp_value\tevent_type\tevent_size\tsize_class\tchrom_arm\tarm_fraction\tchrom_fraction\n";
    ## Process the merged segments ##

    foreach my $line (@merged_events)
    {
        my ($chrom, $chr_start, $chr_stop, $avg_seg_mean, $event_type, $num_segments, $num_mark, $p_value) = split(/\t/, $line);
        my $event_size = $chr_stop - $chr_start + 1;

        ## Save event type ##

        $stats{$event_type}++;

        ## Determine size category ##

        my $size_category = "focal";
        my $event_detail = "";
        my $arm_name = my $chrom_arm_fraction = my $chrom_fraction = "";

        ## Find most-affected chromosome arm ##

        my $p_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tp"}, $chr_start, $chr_stop);
        my $q_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tq"}, $chr_start, $chr_stop);
        my ($q_start, $q_stop) = split(/\t/, $ref_sizes{"$chrom\tq"});
        my $chrom_size = $q_stop;

        if($q_arm_fraction > $p_arm_fraction)
        {
            $arm_name = $chrom . "q";
            $chrom_arm_fraction = $q_arm_fraction;
        }
        else
        {
            $arm_name = $chrom . "p";
            $chrom_arm_fraction = $p_arm_fraction;
        }

        ## Determine fraction of chromosome affected ##

        $chrom_fraction = $event_size / $chrom_size;

        ## If >50% of arm or >25% of chromosome affected, call it large-scale ##

        if(($chrom_arm_fraction && $chrom_arm_fraction >= $size_threshold) || ($chrom_fraction && $chrom_fraction >= ($size_threshold / 2)))
        {
            $size_category = "large-scale";
            $event_detail = join("\t", $chrom . $arm_name, sprintf("%.2f", $chrom_arm_fraction * 100) . "%", sprintf("%.2f", $chrom_fraction * 100) . "%");

        }

        if($num_mark >= $self->min_num_mark)
        {
            ## Save the event ##
            $stats{"$size_category $event_type"}++;

            ## Save which arm ##
            $stats{"$size_category $event_type arms"} .= "," if($stats{"$size_category $event_type arms"});
            $stats{"$size_category $event_type arms"} .= $arm_name;

            $chrom_arm_fraction = sprintf("%.2f", $chrom_arm_fraction * 100) . '%';
            $chrom_fraction = sprintf("%.2f", $chrom_fraction * 100) . '%';

            print EVENTS join("\t", $chrom, $chr_start, $chr_stop, $avg_seg_mean, $num_segments, $num_mark, $p_value, $event_type, $event_size, $size_category, $arm_name, $chrom_arm_fraction, $chrom_fraction) . "\n";
        }
        else
        {
            $stats{"$event_type filtered"}++;
        }

    }

    close(EVENTS);



    ## Print summary stats to file ##
    $stats{'large-scale amplification arms'} = "NA" if(!$stats{'large-scale amplification arms'});
    $stats{'large-scale deletion arms'} = "NA" if(!$stats{'large-scale deletion arms'});

    open(SUMMARY, ">$output_basename.summary.txt") or die "Can't open outfile: $!\n";
    print SUMMARY "segments\tmerged_events\tamps\tlarge-scale\tfocal\tdels\tlarge-scale\tfocal\tamp_regions\tdel_regions\n";

    print SUMMARY join("\t", $stats{'num_amp_del_segments'}, $stats{'num_merged_events'}, $stats{'amplification'}, $stats{'large-scale amplification'}, $stats{'focal amplification'}, $stats{'deletion'}, $stats{'focal deletion'}, $stats{'large-scale deletion'}, $stats{'large-scale amplification arms'}, $stats{'large-scale deletion arms'}) . "\n";

    close(SUMMARY);

    if($self->verbose)
    {
        print $stats{'amplification'} . " classified as amplifications\n";
        print $stats{'amplification filtered'} . " removed by min-num-mark requirement\n";
        print $stats{'large-scale amplification'} . " large-scale amplifications ";
        print "(" . $stats{'large-scale amplification arms'} . ")" if($stats{'large-scale amplification'});
        print "\n";
        print $stats{'focal amplification'} . " focal amplifications\n";

        print $stats{'deletion'} . " classified as deletions\n";
        print $stats{'deletion filtered'} . " removed by min-num-mark requirement\n";
        print $stats{'large-scale deletion'} . " large-scale deletions ";
        print "(" . $stats{'large-scale deletion arms'} . ")" if($stats{'large-scale deletion'});
        print "\n";
        print $stats{'focal deletion'} . " focal deletions\n";
    }
    else
    {
        print join("\t", $stats{'large-scale amplification'}, $stats{'large-scale deletion'}, $stats{'focal amplification'}, $stats{'focal deletion'}, "(" . $stats{'amplification filtered'} . ")", "(" . $stats{'deletion filtered'}) . ")" . "\n";
    }



}




################################################################################################
# Merge Events - merge consecutive events into single ones
#
################################################################################################

sub process_event
{
    my ($event_chrom, $event_chr_start, $event_chr_stop, $event_type, $num_segments, $segment_list) = @_;

    ## Split merged events and calculate averages ##

    my @events = split(/\n/, $segment_list);

    my $total_num_mark = my $seg_mean_sum = my $p_val_sum = my $p_val_num = 0;

    foreach my $line (@events)
    {
        my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\s+/, $line);
        $total_num_mark += $num_mark;
        $seg_mean_sum += $seg_mean;
        $p_val_sum += $p_value if($p_value ne "NA");
        $p_val_num++;
    }

    ## Calculate averages ##

    my $avg_seg_mean = $seg_mean_sum / $num_segments;
    my $avg_p_value = "NA";
    $avg_p_value = $p_val_sum / $p_val_num if($p_val_num);

    ## Calculate region_size ##

    my $event_size = $event_chr_stop - $event_chr_start + 1;

    return(join("\t", $event_chrom, $event_chr_start, $event_chr_stop, $avg_seg_mean, $event_type, $num_segments, $total_num_mark, $avg_p_value));

}



################################################################################################
# Calculate arm fraction - get fraction of arm consumed by copy number event
#
################################################################################################

sub calculate_arm_fraction
{
    my ($ref_sizes, $chr_start, $chr_stop) = @_;
    my ($arm_start, $arm_stop) = split(/\t/, $ref_sizes);
    my $arm_size = $arm_stop - $arm_start + 1;

    ## Check to see if they overlap ##

    if($chr_start <= $arm_stop && $chr_stop >= $arm_start)
    {
        ## Determine bounds of this event on this chromosome arm ##
        my $event_arm_start = $chr_start;
        my $event_arm_stop = $chr_stop;
        $event_arm_start = $arm_start if($event_arm_start < $arm_start);
        $event_arm_stop = $arm_stop if($event_arm_stop > $arm_stop);
        my $event_arm_size = $event_arm_stop - $event_arm_start + 1;

        if($event_arm_size)
        {
            my $fraction_of_arm = $event_arm_size / $arm_size;
        }
        else
        {
            return(0);
        }
    }
    else
    {
        return(0);
    }
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_ref_sizes
{                               # replace with real execution logic.
    my $FileName = shift(@_);

    my %arms = ();

    my $input = new FileHandle ($FileName);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my ($chrom, $chr_start, $chr_stop, $arm_name) = split(/\s+/, $line);
        $arms{$chrom . "\t" . $arm_name} = join("\t", $chr_start, $chr_stop);
    }

    close($input);

    return(%arms);
}


###############################################################################
# commify - add appropriate commas to long integers
###############################################################################

sub commify
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}





################################################################################################
# Execute - the main program logic
#
################################################################################################

sub old_execute
{                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $segments_file = $self->segments_file;
    my $ref_sizes_file = $self->ref_arm_sizes;
    my $output_basename = $self->output_basename;

    ## Get thresholds ##

    my $amp_threshold = $self->amp_threshold;
    my $del_threshold = $self->del_threshold;
    my $size_threshold = $self->size_threshold;

    my %ref_sizes = parse_ref_sizes($ref_sizes_file);
    my %stats = ();
    $stats{'num_segments'} = $stats{'num_amps'} = $stats{'num_dels'} = $stats{'num_neutral'} = 0;
    $stats{'total_bp'} = $stats{'amp_bp'} = $stats{'del_bp'} = $stats{'neutral_bp'} = 0;


    ## Parse file and merge similar events ##
    my $merged_file = "$output_basename.merged.tsv";
    merge_events($segments_file, $merged_file, $amp_threshold, $del_threshold);
    return(0);

    ## Open outfile for amps and dels ##

    open(EVENTS, ">$output_basename.events.tsv") or die "Can't open outfile: $!\n";


    ## Parse the segments file ##

    my $input = new FileHandle ($segments_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\t/, $line);

        if($chrom && $chrom eq "chrom")
        {
            ## Print header
            print EVENTS "$line\tevent_size\tevent_type\n";
        }
        elsif($chrom && $chrom ne "chrom")
        {
            $stats{'num_segments'}++;

            $chrom = "X" if($chrom eq "23");
            $chrom = "Y" if($chrom eq "24");

            ## Determine size ##

            my $event_size = $chr_stop - $chr_start + 1;
            $stats{'total_bp'} += $event_size;

            ## Determine size category ##

            my $size_category = "focal";
            my $event_detail = "";

            my $arm_name = my $chrom_arm_fraction = "";

            my $p_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tp"}, $chr_start, $chr_stop);
            my $q_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tq"}, $chr_start, $chr_stop);
            my ($q_start, $q_stop) = split(/\t/, $ref_sizes{"$chrom\tq"});
            my $chrom_size = $q_stop;

            if($q_arm_fraction > $p_arm_fraction)
            {
                $arm_name = "p";
                $chrom_arm_fraction = $q_arm_fraction;
            }
            else
            {
                $arm_name = "q";
                $chrom_arm_fraction = $p_arm_fraction;
            }

            my $chrom_fraction = $event_size / $chrom_size;

            if(($chrom_arm_fraction && $chrom_arm_fraction >= $size_threshold) || ($chrom_fraction && $chrom_fraction >= ($size_threshold / 2)))
            {
                $size_category = "large-scale";
                $event_detail = join("\t", $chrom . $arm_name, sprintf("%.2f", $chrom_arm_fraction * 100) . "%", sprintf("%.2f", $chrom_fraction * 100) . "%");
            }
            elsif($chrom_fraction > 0.25)
            {
                warn "Event $chrom $chr_start-$chr_stop ($event_size bp) has $p_arm_fraction of p, $q_arm_fraction of q, and $chrom_fraction of $chrom\n";
            }


            ## Determine copy number class ##

            my $copy_class = "neutral";

            if($seg_mean >= $amp_threshold)
            {
                $copy_class = "amplification";
                $stats{'num_amps'}++;
                $stats{'amp_bp'} += $event_size;
                $stats{"$size_category\t$copy_class"}++;
                print EVENTS join("\t", $line, $size_category, $copy_class, $event_detail) . "\n";
            }
            elsif($seg_mean <= $del_threshold)
            {
                $copy_class = "deletion";
                $stats{'num_dels'}++;
                $stats{'del_bp'} += $event_size;
                $stats{"$size_category\t$copy_class"}++;
                print EVENTS join("\t", $line, $size_category, $copy_class, $event_detail) . "\n";
            }
            else
            {
                $stats{'num_neutral'}++;
                $stats{'neutral_bp'} += $event_size;
            }

#                       print join("\t", $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $p_value, $size_category, $copy_class) . "\n";
        }


    }
    close($input);

    ## Determine base pair fractions ##
    my $pct_amp_bp = my $pct_del_bp = my $pct_neutral_bp = "-";
    $pct_amp_bp = sprintf("%.2f", $stats{'amp_bp'} / $stats{'total_bp'} * 100) . '%';
    $pct_del_bp = sprintf("%.2f", $stats{'del_bp'} / $stats{'total_bp'} * 100) . '%';
    $pct_neutral_bp = sprintf("%.2f", $stats{'neutral_bp'} / $stats{'total_bp'} * 100) . '%';


    ## Reset counters ##

    $stats{"large-scale\tamplification"} = 0 if(!$stats{"large-scale\tamplification"});
    $stats{"large-scale\tdeletion"} = 0 if(!$stats{"large-scale\tdeletion"});

    open(SUMMARY, ">$output_basename.summary.txt") or die "Can't open summary file: $!\n";
    print SUMMARY "segments\tneutral\tnum_neutral_bp\tpct_neutral_bp\t";
    print SUMMARY "amplifications\tlarge_scale_amps\tfocal_amps\tnum_amp_bp\tpct_amp_bp\t";
    print SUMMARY "deletions\tlarge_scale_dels\tfocal_dels\tnum_del_bp\tpct_del_bp\t";
    print SUMMARY "\n";

    print SUMMARY join("\t", $stats{'num_segments'}, $stats{'num_neutral'}, $stats{'neutral_bp'}, $pct_neutral_bp) . "\t";
    print SUMMARY join("\t", $stats{'num_amps'}, $stats{"large-scale\tamplification"}, $stats{"focal\tamplification"}, $stats{'amp_bp'}, $pct_amp_bp) . "\t";
    print SUMMARY join("\t", $stats{'num_dels'}, $stats{"large-scale\tdeletion"}, $stats{"focal\tdeletion"}, $stats{'del_bp'}, $pct_del_bp);
    print SUMMARY "\n";


#       print "$lineCounter lines parsed\n";
    print $stats{'num_segments'} . " segments\n";
    print $stats{'num_neutral'} . " ($pct_neutral_bp bp) were neutral\n";

    print $stats{'num_amps'} . " ($pct_amp_bp bp) classified as amplifications\n";
    print "\t" . $stats{"large-scale\tamplification"} . " large-scale events\n";
    print "\t" . $stats{"focal\tamplification"} . " focal events\n";

    print $stats{'num_dels'} . " ($pct_del_bp bp) classified as deletions\n";
    print "\t" . $stats{"large-scale\tdeletion"} . " large-scale events\n";
    print "\t" . $stats{"focal\tdeletion"} . " focal events\n";


    close(SUMMARY);
}




1;


