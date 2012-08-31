
package Genome::Model::Tools::Varscan::ProcessSomatic;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::ProcessSomatic    Process somatic pipeline output
#
#    AUTHOR:     Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#    CREATED:    12/09/2009 by D.K.
#    MODIFIED:   12/29/2009 by D.K.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::ProcessSomatic {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <---
        status_file     => { is => 'Text', doc => "File containing varscan calls, e.g. status.varscan.snp" , is_optional => 0, is_input => 1},
        p_value_for_hc  => { is => 'Number', doc => "P-value threshold for high confidence", is_optional => 1, is_input => 1, default_value => '0.07'},
        max_normal_freq => { is => 'Number', doc => "Maximum normal frequency for HC Somatic", is_optional => 1, is_input => 1, default_value => '5'},
        min_tumor_freq  => { is => 'Number', doc => "Minimum tumor freq for HC Somatic", is_optional => 1, is_input => 1, default_value => '10'},
        report_only     => { is => 'Boolean', doc => "If set to 1, will not produce output files" , is_optional => 1, default_value => 0 },
        skip_if_output_present => { is => 'Boolean', doc => "If set to 1, will not run if output is present" , is_optional => 1, is_input => 1, default_value => 0},
        output_basename => { is => 'Text', doc => 'Base location for output files (e.g. /path/to/results/output.basename', is_optional => 1, is_input => 1 },

        #output filenames
        output_germline       => { is => 'Text', calculate_from => 'output_basename', calculate => q{ $output_basename . '.Germline' }, is_output => 1, },
        output_germline_hc    => { is => 'Text', calculate_from => 'output_germline', calculate => q{ $output_germline . '.hc' }, is_output => 1, },
        output_germline_lc    => { is => 'Text', calculate_from => 'output_germline', calculate => q{ $output_germline . '.lc' }, is_output => 1, },
        output_loh            => { is => 'Text', calculate_from => 'output_basename', calculate => q{ $output_basename . '.LOH' }, is_output => 1, },
        output_loh_hc         => { is => 'Text', calculate_from => 'output_loh', calculate => q{ $output_loh . '.hc' }, is_output => 1, },
        output_loh_lc         => { is => 'Text', calculate_from => 'output_loh', calculate => q{ $output_loh . '.lc' }, is_output => 1, },
        output_somatic        => { is => 'Text', calculate_from => 'output_basename', calculate => q{ $output_basename . '.Somatic' }, is_output => 1, },
        output_somatic_hc     => { is => 'Text', calculate_from => 'output_somatic', calculate => q{ $output_somatic . '.hc' }, is_output => 1, },
        output_somatic_lc     => { is => 'Text', calculate_from => 'output_somatic', calculate => q{ $output_somatic . '.lc' }, is_output => 1, },
        output_other          => { is => 'Text', calculate_from => 'output_basename', calculate => q{ $output_basename . '.other' }, is_output => 1, },
    ],
    has_param => [
        lsf_resource => { value => 'select[tmp>1000] rusage[tmp=1000]'},
    ]
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Process output from Varscan somatic"                 
}

sub help_synopsis {
    return <<EOS
    This command processes output from Varscan somatic (.snp or.indel), classifying variants by somatic status
    (Germline/Somatic/LOH) and by confidence (high/low).
    EXAMPLE:    gmt varscan process-somatic ...
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

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $status_file = $self->status_file;

    $self->output_basename($status_file) unless $self->output_basename;

    if(-e $status_file) {
        if($self->skip_if_output_present && -e $self->output_somatic_hc) {
            ## Skip ##
        } else {
            $self->process_results($status_file);
        }

    } else {
        $self->error_message("Status file $status_file not found!");
        die $self->error_message;
    }

    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub process_results {
    my $self = shift;
    my $variants_file = shift;
    my $file_header = "";

    my $report_only = $self->report_only;

    print "Processing variants in $variants_file...\n";

    my %variants_by_status = ();
    $variants_by_status{'Somatic'} = $variants_by_status{'Germline'} = $variants_by_status{'LOH'} = '';

    ## Parse the variant file ##

    my $input = Genome::Sys->open_file_for_reading($variants_file);
    my $lineCounter = 0;
    my %counters;

    my %filehandles;
    my %temp_paths;
    if(!$report_only) {
        for my $status ('Somatic', 'Germline', 'LOH') {
                my $status_file_accessor = 'output_' . lc($status);
                my $hc_file_accessor = $status_file_accessor . '_hc';
                my $lc_file_accessor = $status_file_accessor . '_lc';

                my ($status_fh, $status_temp) = Genome::Sys->create_temp_file();
                my ($high_confidence_fh, $high_confidence_temp) = Genome::Sys->create_temp_file();;
                my ($low_confidence_fh, $low_confidence_temp) = Genome::Sys->create_temp_file();;

                if($file_header) {
                    for my $fh ($status_fh, $high_confidence_fh, $low_confidence_fh) {
                        $fh->print($file_header,"\n");
                    }
                }

                $filehandles{$status}{status} = $status_fh;
                $temp_paths{$status}{status} = $status_temp;
                $filehandles{$status}{hc} = $high_confidence_fh;
                $temp_paths{$status}{hc} = $high_confidence_temp;
                $filehandles{$status}{lc} = $low_confidence_fh;
                $temp_paths{$status}{lc} = $low_confidence_temp;
        }

        my $other_fh = Genome::Sys->open_file_for_appending($self->output_other);
        $filehandles{other} = $other_fh;
    }

    while (<$input>) {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);

        if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")) {
            $file_header = $line;
            for my $status ('Somatic', 'Germline', 'LOH') {
                for my $fh (values %{ $filehandles{$status} }) {
                    $fh->print($file_header,"\n");
                }
            }
        } else {
            my $somatic_status = "";
            if($lineContents[13] && ($lineContents[13] =~ "Reference" || $lineContents[13] =~ "Somatic" || $lineContents[13] =~ "Germline" || $lineContents[13] =~ "Unknown" || $lineContents[13] =~ "LOH")) {
                $somatic_status = $lineContents[13];
            } elsif($lineContents[12] && ($lineContents[12] =~ "Reference" || $lineContents[12] =~ "Somatic" || $lineContents[12] =~ "Germline" || $lineContents[12] =~ "Unknown" || $lineContents[12] =~ "LOH")) {
                $somatic_status = $lineContents[12];
            } else {
                $self->warning_message("Unable to parse somatic_status from file $variants_file line $lineCounter");
                $somatic_status = "Unknown";
            }

            my $hc_or_lc = $self->process_line($line, $somatic_status, \%filehandles, \%temp_paths);
            $counters{$somatic_status}{$hc_or_lc}++;
        }
    }
    close($input);

    if(!$report_only) {
        for my $status ('Somatic', 'Germline', 'LOH') {
            $filehandles{$status}{status}->close;
            $filehandles{$status}{hc}->close;
            $filehandles{$status}{lc}->close;

            my $status_file_accessor = 'output_' . lc($status);
            my $hc_file_accessor = $status_file_accessor . '_hc';
            my $lc_file_accessor = $status_file_accessor . '_lc';

            Genome::Sys->copy_file($temp_paths{$status}{status}, $self->$status_file_accessor)
            if Genome::Sys->check_for_path_existence($temp_paths{$status}{status});
            Genome::Sys->copy_file($temp_paths{$status}{hc}, $self->$hc_file_accessor)
            if Genome::Sys->check_for_path_existence($temp_paths{$status}{status});
            Genome::Sys->copy_file($temp_paths{$status}{lc}, $self->$lc_file_accessor)
            if Genome::Sys->check_for_path_existence($temp_paths{$status}{status});
        }
        $filehandles{other}->close;
    }

    for my $status (keys %counters) {
        my $total = 0;
        my ($hc, $lc);
        for my $key (keys %{ $counters{$status} }) {
            $total += $counters{$status}{$key};
            if($key eq 'hc') { $hc = $counters{$status}{$key}; }
            if($key eq 'lc') { $lc = $counters{$status}{$key}; }
        }

        print "$total $status\n";
        if($hc) { print "\t$hc high confidence\n"; }
        if($lc) { print "\t$lc low confidence\n"; }
    }

    return 1;
}

sub process_line {
    my $self = shift;
    my $line = shift;
    my $status = shift;
    my $filehandles = shift;
    my $temp_paths = shift;

    my $report_only = $self->report_only;
    my $max_normal_freq = $self->max_normal_freq;
    my $min_tumor_freq = $self->min_tumor_freq;
    my $p_value_for_hc = $self->p_value_for_hc;

    my $hc_or_lc;
    ## Output Germline, Somatic, and LOH ##

    my ($status_fh, $high_confidence_fh, $low_confidence_fh);

    if($status eq "Germline" || $status eq "Somatic" || $status eq "LOH") {
        if(!$report_only) {
            $status_fh = $filehandles->{$status}{status};
            $high_confidence_fh = $filehandles->{$status}{hc};
            $low_confidence_fh = $filehandles->{$status}{lc};
        }

        my @lineContents = split(/\t/, $line);
        my $numContents = @lineContents;

        my $somatic_status = my $p_value = "";
        my $normal_reads2 = my $normal_freq = my $tumor_reads2 = my $tumor_freq = "";

        ## Get Somatic status and p-value ##

        for(my $colCounter = 4; $colCounter < $numContents; $colCounter++) {
            if($lineContents[$colCounter]) {
                my $value = $lineContents[$colCounter];

                if($value eq "Reference" || $value eq "Somatic" || $value eq "Germline" || $value eq "LOH" || $value eq "Unknown") {
                    $normal_reads2 = $lineContents[$colCounter - 7];
                    $normal_freq = $lineContents[$colCounter - 6];

                    $tumor_reads2 = $lineContents[$colCounter - 3];
                    $tumor_freq = $lineContents[$colCounter - 2];

                    $normal_freq =~ s/\%//;
                    $tumor_freq =~ s/\%//;

                    $somatic_status = $value;
                    $p_value = $lineContents[$colCounter + 1];
                    $p_value = $lineContents[$colCounter + 2] if($lineContents[$colCounter + 2] && $lineContents[$colCounter + 2] < $p_value);
                }
            }
        }

        ## Determine somatic status ##
        if(!$report_only) {
            ## Print to master status file ##
            $status_fh->print("$line\n");
            if($status eq "Somatic") {
                if($normal_freq <= $max_normal_freq && $tumor_freq >= $min_tumor_freq && $p_value <= $p_value_for_hc && $tumor_reads2 >= 2) {
                    $high_confidence_fh->print("$line\n");
                    $hc_or_lc = 'hc';
                } else {
                    $low_confidence_fh->print("$line\n");
                    $hc_or_lc = 'lc';
                }
            } elsif($status eq "Germline") {
                if($normal_freq >= $min_tumor_freq && $tumor_freq >= $min_tumor_freq && $p_value <= $p_value_for_hc) {
                    $high_confidence_fh->print("$line\n");
                    $hc_or_lc = 'hc';
                } else {
                    $low_confidence_fh->print("$line\n");
                    $hc_or_lc = 'lc';
                }
            } elsif($status eq "LOH") {
                if($normal_freq >= $min_tumor_freq && $p_value <= $p_value_for_hc) {
                    $high_confidence_fh->print("$line\n");
                    $hc_or_lc = 'hc';
                } else {
                    $low_confidence_fh->print("$line\n");
                    $hc_or_lc = 'lc';
                }
            }
        }

    } else {
        # Print all non-LOH/Somatic/Germline things to the "other" filehandle so we account for all variants that come into the tool
        if(!$report_only) {
            my $other_fh = $filehandles->{other};
            $other_fh->print("$line\n");
        }
        $hc_or_lc = 'n/a';
    }

    return $hc_or_lc;
}


1;
