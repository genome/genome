package Genome::Model::Tools::Cmds::CreateOutputTable;

use warnings;
use strict;
use Genome;
use IO::File;

class Genome::Model::Tools::Cmds::CreateOutputTable {
    is => 'Command',
    has => [
    region_call_dir => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'Directory containing region calls output by gmt cmds individual-region-calls',
    },
    output_file => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        is_output => 1,
        doc => 'Output filename',
    },
    pval_cutoff => {
        type => 'Number',
        is_optional => 1,
        default => '0.05',
        doc => 'P-value cutoff. If P-value is below this value, an Amp or Del call will be made.',
    },
    amp_threshold => {
        type => 'Number',
        is_optional => 1,
        doc => 'If specified, regions with copy number above this value will be called Amp (instead of p-value)',
    },
    del_threshold => {
        type => 'Number',
        is_optional => 1,
        doc => 'If specified, regions with copy number below this value will be called Del (instead of p-value)',
    },
    ]
};

sub help_brief {
    'Make Amp, Del, or Neutral calls from cmds region calls'
}

sub help_detail {
    "This script looks at each region call file created by gmt cmds individual-region-calls (except ROIs.txt) and judges the p-value and sd value to make a call of either Amp, Del, or Neutral. If sd > 0 and p-value < cutoff, the call is Amp. If sd < 0 and p-value < cutoff, the call is Del. Any other case: call is Neutral."
}

sub execute {
    my $self = shift;
    my $region_call_dir = $self->region_call_dir;
    my $outfile = $self->output_file;
    my $pval_cutoff = $self->pval_cutoff;

    #open output filehandle
    my $out_fh = IO::File->new(">$outfile");
    unless ($out_fh) {
        $self->error_message("Failed to create output filehandle");
        die;
    }
    $out_fh->print("CHR\tSTART\tSTOP\tSAMPLE\tCN\tPVAL\tCALL\n");

    ## Open a summary output file ##
    my $out_summary_fh = IO::File->new(">$outfile.summary");
    unless ($out_summary_fh) {
        $self->error_message("Failed to create output filehandle");
        die;
    }
    $out_summary_fh->print("CHR\tSTART\tSTOP\tAMP\tDEL\tNEUTRAL\tNA\n");

    #open region_call_results_dir and parse filenames
    #filenames: chr_start_end_call.csv
    opendir REGIONS,$region_call_dir;
    while (my $filename = readdir REGIONS) {
        next if $filename =~ /^ROI/;
        my $full_path_filename = $region_call_dir . "/" . $filename;
        next if $outfile =~ /$filename/;
        (my $chr = $filename) =~ s/(\w+)_(\d+)_(\d+)_call\.csv/$1/;
        (my $start = $filename) =~ s/(\w+)_(\d+)_(\d+)_call\.csv/$2/;
        (my $stop = $filename) =~ s/(\w+)_(\d+)_(\d+)_call\.csv/$3/;

        my $num_amp = my $num_del = my $num_neutral = my $num_na = 0;

        #open file, read data
        my $in_fh = new IO::File $full_path_filename,"r";
        while (my $line = $in_fh->getline) {
            
            next if $line =~ /samples/;
            (my $sample, my $mean, my $cn, my $sd, my $pvalue, my $fdr) = split /,/,$line;
            $out_fh->print("$chr\t$start\t$stop\t$sample\t$cn\t$pvalue\t");

            ## Genotype based on CN thresholds ##

            if(defined($self->amp_threshold) && defined($self->del_threshold))
            {
                if($cn eq "NA")
                {
                    $out_fh->print("NA\n");
                    $num_na++;
                }
                elsif($cn >= $self->amp_threshold)
                {
                    $out_fh->print("Amp\n");
                    $num_amp++;
                }
                elsif($cn <= $self->del_threshold)
                {
                    $out_fh->print("Del\n");
                    $num_del++;
                }
                else
                {
                    $out_fh->print("Neutral\n");
                    $num_neutral++;
                }
            }
            
            ## Genotype based on p-value ##
            
            else
            {
                if ($pvalue > $pval_cutoff) {
                    $out_fh->print("Neutral\n");
                    $num_neutral++;
                    next;
                }
                if ($pvalue <= $pval_cutoff) {
                    my $sdnew = $sd - 1;
                    if ($sd < 0.0) {
                        $out_fh->print("Del\n");
                        $num_del++;
                        next;
                    }
                    if ($sd > 0.0) {
                        $out_fh->print("Amp\n");
                        $num_amp++;
                        next;
                    }
                    if ($sd == 0.0) {
                        $out_fh->print("Neutral\n");
                        $num_neutral++;
                    }
                    else {
                        die "Could not resolve sd value $sd for sample $sample, chr $chr, range $start - $stop.\n";
                    }
                }#end, if pvalue <= cutoff
                else {
                    die "Could not resolve pvalue $pvalue for sample $sample, chr $chr, range $start - $stop.\n";
                }                
            }

        }
        
        $in_fh->close;
        
        ## Calculate fraction of amp and del ##
        
        my $amp_freq = my $del_freq = "";
        
        if($num_amp || $num_del || $num_neutral)
        {
            my $num_called = $num_amp + $num_del + $num_neutral;
            $amp_freq = sprintf("%.2f", ($num_amp / $num_called * 100)) . '%';
            $del_freq = sprintf("%.2f", ($num_del / $num_called * 100)) . '%';
        }
        
        ## Print summary to outfile ##
        $out_summary_fh->print(join("\t", $chr, $start, $stop, $num_amp, $num_del, $num_neutral, $num_na, $amp_freq, $del_freq) . "\n");
    }
    $out_fh->close;
    return 1;
}
1;
