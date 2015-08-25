package Genome::Model::Tools::CopyNumber::PlotPercentage;

############################################################################################
## AUTHOR: Yang Li
## DATE: July 18 2013
## EMAIL: yl5682@truman.edu
## NOTES: This tools is based on and adapted from PlotSegments written by Chris Miller
############################################################################################

use strict;
use warnings;
use FileHandle;
use File::Basename;
use Genome;
use IO::File;
use POSIX;

require Genome::Sys;

class Genome::Model::Tools::CopyNumber::PlotPercentage {
    is => 'Command',
    has => [
		input => {
			is => 'String',
			is_optional => 0,
			doc => 'The tool takes two types of input files:
I. A group of individual files (--list-input)
The user needs to provide a list of all individual files
e.g /path/to/analysis/BRCA_list.tsv
|/path/to/analysis/BRCA/BRC1.seg.cbs
|/path/to/analysis/BRCA/BRC2.seg.cbs
|/path/to/analysis/BRCA/BRC3.seg.cbs
|/path/to/analysis/BRCA/BRC4.seg.cbs

II. One individual file that contains multiple samples (--single-input)
The file should be in the following format:
[Sample] [Chrom] [Loc_Start] [Loc_End] [Num_Mark] [Mean]
ME001    1       12910       192984    1242       0.1231

Options:
Users can customize the columns to extract information. 
See help for more information.',
		},
		
		list_input => {
			is => 'Boolean',
			is_optional => 1,
			doc => 'Type of input: a list of paths to individual files',
		},
		
		single_input => {
			is => 'Boolean',
			is_optional => 1,
			doc => 'Type of input: a single file of all samples',
		},
		
		col_chrom => {
        	is => 'Int',
        	is_input => 1,
        	is_optional => 1,
        	doc => 'Column index of chromosome',
        },
        
        col_start => {
        	is => 'Int',
        	is_input => 1,
        	is_optional => 1,
        	doc => 'Column index of start location',
        },
        
        col_end => {
        	is => 'Int',
        	is_input => 1,
        	is_optional => 1,
        	doc => 'Column index of end location',
        },
        
        col_mean => {
        	is => 'Int',
        	is_input => 1,
        	is_optional => 1,
        	doc => 'Column index of mean value',
        },
        
		plot_title => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'The title of the plot',
        },
	
		y_label => {
        	is => 'String',
        	is_input => 1,
        	is_optional => 1,
        	doc => 'Y-axis label name',
        	default => 'Percentage',
        },
        
		window_size => {
            is => 'Int',
            is_optional => 1,
            is_input => 1,
            default => 10000,
            doc => 'The size of the window frame',
        },

		baseline => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	default => 0,
        	doc => 'The baseline of the plot',
        },
        
        gain_threshold => {
        	is => 'Float',
        	is_optional => 1,
        	is_input => 1,
        	default => 2.5,
        	doc => 'Gain threshold (copy number state)',
        },
        
        loss_threshold => {
        	is => 'Float',
        	is_optional => 1,
        	is_input => 1,
        	default => 1.5,
        	doc => 'Loss threshold (copy number state)',
        },
        
		genome_build => {
			is => 'Int',
			is_optional => 1,
			is_input => 1,
			default => 37,
			doc => 'The genome build',
		},
		
        sex => {
			is => 'String',
			is_optional => 1,
		    is_input => 1,
			doc => 'sex of the sample - male, female, or autosomes',
		    default => 'male',
		},
		
		plot_width => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The width of the plot',
        	default => 8,
        },
        
        plot_height => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The height of the plot',
        	default => 4,
        },
        
		xmax => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The maximum value of x-axis',
        },
        
        xmin => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The minimum value of x-axis',
        },
        
		ymax => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The maximum value of y-axis',
        	default => 1,
        },
        
        ymin => {
        	is => 'Int',
        	is_optional => 1,
        	is_input => 1,
        	doc => 'The minimum value of y-axis',
        	default => -1,
        },
        
        label_size => {
			is => 'Float',
			is_optional => 1,
		        is_input => 1,
			doc => 'Size of the labels on the plot',
			default => 0.6,
		},
		
		lowres_min => {
			is => 'Integer',
			is_optional => 1,
			doc => 'if lowres is enabled, segments longer than this many bp (and < lowres_max) will be scaled up to the lowres_max value for visibility',
			default => '100000'
		},

		lowres_max => {
			is => 'Integer',
			is_optional => 1,
			doc => 'if lowres is enabled, segments shorter than this many bp (and > lowres_min) will be scaled up to the lowres_max value for visibility',
			default => '5000000'
		},
        
        gain_color => {
			is => 'String',
			is_optional => 1,
		    is_input => 1,
			doc => 'Color of the gain region',
		    default => 'Red',
		},
		
		loss_color => {
			is => 'String',
			is_optional => 1,
		    is_input => 1,
			doc => 'Color of the loss region',
		    default => 'Blue',
		},
		
		output_tsv => {
			is => 'String',
			is_optional => 1,
			is_input => 1,
			doc => 'Raw data',
		},
		
		rcommands_file => {
        	is => 'String',
        	is_optional => 1,
        	is_output => 1,
        	doc => 'R commands file (e.g. for debugging purposes)',
        },

		output_pdf => {
        	is => 'String',
        	is_optional => 0,
        	is_output => 1,
        	doc => 'Plot output in PDF format',
        },
        
        log2_input => {
        	is => 'Boolean',
        	is_optional => 0,
        	is_output => 1,
        	doc => 'Log2 values as input',
        },
        
        log10_input => {
        	is => 'Boolean',
        	is_optional => 0,
        	is_output => 1,
        	doc => 'Log10 values as input',
        },
        
        abs_input => {
        	is => 'Boolean',
        	is_optional => 0,
        	is_output => 1,
        	doc => 'Absolute values as input',
        },
	]
};

sub help_brief {
    "Analyze the percentage of copy numbers in a group, and generate the plot."
}

sub help_detail {
    "Generate a plot of copy number alterations by percentage in a group.\n\nExample:\ngmt copy-number plot-percentage --input copynumber.cn.seg --output-pdf output.pdf --genome-build 37"
}


sub log_base {
    my ($base, $value) = @_;
    $value = 0.000001 if($value == 0);
    return log($value)/log($base);
}

#-------------------------------------
sub getEntrypointsFile{
    my ($sex, $genome_build) = @_;
    #set the appropriate entrypoints file so that we know the
    # chrs and lengths
    my $entrypoints_file = "";
    if($sex eq "male"){
        if($genome_build eq "36"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build36-20130113") . "/entrypoints.male";
        } elsif ($genome_build eq "37"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build37-20130113") . "/entrypoints.male";
        }elsif ($genome_build eq "mm9"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","mouse-mm9-20130113") . "/entrypoints.male";
        }
 
    } elsif ($sex eq "female"){
        if($genome_build eq "36"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build36-20130113") . "/entrypoints.female";
        } elsif ($genome_build eq "37"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build37-20130113") . "/entrypoints.female";
        } elsif ($genome_build eq "mm9"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","mouse-mm9-20130113") . "/entrypoints.female";
        }
    } elsif ($sex eq "autosomes"){
        if($genome_build eq "36"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build36-20130113") . "/entrypoints.autosomes";
        } elsif ($genome_build eq "37"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","human-build37-20130113") . "/entrypoints.autosomes";
        } elsif ($genome_build eq "mm9"){
            $entrypoints_file = Genome::Sys->dbpath("tgi-misc-annotation","mouse-mm9-20130113") . "/entrypoints.autosomes";
        }
    }

    if ($entrypoints_file eq "") {
        if(-e $genome_build){
            print STDERR "Using custom annotation build: $genome_build\n";
            $entrypoints_file = $genome_build;
        } else {
            die "Specify a valid genome build and sex. Only genome builds 36/37 and male/female are currently supported";
        }
    }

    return $entrypoints_file;
}

############################################################################################
## Main
sub execute {
    my $self = shift;
    my $input_file = $self->input;
    my $list_input = $self->list_input;
    my $single_input = $self->single_input;
	my $col_chrom = $self->col_chrom;
	my $col_start = $self->col_start;
	my $col_end = $self->col_end;
	my $col_mean = $self->col_mean;
    my $plot_title = $self->plot_title;
    my $y_label = $self->y_label;
    my $window_size = $self->window_size;
    my $label_size = $self->label_size;
    my $xmax = $self->xmax;
    my $xmin = $self->xmin;
    my $ymax = $self->ymax;
    my $ymin = $self->ymin;
    my $plot_height = $self->plot_height;
    my $plot_width = $self->plot_width;
    my $baseline = $self->baseline;
    my $gain_threshold = $self->gain_threshold;
    my $loss_threshold = $self->loss_threshold;
	my $lowres_min = $self->lowres_min;
	my $lowres_max = $self->lowres_max;
    my $genome_build = $self->genome_build;
   	my $sex = $self->sex;
	my $gain_color = $self->gain_color;
	my $loss_color = $self->loss_color;
    my $output_tsv_file = $self->output_tsv;
    my $rcommands_file = $self->rcommands_file;
    my $output_pdf_file = $self->output_pdf;
    my $log2_input = $self->log2_input;
    my $log10_input = $self->log10_input;
    my $abs_input = $self->abs_input;
	
	## Check for input type
	unless (($log2_input xor $log10_input) xor ($abs_input)){
		die "Please specify the input type. Types include: log2-input, log10-input, and abs-input\n";
	}
	
	## Open input file
	open (my $input, "<$input_file") or die "Cannot open input file $input_file. $!.\n";
	if (!defined($output_tsv_file)){
		my ($tfh,$tfile) = Genome::Sys->create_temp_file;
		unless($tfh) {
		    $self->error_message("Unable to create temporary file $!");
		    die;
		}
		$output_tsv_file=$tfile;
	}
	
	## Open output file
	open (my $output_tsv, ">$output_tsv_file") or die "Cannot open output file $output_tsv_file. $!.\n";
	unless( (defined($list_input)) xor (defined($single_input)) ){
        die $self->error_message("You must specify either the list-input param OR the single-input param, but not both.");
    }

	## Get the details of the input files
	my @files;
	if ($list_input){ # Read the list of input files
		while (my $line = <$input>){
			chomp $line;
			push (@files, $line);
		}
		$col_chrom = defined($col_chrom) ? $col_chrom : 0;
		$col_start = defined($col_start) ? $col_start : 1;
		$col_end = defined($col_end) ? $col_end : 2;
		$col_mean = defined($col_mean) ? $col_mean : 4;
	} elsif ($single_input) { # Read the list file
		push (@files, $input_file);
		$col_chrom = defined($col_chrom) ? $col_chrom : 1;
		$col_start = defined($col_start) ? $col_start : 2;
		$col_end = defined($col_end) ? $col_end : 3;
		$col_mean = defined($col_mean) ? $col_mean : 5;
	}
	
	# Read each individual file into a hash
	my %samples;
	my %input_files;
	foreach my $file (@files){
		open (my $fh, "<$file") or die "Cannot open input file $file. $!.\n";
		my @array;
		while (my $line = <$fh>){
			chomp $line;
			my @line_arr = split ('\t', $line);
			
			next if (!defined($line_arr[0]));
			
			if (defined($single_input)){
				$samples{$line_arr[0]} = "0";
			}
		
			my $rec = {};
			$rec->{chrom} = $line_arr[$col_chrom];
			$rec->{start} = $line_arr[$col_start];
			$rec->{end} = $line_arr[$col_end];
			$rec->{mean} = $line_arr[$col_mean];
		
			push (@array, $rec);
		}
		close ($fh);
		$input_files{basename($file)} = \@array;
	}
	
	my $entrypoints_file = getEntrypointsFile($sex,$genome_build);
	open (my $entrypoints, "<$entrypoints_file") or die "Cannot open entrypoints file. $entrypoints_file. $!\n";
	
	# Read the chromosome length chart
	my %chrmlen;
	while (my $line = <$entrypoints>){
		chomp $line;
		my @line = split('\t', $line);
		$chrmlen{$line[0]} = $line[1];
	}
	
	# Calculate window slice number
	my %window_counts;
	foreach my $key (keys %chrmlen){
		my $window_count = ceil($chrmlen{$key} / $window_size);
		$window_counts{$key} = $window_count;
		#print "Chrom: $key Slices: $window_counts{$key}\n";
	}

	# Initialize a score book
	my %scores;
	foreach my $chrmlenKey (keys %chrmlen){
		my (@gain_array, @loss_array);
		for (my $i = 0; $i<scalar($window_counts{$chrmlenKey}); $i++){
			push (@gain_array, 0);
			push (@loss_array, 0);
		}
		$scores{"gain_$chrmlenKey"} = \@gain_array;
		$scores{"loss_$chrmlenKey"} = \@loss_array;
	}

	# Print hint
	print "[Hint] Increase window size to reduce computing time.\n";
	
	# Compare and score
	foreach my $chrmlenKey (keys %chrmlen){ # Iterate through 23 pairs of chromosomes
		print "Working on chromosome: $chrmlenKey\n";
		my @array;
		for (my $i=0; $i<$window_counts{$chrmlenKey}; $i++){ # Iterate through each window on each chromosome
			my $start = $window_size * $i; # Start pos of a window:  "|"        |
			my $end = $window_size * ($i + 1); # End pos of a window: |        "|"

			foreach my $key ( keys %input_files ){ # Iterate through all files
				#print "File: $key Lines: ".scalar(@{$input_files{$key}})."\n";
				for (my $j = 0; $j < scalar(@{$input_files{$key}}); $j++){ # Run down the lines in each file
			
					my $thisStart = ${@{$input_files{$key}}[$j]}{start};
					my $thisEnd = ${@{$input_files{$key}}[$j]}{end};
					#print "Start: $thisStart End: $thisEnd\n";
				
					next if (!defined(${@{$input_files{$key}}[$j]}{chrom}));
					next if (${@{$input_files{$key}}[$j]}{chrom} ne $chrmlenKey); # Skip the line if it is not the matching chromosome
					next if ($thisStart >= $end); # Skip non-overlapping regions: Start after window (outside) |        | oxxxxxo
					next if ($thisEnd <= $start); # Skip non-overlapping regions: End before window (outside): oxxxxxo  |        |
					next if ((${@{$input_files{$key}}[$j]}{mean} < $gain_threshold) 
							&& (${@{$input_files{$key}}[$j]}{mean} > $loss_threshold)); # Skip non-interesting regions
				
					# Input base conversion
					if ($log2_input) {
						${@{$input_files{$key}}[$j]}{mean} = 2*(2**${@{$input_files{$key}}[$j]}{mean});
					} elsif ($log10_input){
						${@{$input_files{$key}}[$j]}{mean} = 2*(10**${@{$input_files{$key}}[$j]}{mean});
					} elsif ($abs_input){
						;
					}
					
					# |   oxxxxxxo   | (Absolutely inside)
					if (($thisStart >= $start) && ($thisEnd < $end)){
						my $score = ($thisEnd - $thisStart) / $window_size;
						my $window = floor($thisStart / $window_size);

						if (${@{$input_files{$key}}[$j]}{mean} < $loss_threshold) { # loss
							@{$scores{"loss_$chrmlenKey"}}[$window] -= $score; # Write to loss
						}
						if (${@{$input_files{$key}}[$j]}{mean} > $gain_threshold) {
							@{$scores{"gain_$chrmlenKey"}}[$window] += $score; # Write to gain
						}
						
					}
				
					# |   oxxxxxxxxxx|xxxo (Partially inside - head)
					if (( $thisStart >= $start ) && ( $thisEnd >= $end )){
						my $score = ($end - $thisStart) / $window_size;
						my $window = floor($thisStart / $window_size);

						if (${@{$input_files{$key}}[$j]}{mean} < $loss_threshold) { # loss
							@{$scores{"loss_$chrmlenKey"}}[$window] -= $score; # Write to loss
						}
						if (${@{$input_files{$key}}[$j]}{mean} > $gain_threshold) {
							@{$scores{"gain_$chrmlenKey"}}[$window] += $score; # Write to gain
						}
						if ($score > 1 || $score < -1){
							print $score."\n";
						}
						
					} 
				
					# oxxx|xxxxxxxxxo    | (Partially inside - tail)
					if ( $thisStart <= $start && $thisEnd <= $end ){
						my $score = ($thisEnd - $start) / $window_size;
						my $window = floor($thisStart / $window_size);
					
						if (${@{$input_files{$key}}[$j]}{mean} < $loss_threshold) { # loss
							@{$scores{"loss_$chrmlenKey"}}[$window] -= $score; # Write to loss
						}
						if (${@{$input_files{$key}}[$j]}{mean} > $gain_threshold) {
							@{$scores{"gain_$chrmlenKey"}}[$window] += $score; # Write to gain
						}
						
						if ($score > 1 || $score < -1){
							print $score."\n";
						}
					}
				
					# oxxx|xxxxxxxxxx|xxxo (Entire coverage)
					if ( ($thisStart <= $start) && ($thisEnd >= $end) ) {
						my $window = floor(($start) / $window_size);
						if (${@{$input_files{$key}}[$j]}{mean} < $loss_threshold) { # loss
							@{$scores{"loss_$chrmlenKey"}}[$window] -= 1; # Write 1 to loss
						}
						if (${@{$input_files{$key}}[$j]}{mean} > $gain_threshold) {
							@{$scores{"gain_$chrmlenKey"}}[$window] += 1; # Write 1 to gain
						}
					}
				}
			}
		}
	}

	# Write to output files
	foreach my $chrmlenKey (keys %chrmlen){
		for (my $i = 0; $i<scalar($window_counts{$chrmlenKey}); $i++){
	
			my $total_number;
			#calculate total number of samples
			if ($list_input){
				$total_number = scalar keys %input_files; # sample number
			} elsif ($single_input) {
				$total_number = (scalar keys %samples) - 1; # sample number; correct the first row (not a sample)
			}
			
			if (@{$scores{"gain_$chrmlenKey"}}[$i] > $total_number) {
				print @{$scores{"gain_$chrmlenKey"}}[$i]."\n";
			}

			my $percentage_gain = @{$scores{"gain_$chrmlenKey"}}[$i] / $total_number;
			my $percentage_loss = @{$scores{"loss_$chrmlenKey"}}[$i] / $total_number;
			 
			if ($percentage_gain != 0){
				#print $percentage_gain."\n";
				print $output_tsv "$chrmlenKey\t".($i*$window_size)."\t".(($i+1)*$window_size)."\t"."0\t".$percentage_gain."\n";
			}
		
			if ($percentage_loss != 0){
				#print $percentage_loss."\n";
				print $output_tsv "$chrmlenKey\t".($i*$window_size)."\t".(($i+1)*$window_size)."\t"."0\t".$percentage_loss."\n";
			}	
		}

	}
	close($input);
	close($output_tsv);

	############################################### R Section ############################################################
	#open the R file
    if (!defined($rcommands_file)){
		my ($tfh,$tfile) = Genome::Sys->create_temp_file;
		unless($tfh) {
		    $self->error_message("Unable to create temporary file $!");
		    die;
		}
		$rcommands_file=$tfile;
	}
    
    open(R_COMMANDS,">$rcommands_file") || die "can't open $rcommands_file for writing\n";

    #source the R file
    my $dir_name = dirname(__FILE__);
    print R_COMMANDS "source(\"" . $dir_name . "/PlotSegments.R\")\n";
    
    #set up pdf parameters
    my $docwidth = $plot_width;
    my $docheight = $plot_height;
    print R_COMMANDS "pdf(file=\"" . $output_pdf_file . "\",width=" .$docwidth .",height=" . $docheight . ")\n";

	#set up the plotting space
	print R_COMMANDS "par(xaxs=\"i\", xpd=FALSE, mfrow=c(1,1), oma=c(1,1,1,1), mar=c(3,3,1,1))\n";

    #draw the plots for each set of segments
        print R_COMMANDS "plotSegments(";

        #first the core stuff
        print R_COMMANDS "chr=\"ALL\"";

        print R_COMMANDS ", filename=\"" . $output_tsv_file . "\"";
        print R_COMMANDS ", entrypoints=\"" . $entrypoints_file . "\"";

        #then the optional parameters
        if(defined($ymax)){
            print R_COMMANDS ", ymax=" . $ymax;
        }

        if(defined($ymin)){
            print R_COMMANDS ", ymin=" . $ymin;
        }

        if(defined($xmax)){
            print R_COMMANDS ", xmax=" . $xmax;
        }

        if(defined($xmin)){
            print R_COMMANDS ", xmin=" . $xmin;
        }

		print R_COMMANDS ", showNorm=FALSE";
        print R_COMMANDS ", gainThresh=" . $baseline;
        print R_COMMANDS ", lossThresh=" . $baseline;

        print R_COMMANDS ", gainColor=\"" . $gain_color . "\"";
        print R_COMMANDS ", lossColor=\"" . $loss_color . "\"";

        if(defined($baseline)){
            print R_COMMANDS ", baseline=\"" . $baseline . "\"";
        }
        
        if (defined($lowres_min)){
            print R_COMMANDS ", lowResMin=" . $lowres_min;
        }

        if (defined($lowres_max)){
            print R_COMMANDS ", lowResMax=" . $lowres_max;
        }
        
        if (defined($label_size)){
        	print R_COMMANDS ", label_size=$label_size";
        }
        
        print R_COMMANDS ", multiplePlot=FALSE";
        if (defined($plot_title)){
        	print R_COMMANDS ", plotTitle=\"$plot_title\"";
        }
        print R_COMMANDS ", drawLabel=TRUE";
        print R_COMMANDS ", percentagePlot=TRUE";
        
        print R_COMMANDS ", ylabel=\"$y_label\"";
        print R_COMMANDS ")\n";
 

    #close the file out
    print R_COMMANDS "dev.off()\n";
    print R_COMMANDS "q()\n";
    close R_COMMANDS;

    #now run the R command
    my $cmd = "R --vanilla --slave \< $rcommands_file";
    my $return = Genome::Sys->shellcmd(
	cmd => "$cmd",
        );
    unless($return) {
	$self->error_message("Failed to execute: Returned $return");
	die $self->error_message;
    }

}
1;
