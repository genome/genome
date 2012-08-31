
package Genome::Model::Tools::Analysis::Illumina::ModelReadcountErrors;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelReadcountErrors - Align reads with SSAHA2 or other aligner
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Illumina::ModelReadcountErrors {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "List of variants to count for" },
		sample_readcounts	=> { is => 'Text', doc => "Varscan readcounts file(s) for sample of interest (comma-delim)" },
		control_readcounts	=> { is => 'Text', doc => "Varscan readcounts file(s) for control sample(s) (comma-delim)" },
		output_file		=> { is => 'Text', doc => "Output file for adjusted readcounts" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Model errors and adjust calls for deep-readcount experiments"                 
}

sub help_synopsis {
    return <<EOS
This command models errors and adjust calls for deep-readcount experiments
EXAMPLE:	gmt analysis illumina model-readcount-errors --sample-readcounts sample.readcounts --control-readcounts control1.readcounts,control2.readcounts --output-file sample-adjusted.readcounts
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
        my $variants_file = $self->variants_file;
	my $sample_files = $self->sample_readcounts;
        my $control_files = $self->control_readcounts;
	my $output_file = $self->output_file;

        
        ## Step 1: Load readcounts for control files ##
        my %control_coverage = my %control_errors = ();
        my %ref_bases = ();
        my @control_files = split(/\,/, $control_files);
        
        foreach my $control_file (@control_files)
        {
                # Load the readcounts ##
                my $input = new FileHandle ($control_file);
                my $lineCounter = 0;
                
                while (<$input>)
                {
        		chomp;
        		my $line = $_;
        		$lineCounter++;		
			
                        if($lineCounter > 1)
                        {
                                my ($chrom, $position, $ref, $depth, $q20_depth) = split(/\t/, $line);
                                my @lineContents = split(/\t/, $line);
                                my $numContents = @lineContents;
                                for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
                                {
                                        my ($base, $reads) = split(/\:/, $lineContents[$colCounter]);
                                        
                                        my $key = join("\t", $chrom, $position, $base);

                                        $ref_bases{$chrom . "\t" . $position} = $ref;

                                        $control_coverage{$key} = 0 if(!$control_coverage{$key});
                                        $control_coverage{$key} += $q20_depth;
                                        $control_errors{$key} += $reads;
                                }
                        }
                }
	
                close($input);
        }


        ## Step 2: Load readcounts for sample files ##
        my %sample_coverage = my %sample_errors = ();
        my @sample_files = split(/\,/, $sample_files);
        
        foreach my $sample_file (@sample_files)
        {
                # Load the readcounts ##
                my $input = new FileHandle ($sample_file);
                my $lineCounter = 0;
                
                while (<$input>)
                {
        		chomp;
        		my $line = $_;
        		$lineCounter++;		
			
                        if($lineCounter > 1)
                        {
                                my ($chrom, $position, $ref, $depth, $q20_depth) = split(/\t/, $line);
                                my @lineContents = split(/\t/, $line);
                                my $numContents = @lineContents;
                                for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
                                {
                                        my ($base, $reads) = split(/\:/, $lineContents[$colCounter]);
                                        
                                        my $key = join("\t", $chrom, $position, $base);

                                        $ref_bases{$chrom . "\t" . $position} = $ref;

                                        $sample_coverage{$key} = 0 if(!$sample_coverage{$key});
                                        $sample_coverage{$key} += $q20_depth;
                                        $sample_errors{$key} += $reads;
                                }
                        }
                }
	
                close($input);
        }


        ## Open the output file ##
        
        open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
        print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\terrA\terrC\terrG\terrT\tdepth\treads1\treads2\terrors2\tadj_r2\tadj_freq\talt1_base\talt1_freq\talt2_base\talt2_freq\n";


        ## Step 3: Parse the variants file ##
        
	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

                my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
                my $position = $chr_start;

                ## Print the variant ##
                print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                
                if($ref_bases{$chrom . "\t" . $position})
                {
                        my $ref_base = $ref_bases{$chrom . "\t" . $position};

                        ## Get ref reads1 count ##
                        
                        my $ref_key = join("\t", $chrom, $position, $ref);
                        my $reads1 = $sample_errors{$ref_key};
                        
                        ## Determine what the 3 possible variant bases are ##
                        my @bases = ("A", "C", "G", "T");
                        my @var_bases = ();
                        @var_bases = ("C", "G", "T") if($ref_base eq "A");
                        @var_bases = ("A", "G", "T") if($ref_base eq "C");
                        @var_bases = ("A", "C", "T") if($ref_base eq "G");
                        @var_bases = ("A", "C", "G") if($ref_base eq "T");
                        
                        
                        ############## PRINT THE ERROR RATE ##############
                        
                        my $error_string = "";
                        my %errors_by_base = ();
                        
                        ## Calculate error rate for each base ##
                        
                        foreach my $base (sort @bases)
                        {
                                my $key = join("\t", $chrom, $position, $base);
                                $error_string .= "\t" if($error_string);
                                
                                if($control_coverage{$key} && $sample_coverage{$key})
                                {
                                        ## Calculate error rate for this base using control samples ##
                                        my $base_error = $control_errors{$key} / $control_coverage{$key};
                                        $errors_by_base{$base} = $base_error;
#                                        $control_string .= "\t$base\t" . sprintf("%.3f", $base_error * 100) . "\%";
                                        if($ref && $base eq $ref)
                                        {
                                                $error_string .= "NA_ref";
                                        }
                                        else
                                        {
                                                $error_string .= sprintf("%.3f", $base_error * 100) . "\%";
                                        }
                                                
                                }
                                else
                                {
                                        $error_string .= "-";
                                }
                        }
                        
                        print OUTFILE "\t$error_string";


                        ############## PRINT THE VARIANT BASE ##############
                        
                        my $variant_string = "";
                        
                        my $key = join("\t", $chrom, $position, $var);

                        if($control_coverage{$key} && $sample_coverage{$key})
                        {
                                ## Get the error for this base ##
                                my $base_error = $errors_by_base{$var};
                                my $exp_sample_base_errors = sprintf("%d", ($sample_coverage{$key} * $base_error));

                                my $reads2 = $sample_errors{$key};
                                ## Calculate an adjusted reads2 ##
                                        
                                my $adjusted_reads2 = $sample_errors{$key} - $exp_sample_base_errors;
                                $adjusted_reads2 = 0 if($adjusted_reads2 < 0);
                                        
                                ## Calculate adjusted reads1 ##
                                                                                
                                my $adjusted_freq = $adjusted_reads2 / $sample_coverage{$key} * 100;
                                $adjusted_freq = sprintf("%.3f", $adjusted_freq);
                                
                                $variant_string .= "$sample_coverage{$key}\t$reads1\t$reads2\t$exp_sample_base_errors\t$adjusted_reads2\t$adjusted_freq\%";                                
                        	print OUTFILE "\t$variant_string";
                        }
                        else
                        {
                                print OUTFILE "\tNC";
                        }


                        ############## PRINT OTHER BASES ##############
                        
                        my $other_string = "";
                                
                        foreach my $base (sort @var_bases)
                        {
                                my $key = join("\t", $chrom, $position, $base);
                                
                                if($control_coverage{$key} && $sample_coverage{$key} && $base ne $var)
                                {
                                        ## Calculate expected number of errors for this base in sample based on coverage ##
                                        my $base_error = $errors_by_base{$base};
                                        my $exp_sample_base_errors = sprintf("%d", ($sample_coverage{$key} * $base_error));

                                        my $reads2 = $sample_errors{$key};
                                        ## Calculate an adjusted reads2 ##
                                        
                                        my $adjusted_reads2 = $sample_errors{$key} - $exp_sample_base_errors;
                                        $adjusted_reads2 = 0 if($adjusted_reads2 < 0);
                                        
                                        ## Calculate adjusted reads1 ##
                                                                                
                                        my $adjusted_freq = $adjusted_reads2 / $sample_coverage{$key} * 100;
                                        $adjusted_freq = sprintf("%.3f", $adjusted_freq);

                                        if($base ne $ref && $base ne $var)
                                        {
                                        	$other_string .= "\t" if($other_string);
#                                                $other_string .= "$base\:$reads2\:$exp_sample_base_errors\:$adjusted_reads2\:$adjusted_freq\%";
                                                $other_string .= "$base\t$adjusted_freq\%";
                                        }                               
                                }
                        }

                        print OUTFILE "\t$other_string";

                }
                else
                {
                        print OUTFILE "\tNC";
                }
                
                ## End the line ##
                print OUTFILE "\n";
	}
	
	close($input);        

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

