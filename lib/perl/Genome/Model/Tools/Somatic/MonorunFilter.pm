package Genome::Model::Tools::Somatic::MonorunFilter;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Readonly;
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::MonorunFilter {
    is => 'Command',
    has => [
       'variant_file' => {
           type => 'String',
           is_input => 1,
           doc => 'List of variant positions in annotation format',
       },
       'output_file' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
           doc => 'File name in which to write output',
       },
       'filtered_file' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
	   is_optional => 1,
           doc => 'File name in which to write variants that were filtered',
       },       
       'tumor_bam_file' => {
            type => 'String',
            doc => 'Tumor bam file in which to examine reads',
            is_input => 1,
       },
       'min_homopolymer' => {
            type => 'String',
            default => '4',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum number of bases in homopolymer',
       },
       'reference' => {
            type => 'String',
            example_values => ['/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa'],
            is_optional => 1,
            is_input => 1,
            doc => 'Reference sequence to use',
       },
       prepend_chr => {
           is => 'Boolean',
           default => '0',
           is_optional => 1,
           is_input => 1,
           doc => 'prepend the string "chr" to chromosome names. This is primarily used for external/imported bam files.',
       },
       verbose => {
           is => 'Boolean',
           default => '0',
           is_optional => 1,
           is_input => 1,
           doc => 'Print the filtering result for each site.',
       },
       # Make workflow choose 64 bit blades
       lsf_resource => {
            is_param => 1,
            default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1]',
       },
       lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
       },
       skip => {
           is => 'Boolean',
           default => '0',
           is_input => 1,
           is_optional => 1,
           doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
       },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
    ]
};

sub help_brief {
    return "This module removes predicted indels (in annotation format) that occur near homopolymer sequences";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic strand-filter --variant-file somatic.snvs --tumor-bam tumor.bam --output-file somatic.snvs.strandfilter 
EOS
}

sub help_detail {                           
    return <<EOS 
This module removes predicted indels (in annotation format) that occur near homopolymer sequences
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }
    
    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    #test architecture to make sure we can run read count program
    unless (`uname -a` =~ /x86_64/) {
       $self->error_message("Must run on a 64 bit machine");
       die;
    }

    #check on BAM file
    unless(-e $self->tumor_bam_file) {
        $self->error_message("Tumor bam file: " . $self->tumor_bam_file . " does not exist");
        die;
    }

    unless(-e $self->tumor_bam_file . ".bai") {
        $self->error_message("Tumor bam must be indexed");
        die;
    }


    ## Determine the strandedness and read position thresholds ##
    
    my $min_homopolymer = $self->min_homopolymer;
    my $reference = $self->reference;

    ## Reset counters ##
    
    my %stats = ();
    $stats{'num_variants'} = $stats{'num_pass_filter'} = 0;

    ## Open the output file ##
    
    my $ofh = IO::File->new($self->output_file, "w");
    unless($ofh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing. $!");
        die;
    }


    ## Open the variants file ##

    my $input = new FileHandle ($self->variant_file);

    unless($input) {
        $self->error_message("Unable to open " . $self->variant_file . ". $!");
        die;
    }

    ## Reopen file for parsing ##

    $input = new FileHandle ($self->variant_file);
    

    ## Parse the variants file ##

    my $lineCounter = 0;
    
    while (<$input>)
    {
            chomp;
            my $line = $_;
            $lineCounter++;

            $stats{'num_variants'}++;
            
#            if($lineCounter <= 10)
 #           {
                (my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
                
		$ref = uc($ref);
		$var = uc($var);
		
		my $indel_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		my $indel_type = my $indel_size = "";
		
		if($ref eq '-')
		{
		    $indel_type = "INSERTION";
		    $indel_size = length($var);
		}
		elsif($var eq '-')
		{
		    $indel_type = "DELETION";
		    $indel_size = length($ref);
		}
		else
		{
		    $indel_type = "SNV";
		    $indel_size = 1;
		}

		## Set indels to pass filter by default ##
		my $filter_status = "Pass";

		if($indel_size && $indel_size > 0) # <= 2
		{
		    ## Build homopolymer sequence ##
		    
		    my $homoA = 'A' x $min_homopolymer;
		    my $homoC = 'C' x $min_homopolymer;
		    my $homoG = 'G' x $min_homopolymer;
		    my $homoT = 'T' x $min_homopolymer;
		    
		    ## Get the flanking sequence ##
		    
		    my $query_string = "";
		    $chr_start = $chr_start - 5;
		    $chr_stop = $chr_stop + 5;

		    if($self->prepend_chr)
		    {
			$query_string = "chr" . $chrom . ":" . $chr_start . "-" . $chr_stop;
		    }
		    else
		    {
			$query_string = $chrom . ":" . $chr_start . "-" . $chr_stop;
		    }

		    my $sequence = "";

		    if($indel_size <= 2)
		    {
			$sequence = `samtools faidx $reference $query_string | grep -v \">\"`;
			chomp($sequence);
			
			if($sequence)
			{
			    if($sequence =~ $homoA || $sequence =~ $homoC || $sequence =~ $homoG || $sequence =~ $homoT)
			    {
				$filter_status = "Fail";
			    }
    
			}
		    }
		    
		    if($self->verbose)
		    {
        		    print "$indel_key\t$indel_type\t$indel_size\t$sequence\t$filter_status\n";		
		    }
		}

    		if($filter_status eq "Pass")
		{
		    print $ofh "$line\n";
		    $stats{'num_pass_filter'}++;
		}

	    
    }
    
    close($input);

    print $stats{'num_variants'} . " variants\n";
    print $stats{'num_pass_filter'} . " passed the strand filter\n";

    return 1;
}



1;
