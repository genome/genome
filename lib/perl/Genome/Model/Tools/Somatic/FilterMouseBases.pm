package Genome::Model::Tools::Somatic::FilterMouseBases;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Readonly;
use Genome::Info::IUB;
use Genome::Model::Tools::Capture::Helpers 'iupac_to_base';

class Genome::Model::Tools::Somatic::FilterMouseBases {
    is => 'Command',
    has => [
       'variant_file' => {
           type => 'String',
           is_input => 1,
           doc => 'List of variant positions in annotation format',
       },
       'chain_file' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
           doc => 'UCSC liftover chain file for human to mouse',
           example_values => ['/gscuser/dkoboldt/SNPseek/SNPseek2/ucsc/liftOver/hg19ToMm9.over.chain'],
       },
       'human_reference' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
           doc => 'Indexed FASTA of the human genome (defaults to build 37)',
           example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'],
       },
       'mouse_reference' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
           doc => 'Indexed FASTA of the mouse genome ',
           example_values => ['/gscmnt/839/info/medseq/reference_sequences/NCBI-mouse-build37/all_sequences.fa'],
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
	   is_optional => 0,
           doc => 'File name in which to write variants that were filtered',
       },       
       permissive=> {
           is => 'Boolean',
           default => '0',
           is_optional => 1,
           is_input => 1,
           doc => 'If set to 1, only exclude a site if mouse base matches variant base exactly',
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
            default_value => 'long',
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
    return "This module removes predicted SNVs that match the homologous base in mouse";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic filter-mouse-bases --variant-file my.snvs.bed --output-file my.snvs.filtered.bed
EOS
}

sub help_detail {                           
    return <<EOS 
This module removes predicted SNVs that match the homologous base in mouse
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


    ## Reset counters ##
    
    my %stats = ();
    $stats{'num_variants'} = $stats{'num_variants_pass_filter'} = $stats{'num_variants_fail_filter'} = 0;

    ## Open the variants file ##

    my $input = new FileHandle ($self->variant_file);

    unless($input) {
        $self->error_message("Unable to open " . $self->variant_file . ". $!");
        die;
    }
    
    my %mouse_bases = get_mouse_bases($self);

    ## Open pass filter file ##
    
    open(OUTPASS, ">" . $self->output_file) or die "Can't open output file: $!\n";
    open(OUTFAIL, ">" . $self->filtered_file) or die "Can't open filtered file: $!\n";

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
            
	    (my $chrom, my $chr_start, my $chr_stop, my $alleles) = split(/\t/, $line);
	    my ($ref, $cns) = split(/\//, $alleles);
            my $var = iupac_to_base($ref, $cns);

	    my $key = join(":", $chrom, $chr_start, $chr_stop, $ref, $var);
	    
	    if($mouse_bases{$key})
	    {
		my ($mouse_result, $human_ref, $human_var, $mouse_base) = split(/\t/, $mouse_bases{$key});
		
		if($mouse_result eq 'MouseBaseMatchesRef')
		{
		    print OUTPASS join("\t", $line, $mouse_result, $human_ref, $human_var, $mouse_base) . "\n";
		    $stats{'num_variants_pass_filter'}++;
		}
		elsif($self->permissive && $mouse_result ne "MouseBaseMatchesVar")
		{
		    ## IF we're permissive and can't rule out a mouse-derived variant, print as pass ##
		    print OUTPASS join("\t", $line, $mouse_result, $human_ref, $human_var, $mouse_base) . "\n";
		    $stats{'num_variants_pass_filter'}++;
		}
		else
		{
		    print OUTFAIL join("\t", $line, $mouse_result, $human_ref, $human_var, $mouse_base) . "\n";
		    $stats{'num_variants_fail_filter'}++;
		    $stats{'num_variants_fail_filter_' . $mouse_result}++;
		}
		
		print join("\t", $line, $mouse_bases{$key}) . "\n";
	    }
	    else
	    {
		print OUTPASS join("\t", $line, "NoHomologousBase") . "\n";
		$stats{'num_variants_no_homologous_base'}++;
	    }
	    
    }
    
    close($input);

#    print $stats{'num_variants'} . " variants\n";
#   print $stats{'num_variants_pass_filter'} . " passed the strand filter\n";

    foreach my $key (sort keys %stats)
    {
	print "$stats{$key} $key\n";
    }

    return 1;
}



sub get_mouse_bases
{
    my ($self) = @_;
    
    my %mouse_result = ();
    ## Get temp output file ##
    
        ## Open the output file ##

    my $temp_output_file = Genome::Sys->create_temp_file_path;
    my $ofh = Genome::Sys->open_file_for_writing($temp_output_file);
    unless($ofh) {
        $self->error_message("Unable to open temp file for writing.");
        die;
    }
    
        ## Reopen file for parsing ##

    my $input = new FileHandle ($self->variant_file);
    
    ## Print the variants in UCSC bed format ##

    my $lineCounter = 0;
    
    while (<$input>)
    {
            chomp;
            my $line = $_;
            $lineCounter++;
            
    	    (my $chrom, my $chr_start, my $chr_stop, my $alleles) = split(/\t/, $line);
	    my ($ref, $cns) = split(/\//, $alleles);
	    my $var = iupac_to_base($ref, $cns);
            
	    $var = iupac_to_base($ref, $var);
	    my $id_string = join(":", $chrom, $chr_start, $chr_stop, $ref, $var);
	    
	    $chr_start-- if($chr_start == $chr_stop);	## Adjust for bed ##
	    
	    print $ofh join("\t", "chr" . $chrom, $chr_start, $chr_stop, $id_string) . "\n";
	    
    }
    
    close($input);
    
    close($ofh);
    
    ## Run Liftover ##
    print "Running liftOver...\n";
    print "liftOver $temp_output_file " . $self->chain_file . " $temp_output_file.mm9.mapped $temp_output_file.mm9.unmapped\n";
    system("liftOver $temp_output_file " . $self->chain_file . " $temp_output_file.mm9.mapped $temp_output_file.mm9.unmapped");
    print "Done\n";

        ## Reopen file for parsing ##
    

    my $mapped = new FileHandle ("$temp_output_file.mm9.mapped");
    
    ## Parse the mapped variants ##    
    while (<$mapped>)
    {
            chomp;
            my $line = $_;
            $lineCounter++;
            
	    (my $mouse_chrom, my $mouse_chr_start, my $mouse_chr_stop, my $name) = split(/\t/, $line);
            my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\:/, $name);

	    $mouse_chrom =~ s/chr//;
	    
	    		my $mouse_position = $mouse_chr_stop;
		my $position = $chr_stop;

		## Get the human sequence, += 5 bp ##
	
		my $flank_start = $position -5;
		my $flank_stop = $position + 5;
		my $coordinates = $chrom . ":" . $flank_start . "-" . $flank_stop;
		my $human_seq = get_sequence($self->human_reference, $coordinates);

		## Get the same for mouse ##
		$flank_start = $mouse_position -5;
		$flank_stop = $mouse_position + 5;
		$coordinates = $mouse_chrom . ":" . $flank_start . "-" . $flank_stop;
		my $mouse_seq = get_sequence($self->mouse_reference, $coordinates);		


		## Check to see if mouse sequence is in + or - orientation ##

		my $mouse_seq_rev = flip_seq($mouse_seq);
		my $forward_match = match_bases($human_seq, $mouse_seq);
		my $reverse_match = match_bases($human_seq, $mouse_seq_rev);


		## get the reference base from human ##		

		$coordinates = $chrom . ":" . $position . "-" . $position;
		my $human_base = get_sequence($self->human_reference, $coordinates);		

		## get the reference base from mouse ##		

		$coordinates = $mouse_chrom . ":" . $mouse_position . "-" . $mouse_position;
		my $mouse_base = get_sequence($self->mouse_reference, $coordinates);		

		## If mouse is aligned in reverse, flip its base ##
		my $mouse_base_id = 0;
		my $mouse_strand = "";
		if($reverse_match > $forward_match)
		{
			$mouse_base = flip_base($mouse_base);
			$mouse_seq = $mouse_seq_rev;
			$mouse_base_id = $reverse_match;
			$mouse_strand = "-";
		}
		else
		{
			$mouse_base_id = $forward_match;
			$mouse_strand = "+";
		}

		print "$line\n$human_seq\n$mouse_seq ($forward_match)\n$mouse_seq_rev ($reverse_match)\nHuman: $human_base Mouse: $mouse_base\n" if($self->verbose);

		my $comparison_result = "Unknown";

		if($human_base eq $ref)
		{
			if($mouse_base eq $ref)
			{
				$comparison_result = "MouseBaseMatchesRef\t$ref\t$var\t$mouse_base";
			}
			elsif($mouse_base eq $var)
			{
				$comparison_result = "MouseBaseMatchesVar\t$ref\t$var\t$mouse_base";
			}
			else
			{
				$comparison_result = "MouseBaseMatchesNeither\t$ref\t$var\t$mouse_base";
			}
		}
		else
		{
			$comparison_result = "HumanBaseNotRef\t$human_base\t$var\t$mouse_base";
		}
		
		$mouse_result{$name} = $comparison_result;
	    	    
	    ## Get Mouse base 
	    
    }
    
    close($input);
    
    return(%mouse_result);
}



#############################################################
# get_sequence - retrieves the sequence using samtools faidx
#
#############################################################

sub get_sequence
{
	my ($ref_seq, $coordinates) = @_;
	
	my $seq = `samtools faidx $ref_seq $coordinates | grep -v \"\>\"`;
	chomp($seq);
	
	return($seq);
}


#############################################################
# get_sequence - retrieves the sequence using samtools faidx
#
#############################################################

sub match_bases
{
	my ($seq1, $seq2) = @_;
	
	my @bases1 = split(//, $seq1);
	my @bases2 = split(//, $seq2);
	my $num_bases = @bases1;
	
	my $num_matches = 0;
	for(my $baseCounter = 0; $baseCounter < $num_bases; $baseCounter++)
	{
		if($bases1[$baseCounter] && $bases2[$baseCounter] && $bases1[$baseCounter] eq $bases2[$baseCounter])
		{
			$num_matches++;
		}
	}
	
	return($num_matches);
}


#############################################################
# flip_seq: reverse-complement a DNA sequence
#
#############################################################

sub flip_seq
{
	my $sequence = shift(@_);
	# Create reversed sequence #
	my $sequence_rev = "";
	
	my @bases = split(//, $sequence);
	foreach my $base (@bases)
	{
		$sequence_rev = flip_base($base) . $sequence_rev;
	}

	return($sequence_rev);
}


#############################################################
# flip_base: complement a base
#
#############################################################

sub flip_base
{
	my $base = shift(@_);
	$base = uc($base);
	
	if($base eq 'A')
	{
		return("T");
	}
	elsif($base eq 'C')
	{
		return("G");	
	}
	elsif($base eq 'G')
	{
		return("C");
	}
	elsif($base eq 'T')
	{
		return("A");
	}
	elsif($base eq 'a')
	{
		return("t");
	}
	elsif($base eq 'c')
	{
		return("g");	
	}
	elsif($base eq 'g')
	{
		return("c");
	}
	elsif($base eq 't')
	{
		return("a");
	}
	else
	{
		return($base);
	}
}

1;
