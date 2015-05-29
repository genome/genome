
package Genome::Model::Tools::Analysis::Sammy::CompileAll;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CompileAll - Call somatic variants from normal/tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	07/28/2009 by D.K.
#	MODIFIED:	07/28/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Sammy::CompileAll {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		output_dir	=> { is => 'Text', doc => "Output directory for somatic calls", is_optional => 0 },
		sample_name	=> { is => 'Text', doc => "Sample name for file naming purposes", is_optional => 0 },
		somatic_file	=> { is => 'Text', doc => "File name of Sammy-somatic calls", is_optional => 1 },
		p_value	        => { is => 'Text', doc => "P-value threshold for somatic variants [1.0E-06]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compiles somatic variants from normal and tumor BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command compiles somatic variants from Sammy somatic calls
EXAMPLE:	gmt analysis sammy
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

	## Create directory ##
	mkdir($self->output_dir) if(!(-d $self->output_dir));
	my $cmd;
	
	my $p_value = 1.0E-06;
	$p_value = $self->p_value if($self->p_value);

	## Determine dbSNP file ##
	
	my $dbsnp_file = "";
	
	if($self->dbsnp_file)
	{
		$dbsnp_file = $self->dbsnp_file;
	}
	else
	{
		$dbsnp_file = "/gscuser/dkoboldt/SNPseek/SNPseek2/ucsc/snp130.txt.nochrPositions";
	}
	
	## Determine somatic status file ##

	my $somatic_status_file = "";
	if($self->somatic_file)
	{
		$somatic_status_file = $self->somatic_file;
	}
	else
	{
		## infer the name of the file ##
		$somatic_status_file = $self->output_dir . "/" . $self->sample_name . ".snps.compared.status";
	}

	my $somatic_snp_file = $self->output_dir . "/" . $self->sample_name . ".all";
	my $num_somatic = parse_somatic_variants($self, $somatic_status_file, $somatic_snp_file);
	print "$num_somatic Somatic variants\n";

	if($num_somatic > 0)
	{	
		my $formatted_file = $somatic_snp_file . ".formatted";
		my $annotated_file = $formatted_file . ".annotations";
		my $merged_file = $somatic_snp_file . ".annotated";
		
		## Format SNPs for annotation ##

		system("perl ~dkoboldt/src/mptrunk/trunk/Auto454/format_snps_for_annotation.pl $somatic_snp_file $formatted_file");

		my $num_non_dbsnp = `cat $formatted_file | wc -l`;
		chomp($num_non_dbsnp);
		
		print "$num_non_dbsnp Somatic non-dbSNP variants formatted for annotation\n";
		print "Running annotation...\n";
		
        my $variants_obj = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            variant_file => $formatted_file,
            output_file => $annotated_file,
        );
        $variants_obj->execute;

        print "Merging annotations with SNP calls...\n";

        system("perl ~dkoboldt/src/mptrunk/trunk/Auto454/format_snps_with_annotation.pl $somatic_snp_file $annotated_file $merged_file");

    }


    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





################################################################################################
# Execute - the main program logic
#
################################################################################################

sub compile_category
{

}



###############################################################################################
# parse_somatic_variants - parse out the Somatic variants from the status file
#
###############################################################################################

sub parse_somatic_variants
{
    (my $self, my $infile, my $outfile) = @_;

    my $p_value_threshold = 1.0E-06;
    $p_value_threshold = $self->p_value if($self->p_value);

    print "Infile: $infile\nOutfile: $outfile P-value $p_value_threshold\n";

    my $num_somatic = 0;

    open(OUTFILE, ">$outfile") or die "Can't open outfile $outfile: $!\n";

    my $input = new FileHandle ($infile);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;		

        if($line)
        {
            my @lineContents = split(/\t/, $line);
            my $chrom = $lineContents[0];
            my $position = $lineContents[1];
            my $allele1 = $lineContents[2];
            my $allele2 = $lineContents[3];
            my $status = $lineContents[12];
            my $p_value = $lineContents[13];

            if($status && $status ne "Reference" && $p_value <= $p_value_threshold)
            {
                print OUTFILE "$line\n";
                print "$chrom\t$position\t$allele1\t$allele2\t$status\t$p_value\n";

                $num_somatic++;
            }

        }

#		return(0) if($lineCounter > 40);
    }


    close($input);

    close(OUTFILE);

    return($num_somatic);
}



###############################################################################################
# parse_somatic_variants - parse out the Somatic variants from the status file
#
###############################################################################################

sub parse_loh_variants
{
    (my $self, my $infile, my $outfile) = @_;

    my $p_value_threshold = 1.0E-06;
    $p_value_threshold = $self->p_value if($self->p_value);

    print "Infile: $infile\nOutfile: $outfile P-value $p_value_threshold\n";

    my $num_somatic = 0;

    open(OUTFILE, ">$outfile") or die "Can't open outfile $outfile: $!\n";

    my $input = new FileHandle ($infile);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;		

        if($line)
        {
            my @lineContents = split(/\t/, $line);
            my $chrom = $lineContents[0];
            my $position = $lineContents[1];
            my $allele1 = $lineContents[2];
            my $allele2 = $lineContents[3];
            my $status = $lineContents[12];
            my $p_value = $lineContents[13];

            if($status && $status eq "LOH" && $p_value <= $p_value_threshold)
            {
                print OUTFILE "$line\n";
                print "$chrom\t$position\t$allele1\t$allele2\t$status\t$p_value\n";

                $num_somatic++;
            }

        }

#		return(0) if($lineCounter > 40);
    }


    close($input);

    close(OUTFILE);

    return($num_somatic);
}


sub call_sammy
{
    my $classpath = "/gscuser/dkoboldt/Software/Sammy2";
    my $cmd = "java -Xms2000m -Xmx2000m -classpath $classpath Sammy ";
    return($cmd);
}



sub commify
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

1;

