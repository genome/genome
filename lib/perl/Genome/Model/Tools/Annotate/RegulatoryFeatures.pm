package Genome::Model::Tools::Annotate::RegulatoryFeatures;

#####################################################################################################################################
# RegulatoryFeatures - Annotate a list of variants or positions with ENCODE or UCSC regulatory features
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/22/2012 by D.K.
#	MODIFIED:	02/22/2012 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use Genome;

use Genome::Model::Tools::RefCov::ROI::Bed;

my $low  = 20000;
my $high = 250000;
UR::Context->object_cache_size_lowwater($low);
UR::Context->object_cache_size_highwater($high);

class Genome::Model::Tools::Annotate::RegulatoryFeatures {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "A list of variants in annotation format", is_optional => 0, is_input => 1},
		regfeats_file	=> { is => 'Text', doc => "Path to ENCODE regulatory features file", is_optional => 0, is_input => 1, default => '/gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/encode/release-66/AnnotatedFeatures.gff.nochr.tsv'},
		output_file	=> { is => 'Text', doc => "Output file for annotation-results", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Annotates a list of variants with ENCODE regulatory features"                 
}

sub help_synopsis {
    return <<EOS
This command annotates a list of variants with ENCODE regulatory features
EXAMPLE:	gmt annotate regulatory-features --variant-file myvariants.annot --output-file myvariants.annot.regFeats
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
    
    if(!(-e $self->variant_file))
    {
        die "Error: Variant file not found\n";
    }

    if(!(-e $self->regfeats_file))
    {
        die "Error: Regfeats file not found\n";
    }    
    
    my %regfeats = load_regfeats($self->regfeats_file);
    
    open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
    
    my $input = new FileHandle ($self->variant_file);
    my $lineCounter = 0;
	
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;
        
        my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
        $chr_stop =~ s/[^0-9]//g;
        $chr_stop = $chr_start if(!$chr_stop);  ## Hack to allow to work on VCF-like lines ##
        
        my $mb_start = sprintf("%d", $chr_start / 1000000);
        my $mb_stop = sprintf("%d", $chr_stop / 1000000);  
        
        my %variant_features = ();
        
        for(my $mb_pos = $mb_start; $mb_pos <= $mb_stop; $mb_pos++)
        {
            my $region_key = join("\t", $chrom, $mb_pos);
            if($regfeats{$region_key})
            {
                my @regfeats = split(/\n/, $regfeats{$region_key});
                foreach my $regfeat (@regfeats)
                {
                    my ($region_start, $region_stop, $feature) = split(/\t/, $regfeat);
                    if($region_start <= $chr_stop && $region_stop >= $chr_start)
                    {
                        $variant_features{$feature}++;
                    }
                }
            }
        }
        
        ## Compile all variant features into a single result ##
        
        my $variant_feature_list = "";
        
        foreach my $feature (sort keys %variant_features)
        {
            $variant_feature_list .= "," if($variant_feature_list);
            $variant_feature_list .= $feature;
        }

        $variant_feature_list = "-" if(!$variant_feature_list);
        print OUTFILE "$line\t$variant_feature_list\n";
        

    }
    close($input);
    close(OUTFILE);
    return 1;
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub load_regfeats
{                               # replace with real execution logic.
    my $FileName = shift(@_);

    my %regfeats = ();

    my $input = new FileHandle ($FileName);
    my $lineCounter = 0;
	
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;		

        my ($chrom, $chr_start, $chr_stop, $feature) = split(/\t/, $line);
        
        ## determine region start and stop in terms of megabases ##
        
        my $mb_start = sprintf("%d", $chr_start / 1000000);
        my $mb_stop = sprintf("%d", $chr_stop / 1000000);
        
        for(my $mb_pos = $mb_start; $mb_pos <= $mb_stop; $mb_pos++)
        {
            my $region_key = join("\t", $chrom, $mb_pos);
            $regfeats{$region_key} .= "\n" if($regfeats{$region_key});
            $regfeats{$region_key} .= join("\t", $chr_start, $chr_stop, $feature);
        }
    }
	
    close($input);

    return(%regfeats);
}



1;
