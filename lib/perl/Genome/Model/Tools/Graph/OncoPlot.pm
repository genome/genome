package Genome::Model::Tools::Graph::OncoPlot;     # rename this when you give the module file a different name <--

use strict;
use warnings;
use File::Basename;
use File::Spec; 
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Graph::OncoPlot {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        genes      	=> { is => 'String', is_optional => 1, doc => "Comma separated list (no spaces) of genes to analyze" },
        input    	=> { is => 'String', is_optional => 1, doc => "Input file (no header) where the 1st column is sample name; 2nd column is copyCat file; 3rd column is variant file"},
        geneModel   => { is => 'String', is_optional => 1, doc => "Path to a BED file with gene names (4th column should contain gene name)"},
        outFile     => { is => 'String', is_optional => 1, doc => "Path where the output image ought to be placed"},
        variantIndex => { is => 'String', is_optional => 1, doc => "Column number in the variant file that stores the variant type (i.e. 6)"},
        geneNameIndex => { is => 'String', is_optional => 1, doc => "Column number in the variant file that stores the gene name (i.e. 22)"},  
    ], 
};

sub new {
    my ($class) = @_;
    my $self = bless {}, $class;
    return $self; 
}


sub help_brief {
    "A tool to produce a visual summary of gene copy-number and variant status inspired by cBioPortal's OncoPrint",
}

sub help_detail {
	
}

sub execute {
	my $self=shift;
	
	my $command = 'OncoPlot/OncoPrint.R';

	$command .= ' --genes=';
	$command .= $self->genes;

	$command .= ' --gene-name-index=';
	$command .= $self->geneNameIndex;
	
	$command .= ' --variant-index=';
	$command .= $self->variantIndex;


	$command .= ' --input=';
	$command .= $self->input;
		
	$command .= ' --gene-model=';
	$command .= $self->geneModel;	
	
	$command .= ' --out-file=';
	$command .= $self->outFile;
	
	# Figure out where the accessory files are for R
	my $dirName = dirname(__FILE__); 
	my $javaLib = File::Spec->catfile( $dirName, 'OncoPlot', 'Sink.jar' );
	my $rHelpLib = File::Spec->catfile($dirName, 'OncoPlot', 'oncoPrintHelper.R');
	
	$command .= ' --helper-file=';
	$command .= $rHelpLib;

	$command .= ' --java-lib-file=';
	$command .= $javaLib;
	
#	print($command);
	system($command);
	1;
 }

sub execute2 {
	my $self=shift;
	my ($name) = @_;
	
	print("$name is running \n");
 }

