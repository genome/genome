package Genome::Model::phaseSnvs;     # rename this when you give the module file a different name <--

use strict;
use warnings;
use File::Basename;
use File::Spec; 
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::phaseSnvs {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
       distance      	=> { is => 'String', is_optional => 1, doc => "Max Distance between SNV location to check if in phase" },
         command    	=> { is => 'String', is_optional => 1, doc => "Command user wants to run - B for build or C for check"},
        vcfFile   => { is => 'String', is_optional => 1, doc => "Path to VCF file to find SNVs"},
        bamFile     => { is => 'String', is_optional => 1, doc => "Path to bam file"},
        chromosome => { is => 'String', is_optional => 1, doc => "Chromosome of interest"},  
	snvs => { is => 'String', is_optional => 1, doc => "location and alt alleles"},
	relax => { is => 'String', is_optional => 1, doc => "use intsec failed SNVs"},
	sample => { is => 'String', is_optional => 1, doc => "sample name"}, 
    ], 
};

sub new {
    my ($class) = @_;
    my $self = bless {}, $class;
    return $self; 
}


sub help_brief {
     "A tool to find SNVs in phase",
}

sub help_detail {
	
}

sub execute {
	my $self=shift;
	my $dirName = dirname(__FILE__); 
	
	my $command = "java -jar ";
	$command .= File::Spec->catfile( $dirName, 'PhaseSnvs.jar' );

	$command .= ' --m ';
	$command .= $self->distance;

	$command .= ' --c ';
	$command .= $self->command;
	
	$command .= ' --f ';
	$command .= $self->vcfFile;

	$command .= ' --b ';
	$command .= $self->bamFile;
		
	$command .= ' --chr ';
	$command .= $self->chromosome;	
	
	$command .= ' --s ';
	$command .= $self->sample;

	$command .= ' --R ';
	$command .= $self->relax;

	$command .= ' --l ';
	$command .= $self->snvs;

	# Figure out where the accessory files are for R (these are not provided at the command line by the user)
	#my $javaLib = File::Spec->catfile( $dirName, 'OncoPlot', 'Sink.jar' );
	#my $rHelpLib = File::Spec->catfile($dirName, 'OncoPlot', 'oncoPrintHelper.R');
	
	#$command .= ' --helper-file=';
	#$command .= $rHelpLib;

	#$command .= ' --java-lib-file=';
	#$command .= $javaLib;
	
	print($command);
	print("\n");
	system($command); # this executes the command
	1;
 }
