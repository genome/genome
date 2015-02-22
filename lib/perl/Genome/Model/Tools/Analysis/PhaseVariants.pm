package Genome::Model::Tools::Analysis::PhaseVariants;

use strict;
use warnings;
use File::Basename;
use File::Spec;
use Genome;

class Genome::Model::Tools::Analysis::PhaseVariants {
    is => 'Command',
    has => [
        distance    => { is => 'String', is_optional => 1, doc => "Max Distance between SNV location to check if in phase" },
        command     => { is => 'String', is_optional => 1, doc => "Command user wants to run - B for build or C for check"},
        vcf_file    => { is => 'String', is_optional => 1, doc => "Path to VCF file to find SNVs"},
        bam_file    => { is => 'String', is_optional => 1, doc => "Path to bam file"},
        chromosome  => { is => 'String', is_optional => 1, doc => "Chromosome of interest"},
        snvs        => { is => 'String', is_optional => 1, doc => "location and alt alleles"},
        relax       => { is => 'String', is_optional => 1, doc => "use intsec failed SNVs"},
        sample      => { is => 'String', is_optional => 1, doc => "sample name"},
        print_reads => { is => 'String', is_optional => 1, doc => " print reads?"},
        out_file    => { is => 'String', is_optional => 1, doc => " output file"},
    ],
};

sub help_brief {
     "A tool to find SNVs in phase",
}

sub help_detail {

}

sub execute {
    my $self=shift;
    my $dirName = dirname(__FILE__);

    my $command = "java -jar ";
    $command .= File::Spec->catfile( $dirName, 'PhaseVariants.jar' );

    $command .= ' --o ';
    $command .= $self->out_file;

    $command .= ' --p ';
    $command .= $self->print_reads;

    $command .= ' --m ';
    $command .= $self->distance;

    $command .= ' --c ';
    $command .= $self->command;

    $command .= ' --f ';
    $command .= $self->vcf_file;

    $command .= ' --b ';
    $command .= $self->bam_file;

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

    $self->debug_message($command);
    #TODO replace with Genome::Sys->shellcmd(cmd => $command);
    system($command);

    1;
}
