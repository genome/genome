package Genome::Model::Tools::DetectVariants2::Strelka;

use warnings;
use strict;

use Genome;
use Workflow;

#Basic running of Strelka
#cd /gscmnt/gc2142/techd/analysis/strelka/all_chrs/
#cp /gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/strelka/etc/strelka_config_bwa_default.ini .
#/gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/configureStrelkaWorkflow.pl --tumor /gscmnt/gc12001/info/model_data/2886829900/build125238592/alignments/125442070.bam --normal /gscmnt/gc8001/info/model_data/2883343119/build121512916/alignments/121537788.bam --ref /gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa --config strelka_config_bwa_default.ini --output-dir ./results
#time make -j 8 -C /gscmnt/gc2142/techd/analysis/strelka/all_chrs/results
my $DEFAULT_VERSION = '0.4.6.2';

class Genome::Model::Tools::DetectVariants2::Strelka {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces a list of high confidence somatic snps and indels.",
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>100000] span[hosts=1]',
        },
    ],
};

my %STRELKA_VERSIONS = (
    '0.4.6.2' => '/gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/',
);

sub help_synopsis {
    my $self = shift;
    return <<"EOS"

    gmt detect-variant2 strelka something useful...
EOS
}

sub help_detail {                           
    return <<EOS 
    Provide a tumor and normal BAM file and get a list of somatic SNVs and Indels.  
EOS
}

sub _detect_variants {
    my $self = shift;

    $self->status_message("beginning execute");

    my $output_dir = $self->_temp_staging_directory;

    #TODO:
    #Parse the params string specified when the strelka method is called (this will come from the snv_detection_strategy in the processing profile) or by setting
    #In the processing profile, the snv_detection_strategy will look something like this: strelka 0.4.6.2 [params string]
    #Instead actually passing this string through to a strelka command-line command, parse this string, produce a custom config file, and copy that to the working dir instead of the default one

    #It seems that the most compact Strelka config file looks something like this (with bwa default 0.4.6.2 values):
    #[user]
    #isSkipDepthFilters = 0
    #depthFilterMultiple = 3.0
    #snvMaxFilteredBasecallFrac = 0.4
    #snvMaxSpanningDeletionFrac = 0.75
    #indelMaxRefRepeat = 8
    #indelMaxWindowFilteredBasecallFrac = 0.3
    #indelMaxIntHpolLength = 14
    #ssnvPrior = 0.000001
    #sindelPrior = 0.000001
    #ssnvNoise = 0.0000005
    #sindelNoise = 0.000001
    #ssnvNoiseStrandBiasFrac = 0.5
    #minTier1Mapq = 20
    #ssnvQuality_LowerBound = 15
    #sindelQuality_LowerBound = 30
    #isWriteRealignedBam = 0
    #binSize = 25000000
    #extraStrelkaArguments =

    #The order of the options above does not seem to matter
    #Note that all additional arguments to Strelka must be passed to Strelka as a single string using the 'extraStrelkaArguments =' line
    #For example, you could do:
    #extraStrelkaArguments = -used-allele-count-min-qscore 30 -min-qscore 10
    #For a full list of Strelka options, run the strelka binary to view the help docs. 
    #/gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/strelka/bin/strelka
    #Note that there are many possible additional arguments!

    my $params_string = $self->params;
    #If the params are defined, use these to create a Strelka config otherwise use a default string
    if ($params_string){
      #Perform basic sanity checks of the params string
      #TODO: Sanity checking of params may need to be Strelka version specific

    }else{
      #TODO: The default params string may need to be Strelka version specific (implement something similar to the strelka_path)
      #Or perhaps we should not allow it to be left empty...
      $params_string = "isSkipDepthFilters = 0;depthFilterMultiple = 3.0;snvMaxFilteredBasecallFrac = 0.4;snvMaxSpanningDeletionFrac = 0.75;indelMaxRefRepeat = 8;indelMaxWindowFilteredBasecallFrac = 0.3;indelMaxIntHpolLength = 14;ssnvPrior = 0.000001;sindelPrior = 0.000001;ssnvNoise = 0.0000005;sindelNoise = 0.000001;ssnvNoiseStrandBiasFrac = 0.5;minTier1Mapq = 20;ssnvQuality_LowerBound = 15;sindelQuality_LowerBound = 30;isWriteRealignedBam = 0;binSize = 25000000;extraStrelkaArguments =";
    }
    #Prepend the '[user]' field at the beginning of the params string
    $params_string = "[user];" . $params_string;

    my $strelka_working_config_file = $output_dir . "/strelka_config.ini";
    my @strelka_params = split(";", $params_string);

    my $strelka_working_config = IO::File->new(">$strelka_working_config_file");
    $strelka_working_config->print(join("\n", @strelka_params));

    #Run the strelka configuration step that checks your input files and prepares a Makefile
    my $cmd = $self->strelka_path . "/strelka_workflow/configureStrelkaWorkflow.pl" 
                                   . " --tumor " . $self->aligned_reads_input 
                                   . " --normal " . $self->control_aligned_reads_input 
                                   . " --ref " . $self->reference_sequence_input 
                                   . " --config $strelka_working_config_file"
                                   . " --output-dir $output_dir/output";
    Genome::Sys->shellcmd( cmd=>$cmd,
                           input_files=>[$self->aligned_reads_input, $self->control_aligned_reads_input], 
                           output_files=>[$output_dir . "/output/Makefile"], 
                           skip_if_output_is_present=>1, 
                           allow_zero_size_output_files => 0, );

    #Actually make the Makefile (i.e. run the Strelka workflow)
    Genome::Sys->shellcmd( cmd=>"make -j 4 -C $output_dir/output/",
                           input_files=>[$output_dir . "/output/Makefile"], 
                           output_files=>[$output_dir . "/output/results/all.somatic.snvs.vcf", $output_dir . "/output/results/all.somatic.indels.vcf", $output_dir . "/output/results/passed.somatic.snvs.vcf", $output_dir . "/output/results/passed.somatic.indels.vcf"], 
                           skip_if_output_is_present=>1, 
                           allow_zero_size_output_files => 1, );

    #Create symlinks to the unfiltered SNV and Indel results generated by Strelka
    Genome::Sys->create_symlink($output_dir . "/output/results/all.somatic.snvs.vcf", $self->_snv_staging_output);
    Genome::Sys->create_symlink($output_dir . "/output/results/all.somatic.indels.vcf", $self->_indel_staging_output);

    #How is filtering implemented??
    #Scott overrode _run_bed_converter from DetectVariants2::Base below
    #It runs the runs the BED converter twice passing in a parameter that grabs either the passing (hq) or failing (lq) variants from Strelka itself

    $self->status_message("ending execute");
    return 1; 
}


sub _run_bed_converter {
    my $self = shift;
    my $converter = shift;
    my $source = shift;
    
    my $hq_output = $source . '.bed';
    
    my $command1 = $converter->create(
        source => $source,
        output => $hq_output, 
        reference_build_id => $self->reference_build_id,
        limit_variants_to => 'hq', 
    );
    
    unless($command1->execute) {
        $self->error_message('Failed to convert ' . $source . ' to the standard format.');
        return;
    }

    my $lq_output = $hq_output;
    $lq_output =~ s/\.hq\./\.lq\./ or die "file did not match expected pattern *.hq.*!: $hq_output";
    
    my $command2 = $converter->create(
        source => $source,
        output => $lq_output, 
        reference_build_id => $self->reference_build_id,
        limit_variants_to => 'lq', 
    );
    
    unless($command2->execute) {
        $self->error_message('Failed to convert ' . $source . ' to the standard format.');
        return;
    }
    return 1;
}


sub strelka_path {
    my $self = $_[0];
    return $self->path_for_strelka_version($self->version);
}

sub available_strelka_versions {
    my $self = shift;
    return keys %STRELKA_VERSIONS;
}

sub path_for_strelka_version {
    my $class = shift;
    my $version = shift;

    if (defined $STRELKA_VERSIONS{$version}) {
        return $STRELKA_VERSIONS{$version};
    }
    die('No path for strelka version '. $version);
}

sub default_strelka_version {
    die "default strelka version: $DEFAULT_VERSION is not valid" unless $STRELKA_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    if(exists($STRELKA_VERSIONS{$version})){
        return 1;
    }
    return 0;
}

1;
