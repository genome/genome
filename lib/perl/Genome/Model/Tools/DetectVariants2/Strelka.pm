package Genome::Model::Tools::DetectVariants2::Strelka;

use warnings;
use strict;

use Genome;
use Workflow;

use Config::IniFiles;

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

our %STRELKA_VERSIONS = (
    #'0.4.6.2' => '/usr/lib/strelka0.4.6.2/',
    #'0.4.10.1' => '/usr/lib/strelka0.4.10.1/',
    map { 
        /strelka(.*)/;
        ($1 => $_);
    } glob("/usr/lib/strelka*")
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

    $self->debug_message("beginning execute");

    my $output_dir = $self->_temp_staging_directory;

    # Note that all additional arguments to Strelka must be passed to Strelka as a single string using the 'extraStrelkaArguments =' line
    # For example, you could do:
    # extraStrelkaArguments = -used-allele-count-min-qscore 30 -min-qscore 10
    # For a full list of Strelka options, run the strelka binary to view the help docs.
    # /gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/strelka/bin/strelka
    # Note that there are many possible additional arguments!

    # Update the default parameters with those passed in.
    my $strelka_path = $self->strelka_path;
    my $default_config_filename;
    if (-d $strelka_path . '/etc') {
        $default_config_filename = $strelka_path . '/etc/strelka_config_bwa_default.ini';
    }
    else {
        $default_config_filename = join("/", $strelka_path, qw(strelka_workflow strelka etc strelka_config_bwa_default.ini));
    }

    my $working_config_filename = join("/", $output_dir, "strelka_config.ini");
    my %params = parse_params($self->params);
    my $config_file = Config::IniFiles->new(-file=>$default_config_filename);
    for my $key (keys %params) {
        unless($config_file->setval('user', $key, $params{$key})) {
            $self->error_message("$key is an invalid parameter to Strelka.");
            Carp::croak($self->error_message());
        }
    }
    $config_file->WriteConfig($working_config_filename);

    #Run the strelka configuration step that checks your input files and prepares a Makefile
    my $config_path = $self->strelka_path . "/strelka_workflow/configureStrelkaWorkflow.pl";
    unless (-e $config_path) {
        $config_path = $self->strelka_path . "/bin/configureStrelkaWorkflow.pl";
        unless (-e $config_path) {
            die "failed to find configureStrelkaWorkflow.pl under " . $self->strelka_path . " in either bin or strelka_workflow";
        }
    }

    my $cmd = $config_path  . " --tumor " . $self->aligned_reads_input
                            . " --normal " . $self->control_aligned_reads_input
                            . " --ref " . $self->reference_sequence_input
                            . " --config $working_config_filename"
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

    $self->debug_message("ending execute");
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

sub parse_params {
    my ($string) = @_;

    my @kv_pairs = split(";", $string);
    my %result;
    for my $kv_pair (@kv_pairs) {
        my ($key, $value) = split(/\s*=\s*/, $kv_pair);
        $result{$key} = $value;
        }
    return %result;
}

1;

