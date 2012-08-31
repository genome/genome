package Genome::Model::Tools::DetectVariants::Somatic::Pindel;

use warnings;
use strict;

use Genome;
use Workflow;

my $DEFAULT_VERSION = '0.1';
my $PINDEL_COMMAND = 'pindel_64';

class Genome::Model::Tools::DetectVariants::Somatic::Pindel {
    is => ['Genome::Model::Tools::DetectVariants::Somatic'],
    has => [
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of pindel to use",
        },
        detect_indels => { 
            default_value => 1,
            doc => "Whether or not the tool should detect indels. ",
            is_optional => 1,
        },
        reference_sequence_input => {
            is  => 'String',
            is_input=>1,
            is_optional=>1, 
            default => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa', 
            doc => 'The somatic sniper reference file',
        },
        chromosome => {
            is  => 'String',
            is_input=>1,
            is_optional=>1, 
            doc => 'Run pindel with this chromosome. If not set, pindel will run in serial for all chromosomes.',
        },
        control_aligned_reads_input => {
            is => 'Text',
            doc => 'Location of the control aligned reads file to which the input aligned reads file should be compared',
            is_optional => 1,
            shell_args_position => '2',
            is_input => 1,
            is_output => 1,
        },
        # Temporary output files
        _temp_long_insertion_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_LI"); },
            doc => "Where the long insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_short_insertion_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_SI"); },
            doc => "Where the short insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_deletion_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_D"); },
            doc => "Where the deletion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_inversion_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_INV"); },
            doc => "Where the inversion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_tandem_duplication_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_TD"); },
            doc => "Where the tandem duplication output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_breakpoint_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, "all_sequences_BP"); },
            doc => "Where the breakpoing output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        # output files
        long_insertion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_LI"); },
            doc => "Where the long insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        short_insertion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_SI"); },
            doc => "Where the short insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        deletion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_D"); },
            doc => "Where the deletion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        inversion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_INV"); },
            doc => "Where the inversion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        tandem_duplication_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_TD"); },
            doc => "Where the tandem duplication output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        breakpoint_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_BP"); },
            doc => "Where the breakpoint output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
    ],
    # Make workflow choose 64 bit blades
    has_param => [
        lsf_queue => {
            default_value => 'long'
        }, 
        lsf_resource => {
            default_value => "-M 16000000 -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]'",
        },
    ],
    # These are params from the superclass' standard API that we do not require for this class (dont show in the help)
    has_constant_optional => [
        sv_params=>{},
        sv_version=>{},
        snv_params=>{},
        snv_version=>{},
        indel_params=>{}, # pindel does not currently accept any processing parameters that this tool does not provide
        capture_set_input=>{},
    ],
};

# FIXME make this the real deployed path
my %PINDEL_VERSIONS = (
    '0.1' => '/gscmnt/sata921/info/medseq/Pindel_test/' . $PINDEL_COMMAND,
    '0.2' => '/gscmnt/sata921/info/medseq/Pindel_test/maq/'.$PINDEL_COMMAND,     # This version and up works with MAQ aligned bams
    '0.3' => '/gscmnt/sata921/info/medseq/Pindel_test/merged_with_kai/'.$PINDEL_COMMAND,    #this version is merged with changes from KAI
    '0.4' => '/gscmnt/sata921/info/medseq/Pindel_test/version_4/'.$PINDEL_COMMAND,
    '0.5' => '/gscmnt/sata921/info/medseq/Pindel_test/pindel_v0.5/'.$PINDEL_COMMAND,
);

sub help_brief {
    "Discovers somatic indels when provided a control bam (usually a normal sample) and a comparison bam (usually a tumor sample).";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants somatic pindel --aligned-reads-input tumor.bam --control-aligned-reads-input normal.bam --output-directory pindel 
EOS
}

sub help_detail {                           
    return <<EOS 
    Provide a tumor and normal BAM file and get a list of somatic indels.  
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    if ($self->chromosome) {
        mkdir $self->output_directory . '/' . $self->chromosome;
        $self->output_directory($self->output_directory . '/' . $self->chromosome);
    }

    return $self;
}

# The config file that is internally generated to store bams, average insert size, and tag
sub _config_file {
    my $self = shift;
    return $self->_temp_staging_directory . "/pindel.config";
}

# TODO hardcoded for now, but look at bam headers soon or go back to the old method of 
# calculating via a model id and instrument data as in gmt pindel run-pindel
sub _calculate_average_insert_size { 
    return "400";
}

sub _output_basename_for_chrom {
    my $self = shift;
    my $chromosome = shift;
    return join("/", $self->_temp_staging_directory, $chromosome);
}

# Generate the config file that is required by pindel. It needs one line each for normal and tumor. The format is tab delimited and contains (per line): bam_name, avg_insert_size, tag (normal or tumor)
sub _generate_config_file { 
    my $self = shift;

    my $config_path = $self->_config_file;
    my $config_fh = Genome::Sys->open_file_for_writing($config_path);
    unless ($config_fh) {
        $self->error_message("Could not open $config_path for writing");
        die;
    }
    if($self->control_aligned_reads_input){
        my $normal_line = join("\t", ($self->control_aligned_reads_input, $self->_calculate_average_insert_size, "normal") );
        $config_fh->print("$normal_line\n");
    }

    my $tumor_line = join("\t", ($self->aligned_reads_input, $self->_calculate_average_insert_size, "tumor") );
    $config_fh->print("$tumor_line\n");

    $config_fh->close;

    unless (-s $self->_config_file) {
        $self->error_message("Config file: " . $self->_config_file . " does not exist or has zero size");
        die;
    }

    return 1;
}

sub _detect_variants {
    my $self = shift;
    # test architecture to make sure we can run
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }
    $self->_generate_config_file;

    my $result;
    if ($self->chromosome) {
        $result  = $self->_run_pindel_for_chromosome($self->chromosome);

        ## this is a hack, because the rest of the DetectVariants wants things to
        ## be named like all_sequences, there's probably a cleaner way to do this
        ## Eric
        my $chromosome = $self->chromosome;
        for my $target_file ($self->_temp_short_insertion_output, $self->_temp_long_insertion_output, $self->_temp_deletion_output, 
                             $self->_temp_tandem_duplication_output, $self->_temp_inversion_output, $self->_temp_breakpoint_output) {
            my $chunk_file = $target_file;
            $chunk_file =~ s/all_sequences/$chromosome/;
            unless (-e $chunk_file) {
                $self->error_message("Output for chromosome $chromosome: file $chunk_file does not appear to exist"); 
                die;
            }

            rename($chunk_file,$target_file) or die $!;
        }

        # Put the insertions and deletions where the rest of the pipe expects them 
        my $files_to_cat = join(" ", ($self->_temp_short_insertion_output, $self->_temp_deletion_output) );

        my $cmd = "cat $files_to_cat > " . $self->_indel_staging_output;
        unless (system($cmd) == 0) {
            $self->error_message("Problem running $cmd");
            die;
        }
    } else {
        $result = $self->_run_pindel($self->_indel_staging_output);
    }
    
    return $result; 
}

sub _run_pindel {
    my $self = shift;
    for my $chromosome (1..22, "X", "Y") {
        unless ( $self->_run_pindel_for_chromosome($chromosome) ) {
            $self->error_message("Failed to run pindel for chromosome $chromosome");
            die;
        }
        # Append this chromosome's output files to the rest
        for my $target_file ($self->_temp_short_insertion_output, $self->_temp_long_insertion_output, $self->_temp_deletion_output, 
                             $self->_temp_tandem_duplication_output, $self->_temp_inversion_output, $self->_temp_breakpoint_output) {
            my $chunk_file = $target_file;
            $chunk_file =~ s/all_sequences/$chromosome/;
            unless (-e $chunk_file) {
                $self->error_message("Output for chromosome $chromosome: file $chunk_file does not appear to exist"); 
                die;
            }

            my $cmd = "cat $chunk_file >> $target_file";
            unless (system($cmd) == 0) {
                $self->error_message("Problem running $cmd");
                die;
            }
        }

        # Clean up the per chromosome files
        my $output_basename = $self->_output_basename_for_chrom($chromosome);
        unlink glob $output_basename . "_*";
    }

    # Put the insertions and deletions where the rest of the pipe expects them 
    my $files_to_cat = join(" ", ($self->_temp_short_insertion_output, $self->_temp_deletion_output) );
    my $cmd = "cat $files_to_cat > " . $self->_indel_staging_output;
    unless (system($cmd) == 0) {
        $self->error_message("Problem running $cmd");
        die;
    }

    return 1;
}

sub _run_pindel_for_chromosome {
    my $self = shift;
    my $chromosome = shift;

    my $reference_sequence_for_chrom = $self->reference_sequence_input;
    $reference_sequence_for_chrom =~ s/all_sequences/$chromosome/;
    $reference_sequence_for_chrom =~ s/fasta$/fa/;
    unless (-e $reference_sequence_for_chrom) {
        $self->error_message("Reference sequence file $reference_sequence_for_chrom does not exist");
        die;
    }
    my $output_basename = $self->_output_basename_for_chrom($chromosome);
    my $cmd = $self->pindel_path . " -f $reference_sequence_for_chrom" . " -i " . $self->_config_file . " -o $output_basename" . " -c $chromosome" . " -b /dev/null";
    my $result = Genome::Sys->shellcmd( cmd=>$cmd, input_files=>[$self->aligned_reads_input]);

    unless($result) {
        $self->error_message("Running pindel failed with command $cmd");
        die;
    }

    return $result;
}

sub pindel_path {
    my $self = $_[0];
    return $self->path_for_pindel_version($self->version);
}

sub available_pindel_versions {
    my $self = shift;
    return keys %PINDEL_VERSIONS;
}

sub path_for_pindel_version {
    my $class = shift;
    my $version = shift;

    if (defined $PINDEL_VERSIONS{$version}) {
        return $PINDEL_VERSIONS{$version};
    }
    die('No path for pindel version '. $version);
}

sub default_pindel_version {
    die "default pindel version: $DEFAULT_VERSION is not valid" unless $PINDEL_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub _verify_inputs {
    my $self = shift;
    
    my $ref_seq_file = $self->reference_sequence_input;
    unless(Genome::Sys->check_for_path_existence($ref_seq_file)) {
        $self->error_message("reference sequence input $ref_seq_file does not exist");
        return;
    }

    my $aligned_reads_file = $self->aligned_reads_input;
    unless(Genome::Sys->check_for_path_existence($aligned_reads_file)) {
        $self->error_message("aligned reads input $aligned_reads_file was not found.");
        return;
    }
    
    return 1;
}

1;
