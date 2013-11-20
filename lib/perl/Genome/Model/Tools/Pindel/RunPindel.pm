package Genome::Model::Tools::Pindel::RunPindel;

use warnings;
use strict;

use Genome;
use Workflow;

my $DEFAULT_VERSION = '0.5';
my $PINDEL_COMMAND = 'pindel';

use Genome::Utility::Text;

class Genome::Model::Tools::Pindel::RunPindel {
    is => ['Command'],
    has => [
        aligned_reads_input => {
            is => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'This is the bam that you want pindel to call indels in',
        },
        output_directory => {
            is => 'String',
            is_optional => 0,
            is_input => 1,
            is_output => 1,
            doc => 'Where the output should go.'
        },
        version => {
            is => 'Version',
            is_optional => 0,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of pindel to use",
        },
        window_size => {
            is => 'String',
            doc => 'This is the size (in millions of bases) of the window which pindel will load in bam reads by. Default is 10',
            default => '10',
            is_input => 1,
        },
        detect_indels => { 
            default_value => 1,
            doc => "Whether or not the tool should detect indels. ",
            is_optional => 1,
        },
        reference_build_id => {
            is => 'Text',
            doc => 'The build-id of a reference sequence build',
            is_input => 1,
        },
        reference_sequence_input => {
            calculate_from => ['reference_build_id'],
            calculate => q{ Genome::Model::Build->get($reference_build_id)->cached_full_consensus_path('fa') },
            doc => 'Location of the reference sequence file',
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
        },
        _temp_long_insertion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_LI"); },
            doc => "Where the long insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_short_insertion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_SI"); },
            doc => "Where the short insertion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_deletion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_D"); },
            doc => "Where the deletion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_inversion_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_INV"); },
            doc => "Where the inversion output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_tandem_duplication_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_TD"); },
            doc => "Where the tandem duplication output should be in the temp staging directory. This is contingent upon code in pindel and will change if the pindel binary does.",
        },
        _temp_breakpoint_output => {
            calculate_from => ['output_directory'],
            calculate => q{ join("/", $output_directory, "all_sequences_BP"); },
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
            default_value => 'apipe'
        }, 
        lsf_resource => {
            default_value => "-M 12000000 -R 'select[type==LINUX64 && mem>12000] rusage[mem=12000]'",
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
    '0.2'    => $ENV{GENOME_SW} . '/pindel/0.2/pindel',
    '0.2.4d' => $ENV{GENOME_SW} .'/pindel/0.2.4d/pindel',
    '0.2.4t' => '/usr/bin/pindel0.2.4t',
    '0.2.5' => '/usr/bin/pindel0.2.5',
    '0.3'    => $ENV{GENOME_SW} .'/pindel/0.3/pindel',
    '0.4'    => $ENV{GENOME_SW} .'/pindel/0.4/pindel',
    '0.5'    => $ENV{GENOME_SW} .'/pindel/0.5/pindel',
);

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
        my $sanitized_chromosome = Genome::Utility::Text::sanitize_string_for_filesystem($self->chromosome);
        mkdir $self->output_directory . '/' . $sanitized_chromosome;
        $self->output_directory($self->output_directory . '/' . $sanitized_chromosome);
    }

    return $self;
}

sub execute {
    my $self = shift;

    unless($self->_verify_inputs) {
        die $self->error_message('Failed to verify inputs.');
    }

    unless($self->_create_directories) {
        die $self->error_message('Failed to create directories.');
    }

    unless ($self->_detect_variants) {
        die $self->error_message("Failed to get a return value from _detect_variants.");
    }

    return 1;
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
    if(defined($self->control_aligned_reads_input)){
        my $control_aligned_reads_file = $self->control_aligned_reads_input;
        unless(Genome::Sys->check_for_path_existence($control_aligned_reads_file)) {
            $self->error_message("control_aligned reads input $control_aligned_reads_file was not found.");
            return;
        }
    }

    $self->status_message("Completed verify_inputs step.");

    return 1;
}


sub _create_directories {
    my $self = shift;

    my $output_directory = $self->output_directory;
    unless (-d $output_directory) {
        eval {
            Genome::Sys->create_directory($output_directory);
        };

        if($@) {
            $self->error_message($@);
            return;
        }

        $self->status_message("Created directory: $output_directory");
        chmod 02775, $output_directory;
    }

    #$self->_temp_staging_directory(Genome::Sys->create_temp_directory);
    #$self->_temp_scratch_directory(Genome::Sys->create_temp_directory);

    return 1;
}

# The config file that is internally generated to store bams, average insert size, and tag
sub _config_file {
    my $self = shift;
    return $self->output_directory . "/pindel.config";
}

# TODO hardcoded for now, but look at bam headers soon or go back to the old method of 
# calculating via a model id and instrument data as in gmt pindel run-pindel
sub _calculate_average_insert_size { 
    return "400";
}

sub _output_basename_for_chrom {
    my $self = shift;
    my $chromosome = shift;
    return join("/", $self->output_directory, Genome::Utility::Text::sanitize_string_for_filesystem($chromosome));
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

    my $result = undef;
    if ($self->chromosome) {
        $result  = $self->_run_pindel_for_chromosome($self->chromosome);

        ## this is a hack, because the rest of the DetectVariants wants things to
        ## be named like all_sequences, there's probably a cleaner way to do this
        ## Eric
        my $chromosome = $self->chromosome;
        my $sanitized_chromosome = Genome::Utility::Text::sanitize_string_for_filesystem($chromosome);
        for my $target_file ($self->_temp_short_insertion_output, $self->_temp_long_insertion_output, $self->_temp_deletion_output, 
                             $self->_temp_tandem_duplication_output, $self->_temp_inversion_output, $self->_temp_breakpoint_output) {
            my $chunk_file = $target_file;
            $chunk_file =~ s/all_sequences/$sanitized_chromosome/;
            unless (-e $chunk_file) {
                $self->error_message("Output for chromosome $chromosome: file $chunk_file does not appear to exist"); 
                die;
            }

            File::Copy::move($chunk_file,$target_file) or die $!;
        }

        # Put the insertions and deletions where the rest of the pipe expects them 
        my $files_to_cat = join(" ", ($self->_temp_short_insertion_output, $self->_temp_deletion_output) );

        my $cmd = "cat $files_to_cat > " . $self->output_directory."/indels.hq";
        unless (system($cmd) == 0) {
            $self->error_message("Problem running $cmd");
            die;
        }
    } else {
        $result = $self->_run_pindel($self->output_directory);
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
    my $cmd = "cat $files_to_cat > " . $self->output_directory;
    unless (system($cmd) == 0) {
        $self->error_message("Problem running $cmd");
        die;
    }

    return 1;
}

sub _run_pindel_for_chromosome {
    my $self = shift;
    my $chromosome = shift;

    my $sanitized_chromosome = Genome::Utility::Text::sanitize_string_for_filesystem($chromosome);
    my $reference_sequence_for_chrom = $self->reference_sequence_input;
    $reference_sequence_for_chrom =~ s/all_sequences/$sanitized_chromosome/;
    $reference_sequence_for_chrom =~ s/fasta$/fa/;
    $reference_sequence_for_chrom =~ s/\/opt\/fscache//;
    unless (-e $reference_sequence_for_chrom) {
        $self->status_message("No per-chromosome reference fasta found, falling back to all_sequences.fa");
        $reference_sequence_for_chrom = $self->reference_sequence_input;
        $reference_sequence_for_chrom =~ s/fasta$/fa/;
        $reference_sequence_for_chrom =~ s/\/opt\/fscache//;
        unless(-s $reference_sequence_for_chrom){
            die $self->error_message("Unable to locate either a per chromosome fasta or the all_sequences fasta, shutting down.");
        }
    }
    my $output_basename = $self->_output_basename_for_chrom($chromosome);
    my $window_size = $self->window_size;
    my $cmd = $self->pindel_path . " -f ".$reference_sequence_for_chrom. " -i " . $self->_config_file . " -o ". $output_basename . " -c '".$chromosome . "' -w ".$window_size." -b /dev/null";
    
    my $result;
    if(defined($self->control_aligned_reads_input)){
        $result = Genome::Sys->shellcmd( cmd=>$cmd, input_files=> [$self->aligned_reads_input, $self->control_aligned_reads_input]);
    } 
    else {
        $result = Genome::Sys->shellcmd( cmd=>$cmd, input_files=>[$self->aligned_reads_input]);
    }

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

1;
