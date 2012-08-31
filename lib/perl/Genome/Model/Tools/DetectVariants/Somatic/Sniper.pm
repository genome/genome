package Genome::Model::Tools::DetectVariants::Somatic::Sniper;

use warnings;
use strict;

use Genome;
use Workflow;

my $DEFAULT_VERSION = '0.7.3';
my $SNIPER_COMMAND = 'bam-somaticsniper';

class Genome::Model::Tools::DetectVariants::Somatic::Sniper {
    is => ['Genome::Model::Tools::DetectVariants::Somatic'],
    has => [
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of sniper to use",
        },
        detect_snvs => { 
            default_value => 1,
            doc => "Whether or not the tool should detect snps.  If set to false, the tool will still discover snps and indels at the same time, but will throw away any snps it detects",
            is_optional => 1,
        },
        detect_indels => { 
            default_value => 1,
            doc => "Whether or not the tool should detect indels.  If set to false, the tool will still discover snps and indels at the same time, but will throw away any indels it detects",
            is_optional => 1,
        },
        snv_params => {
            is => 'Text',
            default => '-q 1 -Q 15',
            is_input=>1, 
            doc => "Parameters for running bam-somaticsniper for snps. Since it discovers both snps and indels in one run, providing different parameters for snps and indels causes bam-somatisniper to run twice.",
        },
        indel_params => {
            is => 'Text',
            default => '-q 1 -Q 15',
            is_input=>1, 
            doc => "Parameters for running bam-somaticsniper for indels. Since it discovers both snps and indels in one run, providing different parameters for snps and indels causes bam-somatisniper to run twice.",
        },
        reference_sequence_input => {
            is  => 'String',
            is_input=>1,
            is_optional=>1, 
            default => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa', 
            doc => 'The somatic sniper reference file',
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 1,
            doc => 'enable this flag to skip this step if the output_file is already present. Useful for pipelines.',
        },
    ],
    # Make workflow choose 64 bit blades
    has_param => [
        lsf_queue => {
            default_value => 'long'
        }, 
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>100000] span[hosts=1]',
        },
    ],
    # These are params from the superclass' standard API that we do not require for this class (dont show in the help)
    has_constant_optional => [
        sv_params=>{},
        capture_set_input=>{},
    ],
};

my %SNIPER_VERSIONS = (
    '0.7' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7/' . $SNIPER_COMMAND,
    '0.7.1' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.1/' . $SNIPER_COMMAND,
    '0.7.2' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.2/' . $SNIPER_COMMAND,
    '0.7.3' => $ENV{GENOME_SW} . '/samtools/sniper/somatic_sniper-v0.7.3/' . $SNIPER_COMMAND,
    '0.7.4' => "/usr/bin/${SNIPER_COMMAND}0.7.4",
);

sub help_brief {
    "Produces a list of high confidence somatic snps and indels.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic sniper --aligned-reads-input tumor.bam --control-aligned-reads-input normal.bam --output-directory sniper
gmt somatic sniper --aligned-reads tumor.bam --control normal.bam --out sniper --quality 25
EOS
}

sub help_detail {                           
    return <<EOS 
    Provide a tumor and normal BAM file and get a list of somatic snps.  
EOS
}

sub _should_skip_execution {
    my $self = shift;
    
    if (($self->skip_if_output_present)&&(-s $self->snv_output)&&(-s $self->indel_output)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }
    
    return $self->SUPER::_should_skip_execution;
}

sub _detect_variants {
    my $self = shift;

    $self->status_message("beginning execute");

    # Validate files
    unless ( Genome::Sys->validate_file_for_reading($self->aligned_reads_input) ) {
        $self->error_message("Could not validate tumor file:  ".$self->aligned_reads_input );
        die;
    } 

    unless ( Genome::Sys->validate_file_for_reading($self->control_aligned_reads_input) ) {
        $self->error_message("Could not validate normal file:  ".$self->control_aligned_reads_input );
        die;
    } 

    # Run sniper C program... run twice if we get different sets of params for snps and indels
    my $snv_params = $self->snv_params || "";
    my $indel_params = $self->indel_params || "";
    my $result;
    if ( ($self->detect_snvs && $self->detect_indels) && ($snv_params eq $indel_params) ) {
        $result = $self->_run_sniper($snv_params, $self->_snv_staging_output, $self->_indel_staging_output);
    } else {
        # Run twice, since we have different parameters. Detect snps and throw away indels, then detect indels and throw away snps
        if ($self->detect_snvs && $self->detect_indels) {
            $self->status_message("Snp and indel params are different. Executing sniper twice: once each for snps and indels with their respective parameters");
        }
        my ($temp_fh, $temp_name) = Genome::Sys->create_temp_file();

        if ($self->detect_snvs) {
            $result = $self->_run_sniper($snv_params, $self->_snv_staging_output, $temp_name);
        }
        if ($self->detect_indels) {
            if($self->detect_snps and not $result) {
                $self->status_message('Sniper did not report success for snp detection. Skipping indel detection.')
            } else {
                $result = $self->_run_sniper($indel_params, $temp_name, $self->_indel_staging_output);
            }
        }
    }
    
    #Manually check for $self->_indel_staging_output as there might not be any indels and shellcmd()
    # chokes unless either all are present or all are empty.
    #(This means shellcmd() can check for the SNPs file on its own and still work given an empty result.)
    #Varied the warning text slightly so this message can be disambiguated from shellcmd() output in future debugging
    unless(-s $self->_indel_staging_output) {
        #Touch the file to make sure it exists
        my $fh = Genome::Sys->open_file_for_writing($self->_indel_staging_output);
        unless ($fh) {
            $self->error_message("failed to touch " . $self->_indel_staging_output . "!: " . Genome::Sys->error_message);
            die;
        }
        $fh->close;
        
        $self->warning_message("ALLOWING zero size output file " . $self->_indel_staging_output);
    }

    $self->status_message("ending execute");
    return $result; 
}

sub _run_sniper {
    my ($self, $params, $snp_output, $indel_output) = @_;
    
    my $cmd = $self->sniper_path . " " . $params . " -f ".$self->reference_sequence_input." ".$self->aligned_reads_input." ".$self->control_aligned_reads_input ." " . $snp_output . " " . $indel_output; 
    my $result = Genome::Sys->shellcmd( cmd=>$cmd, input_files=>[$self->aligned_reads_input,$self->control_aligned_reads_input], output_files=>[$snp_output], skip_if_output_is_present=>0, allow_zero_size_output_files => 1, );

    return $result;
}

sub sniper_path {
    my $self = $_[0];
    return $self->path_for_sniper_version($self->version);
}

sub available_sniper_versions {
    my $self = shift;
    return keys %SNIPER_VERSIONS;
}

sub path_for_sniper_version {
    my $class = shift;
    my $version = shift;

    if (defined $SNIPER_VERSIONS{$version}) {
        return $SNIPER_VERSIONS{$version};
    }
    die('No path for bam-somaticsniper version '. $version);
}

sub default_sniper_version {
    die "default bam-somaticsniper version: $DEFAULT_VERSION is not valid" unless $SNIPER_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}


1;
