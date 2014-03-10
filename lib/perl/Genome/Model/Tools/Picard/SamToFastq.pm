
package Genome::Model::Tools::Picard::SamToFastq;

use strict;
use warnings FATAL => 'all';

use Genome;
use version;

class Genome::Model::Tools::Picard::SamToFastq {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input => {
            is  => 'String',
            doc => 'Input SAM/BAM file to extract reads from. Required.',
        },
        fastq => {
            is          => 'String',
            doc         => 'Output fastq file (single-end fastq or, if paired, first end of the pair fastq).',
        },
        fastq2 => {
            is          => 'String',
            doc         => 'Output fastq file (if paired, second end of the pair fastq).',
            is_optional => 1,
        },
        fragment_fastq => {
            is          => 'String',
            doc         => 'Output fastq file for bams which contain a mix of fragments & pairs -- required if paired',
            is_optional => 1,
        },
        no_orphans      => {
            is => 'Boolean',
            doc => 'Do not warn on orphaned reads (good reads, but whose mates were marked as failing quality filtering)',
            default_value => 0,
        },
        read_group_id => {
            is          => 'String',
            doc         => 'Limit to a single read group',
            is_optional => 1,
        }
    ],
};

sub help_brief {
    'Tool to create FASTQ file from SAM/BAM using Picard with added support for mixed paired/fragment BAMs';
}

sub help_detail {
    return <<EOS
    Tool to create FASTQ file from SAM/BAM using Picard.  Based on Picard's "SamToFastq" (see `gmt picard standard-sam-to-fastq`)
EOS
}

sub samtools_version { return 'r982'; }

sub execute {
    my $self = shift;

    my $picard_version = $self->use_version;

    if ($self->fastq2 && !$self->fragment_fastq) {
        $self->error_message("you must specify a fragment fastq file output if you are using pairs!");
        return;
    }

    my $input_file = $self->input;
    my $unlink_input_bam_on_end = 0;
    my $bam_read_count;

    if (defined $self->read_group_id) {
        $unlink_input_bam_on_end = 1;

        my $temp = Genome::Sys->base_temp_directory;
        my $sorted_temp_bam_file = $temp . "/temp_rg.sorted." . $$ . ".bam";
        my $extract_rg_cmd = Genome::Model::Tools::Sam::ExtractReadGroup->create(input=>$input_file, output=>$sorted_temp_bam_file, read_group_id=>$self->read_group_id);
        unless ($extract_rg_cmd->execute()) {
            $self->error_message("Failed to extract read group.");
            return;
        }
        $input_file = $sorted_temp_bam_file;
        $bam_read_count = $extract_rg_cmd->read_count;
    }

    my $picard_dir = $self->picard_path;
    my $picard_jar_path = $picard_dir . "/sam-".$picard_version.".jar";
    my $sam_jar_path = $picard_dir . "/picard-".$picard_version.".jar";

    my $picard_api_change_version = version->parse('1.77');
    my $picard_use_version = version->parse($picard_version);
    my $tool_jar_path = '';

    if($picard_use_version >= $picard_api_change_version){
      $tool_jar_path = $self->class->base_dir . "/GCSamToFastqPicard177.jar";
    }else{
      $tool_jar_path = $self->class->base_dir . "/GCSamToFastq.jar";
    }

    my $cp = join ":", ($picard_jar_path, $sam_jar_path, $tool_jar_path);

    my $jvm_options = $self->additional_jvm_options || '';
    my $java_vm_cmd = $cp . ' edu.wustl.genome.samtools.GCSamToFastq ';


    my $args = '';

    $args .= ' INPUT=' . "'" . $input_file . "'";
    $args .= ' FASTQ=' . "'" . $self->fastq . "'";
    $args .= ' SECOND_END_FASTQ=' . "'" . $self->fastq2 . "'" if ($self->fastq2);
    $args .= ' FRAGMENT_FASTQ=' . "'" . $self->fragment_fastq. "'" if ($self->fragment_fastq);
    $args .= ' NO_ORPHAN=true' if ($self->no_orphans);

    $java_vm_cmd .= $args;

    print $java_vm_cmd . "\n";

    $self->run_java_vm(
        cmd          => $java_vm_cmd,
        input_files  => [ $input_file ],
        skip_if_output_is_present => 0,
    );

    # VERIFY READ COUNTS: INPUT BAM v. FASTQS
    if ( not defined $bam_read_count ) {
        $bam_read_count = $self->_read_count_for_bam($input_file);
        return if not $bam_read_count;
    }
    my @output_files = ($self->fastq);
    push @output_files, $self->fastq2 if $self->fastq2;
    push @output_files, $self->fragment_fastq if $self->fragment_fastq;

    my $fastq_read_count = $self->_read_count_for_fastq(@output_files);
    return if not $fastq_read_count;
    $self->debug_message("VERIFY READ COUNTS: INPUT BAM v. OUTPUT FASTQ(s)");
    $self->debug_message("$bam_read_count reads in INPUT BAM: $input_file");
    $self->debug_message("$fastq_read_count reads in OUTPUT FASTQ(s): ".join(' ', @output_files));

    # RM INPUT BAM
    unlink $input_file if ($unlink_input_bam_on_end && $self->input ne $input_file);

    # COMPARE READ COUNTS
    if ( $bam_read_count ne $fastq_read_count ) {
        $self->error_message("Different number of reads in BAM ($bam_read_count) and FASTQ ($fastq_read_count)");
        return;
    }

    return 1;
}

sub _read_count_for_bam {
    my ($self, $bam) = @_;

    Carp::confess('No bam to get read count!') if not $bam;

    my $tmpdir = Genome::Sys->base_temp_directory;
    my $flagstat_file = $tmpdir.'/flagstat';
    unlink $flagstat_file;
    my $gmt = Genome::Model::Tools::Sam::Flagstat->create(
        bam_file => $bam,
        output_file => $flagstat_file,
        use_version => $self->samtools_version,
    );
    if ( not $gmt ) {
        $self->error_message('Failed to create gmt sam flagstat!');
        return;
    }
    $gmt->dump_status_messages(1);
    my $ok = $gmt->execute;
    if ( not $ok ) {
        $self->error_message('Failed to execute gmt sam flagstat!');
        return;
    }

    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
    if ( not $flagstat ) {
        $self->error_message('Failed to get metrics from flagstat file: '.$flagstat_file);
        return;
    }

    #It seems this picard tool will only return reads passing QC. No QC
    #failed reads will be put in fastq files.
    if ( not defined $flagstat->{reads_marked_passing_qc} ) {
        $self->error_message('No reads_marked_passing_qc from flagstat file!');
        return;
    }

    return $flagstat->{reads_marked_passing_qc};
}

sub _read_count_for_fastq {
    my ($self, @fastqs) = @_;

    Carp::confess('No fastq to get read count!') if not @fastqs;

    my $read_count;
    for my $fastq ( @fastqs ) {
        next if not -s $fastq;
        my $line_count = `wc -l < $fastq`;
        if ( $? or not $line_count ) {
            $self->error_message("Line count on fastq ($fastq) failed : $?");
            return;
        }

        chomp $line_count;
        if ( ($line_count % 4) != 0 ) {
            $self->error_message("Line count ($line_count) on fastq ($fastq) not divisble by 4.");
            return;
        }
        $read_count += $line_count / 4;
    }

    return $read_count;
}

1;
__END__

