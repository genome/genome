
package Genome::Model::Tools::Sam::InferPairedStatus;

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Model::Tools::Sam::InferPairedStatus {
    is  => 'Genome::Model::Tools::Sam',
    has_input => [
        input => {
            is  => 'String',
            doc => 'Input SAM/BAM file to extract reads from. Required.',
        },
    ],
    has_output => [
        is_paired_end   => { is => 'Boolean',
                        is_optional => 1 
        },
    ],
};

sub help_brief {
    'Tool to infer if a BAM is paired end or not'
}

sub help_detail {
    return <<EOS
    Tool to infer if a BAM is paired end or not by taking a quick sample of the reads, and using flagstat.  NOT VALID on merged BAMs which may contain read groups from different sources.
EOS
}

sub samtools_version { return 'r982'; }

sub execute {
    my $self = shift;
    my $bam = $self->input;
    Carp::confess('No bam to get read count!') if not $bam;

    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $temp = Genome::Sys->base_temp_directory;
    my $temp_bam_file = $temp . "/temp_smp." . $$ . ".bam";

    my $header_line_count = `samtools view -H $bam | wc -l`;
    if (!$header_line_count || ! $header_line_count > 0) {
        $self->error_message("Could not determine BAM header size");
        return;
    }

    # grab a chunk of reads from the BAM to sample
    my $sampler_length = $header_line_count + 10_000;

    my $samtools_sampler_cmd = sprintf(
        "%s view -h %s | head -n %s | %s view -S -b -o %s -",
        $samtools_path,
        $bam, 
        $sampler_length,
        $samtools_path,
        $temp_bam_file,
    );

    Genome::Sys->shellcmd(
        cmd=>$samtools_sampler_cmd, 
        output_files=>[$temp_bam_file],
        skip_if_output_is_present=>0,
    );

    my $tmpdir = Genome::Sys->base_temp_directory;
    my $flagstat_file = $tmpdir.'/flagstat';
    unlink $flagstat_file;
    my $gmt = Genome::Model::Tools::Sam::Flagstat->create(
        bam_file => $temp_bam_file,
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

    $self->is_paired_end($flagstat->{'reads_paired_in_sequencing'} > 0);
    my $friendly_paired_status = $self->is_paired_end? 'PAIRED' : 'FRAGMENT';
    $self->debug_message("Inferred paired status: $friendly_paired_status");

    return 1;
}

1;
__END__

