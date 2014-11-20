package Genome::Model::Tools::Sam::ExtractReadGroupBam;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Sam::ExtractReadGroupBam {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        source_bam => {
            is  => 'FilePath',
            doc => 'The bam file to create per read group bam, like a merged bam',
        },
        output_bam => {
            is  => 'FilePath',
            doc => 'The output bam file created from source bam, like recreating a per lane bam',
        },
        read_group_id => {
            is  => 'String',
            doc => 'The read group id to create per read group bam, like a instrument data id',
        },
        flagstat => {
            is  => 'Boolean',
            doc => 'Option to run flagstat after',
            is_optional   => 1,
            default_value => 0,
        },
        index_bam => {
            is  => 'Boolean',
            doc => 'Option to run bam index after',
            is_optional   => 1,
            default_value => 0,
        },
        md5 => {
            is  => 'Boolean',
            doc => 'Option to run md5sum on bam file',
            is_optional   => 1,
            default_value => 0,
        },
    ],
};


sub help_brief {
    'Tool to create a per read group bam from a source bam.';
}

sub help_detail {
    return <<EOS
    Tool to create a per read group bam from a source bam.
EOS
}

sub execute {
    my $self = shift;

    my $samtools   = $self->samtools_path;
    my $source_bam = $self->source_bam;
    my $out_bam    = $self->output_bam;
    my $rg_info    = "<(echo '" .$self->read_group_id ."')";

    Genome::Sys->validate_file_for_reading($source_bam);
    unless (-s $source_bam) {
        die $self->error_message("source_bam $source_bam is not valid");
    }

    Genome::Sys->validate_file_for_writing($out_bam);
    
    my $cmd = $samtools .' view -h -R ' .$rg_info .' -b -o ' .$out_bam .' '. $source_bam;

    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$out_bam],
        skip_if_output_is_present => 0,
    );

    if ($self->md5) {
        my $md5_content = Genome::Sys->md5sum($out_bam);
        my $basename    = basename($out_bam);
        $md5_content   .= "\t$basename\n";
        my $fh = Genome::Sys->open_file_for_writing($out_bam.'.md5');
        $fh->print($md5_content);
        $fh->close;
    }

    if ($self->flagstat) {
        my $flagstat = Genome::Model::Tools::Sam::Flagstat->create(
            bam_file    => $out_bam,
            output_file => $out_bam.'.flagstat', 
            use_version => $self->use_version,  
        );
        
        unless ($flagstat->execute) {
            die $self->error_message("Fail to run flagstat on $out_bam");
        }
    }

    if ($self->index_bam) {
        my $index = Genome::Model::Tools::Sam::IndexBam->create(
            bam_file       => $out_bam,
            bam_index_file => $out_bam.'.bai', 
            use_version    => $self->use_version,  
        );

        unless ($index->execute) {
            die $self->error_message("Fail to run bam index on $out_bam");
        }
    }
    
    return 1;
}

1;
