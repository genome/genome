package Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam;

use strict;
use warnings;

use Genome;
use File::Compare;
use File::Basename;

class Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam {
    is  => 'Command::V2',
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
        samtools_version => {
            is  => 'Version',
            doc => 'samtools version to use',
            is_optional => 1,
            default_value => Genome::Model::Tools::Sam->default_samtools_version,
        },
        picard_version => {
            is  => 'Version',
            doc => 'Picard version to replace bam header',
            is_optional => 1,
            default_value => Genome::Model::Tools::Picard->default_picard_version,
        },
    ],
    has_transient_optional => [
        _temp_out_bam => {
            is  => 'FilePath',
            doc => 'temp output bam path',
        },
    ],
};


sub help_brief {
    'Tool to create a per lane bam from a source bam.';
}

sub help_detail {
    return <<EOS
    Tool to create a per lane bam from a source bam. The new bam need be reheadered and compare flagstat with original per lane bam. The new bam also need to revert markdup flag.
EOS
}

sub execute {
    my $self = shift;

    my $source_bam = $self->source_bam;
    my $out_bam    = $self->output_bam;
    
    unless (-s $source_bam) {
        die $self->error_message("source_bam $source_bam is not valid");
    }

    unless (Genome::Sys->validate_file_for_writing($out_bam)) {
        die $self->error_message("output_bam $out_bam is not writable");
    }
    
    $self->_temp_out_bam(Genome::Sys->create_temp_file_path(basename $out_bam));
    my ($temp_bam, $no_markdup_bam) = map{Genome::Sys->create_temp_file_path .".$_"}qw(bam no_markdup.bam);

    $self->debug_message('Extract read group bam');
    $self->_extract_readgroup_bam($temp_bam);

    $self->debug_message('Compare flagstat');
    $self->_compare_flagstat($temp_bam);

    $self->debug_message('Revert bam Markdup tag');
    $self->_revert_markdup($temp_bam, $no_markdup_bam);
    
    $self->debug_message('Reheader bam');
    $self->_reheader_bam($no_markdup_bam);
    
    $self->debug_message('Create Bam index');
    $self->_create_bam_index;

    $self->debug_message('Create Bam md5');
    $self->_create_bam_md5;

    $self->debug_message('Move outputs over');
    $self->_move_outputs;
    
    return 1;
}


sub _extract_readgroup_bam {
    my ($self, $temp_bam) = @_;

    my $extract = Genome::Model::Tools::Sam::ExtractReadGroup->create(
        input         => $self->source_bam,
        output        => $temp_bam,
        read_group_id => $self->read_group_id,
        use_version   => $self->samtools_version,
    );
    
    unless ($extract->execute) {
        die $self->error_message('Failed to run sam ExtractReadGroup');
    }
}


sub _compare_flagstat {
    my ($self, $temp_bam) = @_;

    my $out_dir = dirname $self->output_bam;

    my @flagstats = glob("$out_dir/*.flagstat");
    unless (@flagstats and @flagstats == 1) {
        die $self->error_message("Failed to get 1 flagstat file from original per lane bam location");
    }

    my $temp_flagstat = Genome::Sys->create_temp_file_path;

    my $flagstat = Genome::Model::Tools::Sam::Flagstat->create(
        bam_file    => $temp_bam,
        output_file => $temp_flagstat, 
        use_version => $self->samtools_version,  
    );

    unless ($flagstat->execute) {
        die $self->error_message("Fail to run flagstat on $temp_bam");
    }
    
    unless (compare($temp_flagstat, $flagstats[0]) == 0) {
        die $self->error_message('The bam flagstat from the extracting is different from the original per lane bam flagstat');
    }
}


sub _revert_markdup {
    my ($self, $temp_bam, $no_markdup_bam) = @_;

    Genome::Sys->shellcmd(
        cmd => "/usr/bin/seq-grind0.1 revert picard-mark-dup -i $temp_bam -o $no_markdup_bam",
        input_files  => [$temp_bam],
        output_files => [$no_markdup_bam],
        skip_if_output_is_present => 0,
    );
}


sub _reheader_bam {
    my ($self, $temp_bam) = @_;

    my $out_dir = dirname $self->output_bam;
    my @header_files = glob("$out_dir/*.header");

    unless (@header_files and @header_files == 1) {
        die $self->error_message("Failed to get 1 header file from original per lane bam location");
    }

    my $header_content = Genome::Sys->read_file($header_files[0]);
    if ($header_content =~ /(SO:unsorted)/) {
        $header_content =~ s/$1/SO:coordinate/;
    }
    
    my $header_file = Genome::Sys->create_temp_file_path; 
    Genome::Sys->write_file($header_file, $header_content);

    my $reheader = Genome::Model::Tools::Picard::ReplaceSamHeader->create(
        input_file  => $temp_bam,
        output_file => $self->_temp_out_bam, 
        header_file => $header_file,
        use_version => $self->picard_version,  
    );

    unless ($reheader->execute) {
        die $self->error_message("Fail to run picard bam reheader on $temp_bam");
    }
}


sub _create_bam_index {
    my $self = shift;
    my $out_bam = $self->_temp_out_bam;

    my $index = Genome::Model::Tools::Sam::IndexBam->create(
        bam_file       => $out_bam,
        bam_index_file => $out_bam.'.bai', 
        use_version    => $self->samtools_version,  
    );

    unless ($index->execute) {
        $self->error_message("Fail to run bam index on $out_bam");
    }
}


sub _create_bam_md5 {
    my $self = shift;
    my $out_bam = $self->_temp_out_bam;

    my $md5_content = Genome::Sys->md5sum($out_bam);
    my $basename    = basename $out_bam;
    $md5_content   .= "\t$basename\n";

    my $fh = Genome::Sys->open_file_for_writing($out_bam.'.md5');
    $fh->print($md5_content);
    $fh->close;
}


sub _move_outputs {
    my $self = shift;
    my $tmp_out = $self->_temp_out_bam;
    my $out_dir = dirname $self->output_bam;

    for my $type ('', '.bai', '.md5') {
        my $tmp_file = $tmp_out . $type;
        my $out_file = $out_dir . '/' . basename($tmp_file);
        Genome::Sys->move_file($tmp_file, $out_file);
    }
}


1;
