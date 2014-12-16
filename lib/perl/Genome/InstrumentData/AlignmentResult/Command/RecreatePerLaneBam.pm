package Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam;

use strict;
use warnings;

use Genome;
use File::Compare;
use File::Basename;

class Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam {
    is  => 'Command::V2',
    has => [
        merged_bam => {
            is  => 'FilePath',
            doc => 'The merged bam to be used to create per lane bam',
        },
        per_lane_bam => {
            is  => 'FilePath',
            doc => 'The output per lane bam',
        },
        instrument_data_id => {
            is  => 'String',
            doc => 'The per lane instrument data id used as read group id in merged bam',
        },
        samtools_version => {
            is  => 'Version',
            doc => 'samtools version to use',
            valid_values => [Genome::Model::Tools::Sam->available_samtools_versions],
        },
        picard_version => {
            is  => 'Version',
            doc => 'Picard version to replace bam header',
            valid_values => [Genome::Model::Tools::Picard->available_picard_versions],
        },
    ],
    has_transient_optional => [
        _temp_out_bam => {
            is  => 'FilePath',
            doc => 'temp output bam path',
        },
        _output_dir => {
            is  => 'Text',
            doc => 'the directory of per lane bam',
        },
    ],
};


sub help_brief {
    'Command to recreate a per lane bam from a merged bam.';
}

sub help_detail {
    return <<EOS
    This command will recreate a per lane bam from the merged bam if this per lane bam is needed again for new bam merge. When removing per lane bam, we mean to remove all_sequences.bam along with the .bai and .md5 files, but keep .flagstat and make all_sequences.bam.header in the directory. When recreating the per lane bam, we run flagstat on the temp extract read-group bam (made from merged bam) and compare this flagstat with the original all_sequences.bam.flagstat to make sure they are exactly same. Then we revert markdup tag in this bam and reheader this bam using all_sequences.bam.header created previously in order to make this bam resembling the original per lane bam as much as possible. 
EOS
}

sub execute {
    my $self = shift;

    my $merged_bam   = $self->merged_bam;
    my $per_lane_bam = $self->per_lane_bam;
    
    unless (-s $merged_bam) {
        die $self->error_message("merged_bam $merged_bam is not valid");
    }

    unless (Genome::Sys->validate_file_for_writing($per_lane_bam)) {
        die $self->error_message("per_lane_bam $per_lane_bam is not writable");
    }
    
    $self->_temp_out_bam(Genome::Sys->create_temp_file_path(basename $per_lane_bam));
    $self->_output_dir(dirname $self->per_lane_bam);

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
        input         => $self->merged_bam,
        output        => $temp_bam,
        read_group_id => $self->instrument_data_id,
        use_version   => $self->samtools_version,
    );
    
    unless ($extract->execute) {
        die $self->error_message('Failed to run sam ExtractReadGroup');
    }
}


sub _compare_flagstat {
    my ($self, $temp_bam) = @_;
    my $out_dir = dirname $self->per_lane_bam;

    my @flagstats = glob($self->_output_dir ."/*.flagstat");
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
    
    #Temporarily place seq-grind0.1 here, it will be immediately
    #updated once it is deployed
    Genome::Sys->shellcmd(
        cmd => "/gscuser/fdu/bin/seq-grind0.1 revert picard-mark-dup -i $temp_bam -o $no_markdup_bam",
        input_files  => [$temp_bam],
        output_files => [$no_markdup_bam],
        skip_if_output_is_present => 0,
    );
}


sub _reheader_bam {
    my ($self, $temp_bam) = @_;

    my @header_files = glob($self->_output_dir ."/*.header");
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

    Genome::Sys->write_file($out_bam.'.md5', $md5_content);
}


sub _move_outputs {
    my $self = shift;
    my $tmp_out = $self->_temp_out_bam;

    for my $type ('', '.bai', '.md5') {
        my $tmp_file = $tmp_out . $type;
        my $out_file = File::Spec->join($self->_output_dir, basename $tmp_file);
        Genome::Sys->move_file($tmp_file, $out_file);
    }
}


1;
