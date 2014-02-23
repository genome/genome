package Genome::Model::Tools::Sam::SamToBam;

use strict;
use warnings;

use Genome;
use File::Copy;
use File::Basename;


class Genome::Model::Tools::Sam::SamToBam {
    is  => 'Genome::Model::Tools::Sam',
    has => [ 
        sam_file    => { 
            is  => 'String',      
            doc => 'name of sam file',
        },
    ],
    has_optional => [
        bam_file    => {
            is  => 'String',
            doc => 'Name of output bam file (default: use base name of input sam file -- e.g. foo.sam -> foo.bam)'
        },
        ref_list    => {
            is  => 'String',
            doc => 'path to a tab delimited file containing each contig name in the reference and its length',
        },
        index_bam   => {
            is  => 'Boolean',
            doc => 'flag to index bam file, default yes',
            default => 1,
        },
        is_sorted => {
            is  => 'Boolean',
            doc => 'flag bam file is sorted, default yes',
            default => 1,
        },
        fix_mate    => {
            is  => 'Boolean',
            doc => 'fix mate info problem in sam/bam, default no',
            default => 0,
        },
        keep_sam    => {
            is  => 'Boolean',
            doc => 'flag to keep sam file, default no',
            default => 0,
        },
    ],
};


sub help_brief {
    'create bam file from sam file';
}


sub help_detail {
    return <<EOS 
This tool makes bam file from samfile with options to index bam file, fix mate pair info.
EOS
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->error_message('Sam file('.$self->sam_file.') does not exist') and return unless -e $self->sam_file;
	if ($self->ref_list) {
	    $self->error_message('Ref list('.$self->ref_list.') does not have size') and return unless -s $self->ref_list;
	}

    return $self;
}


sub execute {
    my $self = shift;

    my $samtools = $self->samtools_path;
    my $sam_file = $self->sam_file;
	my $sam_read_count = $self->read_count($sam_file);
    
    my ($root_name) = basename $sam_file =~ /^(\S+)\.sam/;
    
    my $sam_dir  = dirname $sam_file;
    my $bam_file = $self->bam_file || $sam_dir . "/$root_name.bam";
    
    my $cmd;
	$cmd .= sprintf('%s view -b -o %s', $samtools, $bam_file);
    if ($self->ref_list) {
        $cmd .= sprintf(' -t %s %s', $self->ref_list, $sam_file);
    } else {
        $cmd .= sprintf(' -S %s', $sam_file);
    }
    $self->debug_message("SamToBam conversion command: $cmd");
    
    my $rv  = Genome::Sys->shellcmd(
        cmd => $cmd, 
        output_files => [$bam_file],
        skip_if_output_is_present => 0,
    );
        
    $self->error_message("Converting to Bam command: $cmd failed") and return unless $rv == 1;

	my $bam_read_count = $self->read_count($bam_file);
	unless ( $bam_read_count == $sam_read_count ) {
		$self->error_message("Read counts differ after SAM to BAM conversion! (BAM: $bam_read_count vs. SAM: $sam_read_count)");
		return;
	}
     
    #watch out disk space, for now hard code maxMemory 2000000000
    if ($self->fix_mate) {
        my $tmp_file = $bam_file.'.sort';
        #402653184 bytes = 3 Gb 
        $rv = system "$samtools sort -n -m 402653184 $bam_file $tmp_file";
        $self->error_message("Sort by name failed") and return if $rv or !-s $tmp_file.'.bam';

        $rv = system "$samtools fixmate $tmp_file.bam $tmp_file.fixmate";
        $self->error_message("fixmate failed") and return if $rv or !-s $tmp_file.'.fixmate';
        unlink "$tmp_file.bam";

        $rv = system "$samtools sort -m 402653184 $tmp_file.fixmate $tmp_file.fix";
        $self->error_message("Sort by position failed") and return if $rv or !-s $tmp_file.'.fix.bam';
        
        unlink "$tmp_file.fixmate";
        unlink $bam_file;

        move "$tmp_file.fix.bam", $bam_file;
    }else{
        $self->is_sorted(0);
    }

    if ($self->index_bam) {
        unless ($self->is_sorted) {
            #This may not work for large genome BAM files
            my $tmp_bam = Genome::Sys->create_temp_file_path();
            unless (Genome::Model::Tools::Sam::SortBam->execute(
                file_name => $bam_file,
                output_file => $tmp_bam,
                use_version => $self->use_version,
            )) {
                $self->error_message('Failed to sort bam file '. $bam_file);
                die($self->error_message);
            }
            $tmp_bam .= '.bam';
            unless (unlink($bam_file)) {
                $self->error_message('Failed to remove unsorted bam file '. $bam_file .":  $!");
                die($self->error_message);
            }
            unless (move($tmp_bam,$bam_file)) {
                $self->error_message('Failed to move the sorted bam file from '. $tmp_bam .' to '. $bam_file .":  $!");
                die($self->error_message);
            }
        }
        unless (Genome::Model::Tools::Sam::IndexBam->execute(
            bam_file => $bam_file,
            use_version => $self->use_version,
        )) {
            $self->error_message('Failed to index BAM file '. $bam_file);
            die($self->error_message);
        }
    }

    unlink $sam_file unless $self->keep_sam;
    return 1;
}

1;
