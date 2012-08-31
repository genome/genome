package Genome::Model::Tools::Sam::StreamToBam;

use strict;
use warnings;

use Genome;
use File::Copy;
use File::Basename;


class Genome::Model::Tools::Sam::StreamToBam {
    is  => 'Genome::Model::Tools::Sam',
    has => [ 
        bam_file    => {
            is  => 'String',
            doc => 'Name of output bam file (default: use base name of input sam file -- e.g. foo.sam -> foo.bam)'
        },
        index_bam   => {
            is  => 'Boolean',
            doc => 'flag to index bam file, default yes',
            default => 1,
        },
        fix_mate    => {
            is  => 'Boolean',
            doc => 'fix mate info problem in sam/bam, default no',
            default => 0,
        },
        sam_input_stream => {
            is => 'IO::File',
            doc => 'bam stream to read from'
        },
        ref_list    => {
            is  => 'String',
            doc => 'path to a tab delimited file containing each contig name in the reference and its length',
        },
    ]
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

    return $self;
}


sub execute {
    my $self = shift;

    unless (defined $self->ref_list and -s $self->ref_list) {
        $self->error_message('Ref list does not exist');
        return;
    }

    my $samtools = $self->samtools_path;
    
    my $bam_file = $self->bam_file;
    
    my $cmd = sprintf('%s view -bt %s -o %s -', $samtools, $self->ref_list, $bam_file);
    $self->status_message("SamToBam conversion command: $cmd");
    
    my $bam_fh = IO::File->new("|$cmd");
    unless ($bam_fh) {
        $self->error_message("Failed to open pipe to sam writer!");
        return;
    } 

    my $sam_fh = $self->sam_input_stream;

    my $lines_read = 0;

    while (my $line = <$sam_fh>) {
        print $line;
        print $bam_fh $line;
    }

    $bam_fh->close;

    my $rv;
     
    #watch out disk space, for now hard code maxMemory 2000000000
    if ($self->fix_mate) {
    
        $self->status_message("Fix-mate running.  First sorting by name to gather up mates");
        my $tmp_file = $bam_file.'.sort';
        #402653184 bytes = 3 Gb 
        Genome::Sys->shellcmd(cmd=>"samtools sort -n -m 402653184 $bam_file $tmp_file",
                                              skip_if_output_is_present=>0,
                                              output_files=>["$tmp_file.bam"]);
    
        $self->status_message("Removing the unsorted orignal bam file $bam_file");
        # remove the unsorted bam file
        unlink($bam_file);

        
        $self->status_message("Running fixmate on the sorted file");
        Genome::Sys->shellcmd(cmd=>"$samtools fixmate $tmp_file.bam $tmp_file.fixmate",
                                              skip_if_output_is_present=>0,
                                              output_files=>["$tmp_file.fixmate"]);
        $self->status_message("Removing the name-sorted, but not fixmated bam file $tmp_file.bam");
        unlink "$tmp_file.bam";

        $self->status_message("Now restoring the original sort order");
        Genome::Sys->shellcmd(cmd=>"$samtools sort -m 402653184 $tmp_file.fixmate $tmp_file.fix",
                                              skip_if_output_is_present=>0,
                                              output_files=>["$tmp_file.fix.bam"]);
        $self->error_message("Sort by position failed") and return if $rv or !-s $tmp_file.'.fix.bam';
        
        $self->status_message("Now removing the fixmated, namesorted bam");
        unlink "$tmp_file.fixmate";

        move "$tmp_file.fix.bam", $bam_file;
    }

    if ($self->index_bam) {
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


        unless (Genome::Model::Tools::Sam::IndexBam->execute(
            bam_file => $bam_file,
            use_version => $self->use_version,
        )) {
            $self->error_message('Failed to index BAM file '. $bam_file);
            die($self->error_message);
        }
    }

    return 1;
}

1;
