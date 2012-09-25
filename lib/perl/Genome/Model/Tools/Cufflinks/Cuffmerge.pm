package Genome::Model::Tools::Cufflinks::Cuffmerge;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cufflinks::Cuffmerge {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        input_gtf_paths => {
            doc => 'A list of GTF paths used as the input(s)',
        },
        output_directory => {
            is => 'Text',
            doc => 'Directory where merged assembly will be written',
            is_optional => 1,
        },
        reference_gtf_path => {
            is => 'Text',
            doc => 'a set of known mRNAs to use as a reference for assessing the accuracy of mRNAs or gene models given in <input.gtf>',
            is_optional => 1,
        },
        reference_fasta_path => {
            is => 'Text',
            doc => 'Can be a multi-fasta file with all the genomic sequences or a directory containing multiple single-fasta files (one file per contig); lower case bases will be used to classify input transcripts as repeats',
            is_optional => 1,
        },
        min_isoform_fraction => {
            is => 'Number',
            doc => 'Discard isoforms with abundance below this. Cuffmerge default of 0.05',
            is_optional => 1,
        },
        num_threads => {
            is => 'Integer',
            doc => 'Use this many threads to merge assemblies.',
            is_optional => 1,
        },
        keep_tmp => {
            is => 'Boolean',
            doc => 'Keep all intermediate files during merge',
            default_value => 0,
            is_optional => 1,
        },
    ],
    has_output_optional => [
    ],
};

sub help_synopsis {
    return <<EOS
gmt cufflinks cuffmerge --input-gtf-paths=?  ...
EOS
}

sub help_brief {
    "A wrapper around the cuffmerge command";
}

sub help_detail {
    return <<EOF

Cufflinks includes a script called cuffmerge that you can use to merge together several Cufflinks assemblies. It handles also handles running Cuffcompare for you, and automatically filters a number of transfrags that are probably artfifacts. If you have a reference GTF file available, you can provide it to the script in order to gracefully merge novel isoforms and known isoforms and maximize overall assembly quality. The main purpose of this script is to make it easier to make an assembly GTF file suitable for use with Cuffdiff

More information about cuffcompare can be found at http://cufflinks.cbcb.umd.edu/.
EOF
}

sub execute {
    my $self = shift;

    my $cmd = $self->cuffmerge_path;

    my @input_files;
    my @output_files;

    my $output_directory = $self->output_directory;
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }
    $cmd .= ' -o '. $output_directory;
    
    my $input_gtfs = $self->input_gtf_paths;
    my $input_gtf_string;
    my @input_gtfs;
    if (ref($input_gtfs) eq 'ARRAY') {
        @input_gtfs = @{$input_gtfs};
        push @input_files, @input_gtfs;
        $input_gtf_string = join(' ', @{$input_gtfs});
    } else {
        @input_gtfs = split(',',$input_gtfs);
        $input_gtf_string = join(' ', @input_gtfs);
        push @input_files, @input_gtfs;
    }

    my $transcript_gtf_fof = $output_directory .'/cufflinks_GTF_list.txt';
    my $transcript_gtf_fof_fh = Genome::Sys->open_file_for_writing($transcript_gtf_fof);
    for my $input_file (@input_files) {
        unless (-e $input_file) {
            $self->error_message('Failed to find cufflinks transcript GTF file: '. $input_file);
            return;
        }
        print $transcript_gtf_fof_fh $input_file ."\n";
    }
    $transcript_gtf_fof_fh->close;
    push @input_files, $transcript_gtf_fof;
    

    if ($self->reference_gtf_path) {
        $cmd .= ' --ref-gtf '. $self->reference_gtf_path;
        push @input_files, $self->reference_gtf_path;
    }
    if ($self->reference_fasta_path) {
        $cmd .= ' --ref-sequence '. $self->reference_fasta_path;
        push @input_files, $self->reference_fasta_path;
    }

    if ($self->keep_tmp) {
        $cmd .= ' --keep-tmp';
    }

    # TODO: Define expected output files
    
    if ($self->num_threads) {
        $cmd .= ' --num-threads '. $self->num_threads;
    }
    
    if ($self->min_isoform_fraction) {
        $cmd .= ' --min-isoform-fraction '. $self->min_isoform_fraction;
    }

    $cmd .= ' '. $transcript_gtf_fof;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        # TODO: uncomment once output_files are defined
        #output_files => \@output_files,
    );
    return 1;
}

1;
