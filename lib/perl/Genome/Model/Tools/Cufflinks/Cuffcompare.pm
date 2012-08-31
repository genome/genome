package Genome::Model::Tools::Cufflinks::Cuffcompare;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cufflinks::Cuffcompare {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        input_gtf_paths => {
            doc => 'A list of GTF paths used as the input(s)',
        },
        output_prefix => {
            is => 'Text',
            doc => 'The prefix used to generate the path for the combined GTF and tracking files',
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
        include_contained => {
            is => 'Boolean',
            doc => 'Synonym for -C command line flag: include the "contained" transcripts in the .combined.gtf file',
            default_value => 0,
            is_optional => 1,
        },
        generic_gtf_input => {
            is => 'Boolean',
            doc => 'Synonym for -G command line flag: generic GFF input file(s) (do not assume Cufflinks GTF)',
            default_value => 0,
            is_optional => 1,
        },
        generate_tracking_files => {
            is => 'Boolean',
            doc => 'Synonym for -T command line flag: generate .tmap and .refmap files for each input file',
            default_value => 1,
            is_optional => 1,
        },
        verbose => {
            is => 'Boolean',
            doc => 'Synonym for -V command line flag: verbose processing mode (showing all GFF parsing warnings)',
            default_value => 0,
            is_optional => 1,
        },
        transcript_prefix => {
            is => 'Text',
            doc => 'Synonym for -p command line option: the name prefix to use for consensus transcripts in the <outprefix>.combined.gtf file (default: \'TCONS\')',
            is_optional => 1,
        },
        max_distance => {
            is => 'Integer',
            doc => 'Synonym for -d command line option: max distance (range) for grouping transcript start sites (100)',
            is_optional => 1,
        },
        discard_single_exon => {
            is => 'Text',
            doc => 'discard (ignore) transcripts (-N) or transcripts and transfrags(-M) with single exons',
            valid_values => ['none','transcripts', 'transfrags_and_transcripts'],
            default_value => 'none',
            is_optional => 1,
        },
        reduce_transcripts => {
            is => 'Text',
            doc => 'Synonm for -R command line flag:  for -r option, reduce the set of reference transcripts to only those found to overlap any of the input loci',
            default_value => 0,
            is_optional => 1,
        },
    ],
    has_output_optional => [
        combined_gtf_path => {},
        loci_path => {},
        stats_path => {},
        tracking_path => {},
        refmap_paths => {},
        tmap_paths => {},
    ],
};

sub help_synopsis {
    return <<EOS
gmt cufflinks cuffcompare --input-gtf-paths=? ...
EOS
}

sub help_brief {
    "A wrapper around the cuffcompare command";
}

sub help_detail {
    return <<EOF
Cuffcompare provides classification, reference annotation mapping and various statistics for Cufflinks transfrags.  Cuffcompare clusters and tracks transfrags across multiple samples, writing matching transcripts (intron chains) into <outprefix>.tracking, and a GTF file <outprefix>.combined.gtf containing a nonredundant set of transcripts across all input files (with a single representative transfrag chosen for each clique of matching transfrags across samples).

More information about cuffcompare can be found at http://cufflinks.cbcb.umd.edu/.
EOF
}

sub execute {
    my $self = shift;

    my $cmd = $self->cuffcompare_path;

    my @input_files;
    my @output_files;

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

    if ($self->reference_gtf_path) {
        $cmd .= ' -r '. $self->reference_gtf_path;
        push @input_files, $self->reference_gtf_path;
    }

    if ($self->reduce_transcripts) {
        $cmd .= ' -R';
    }
    my $output_prefix = 'cuffcmp';
    if ($self->output_prefix) {
        $cmd .= ' -o '. $self->output_prefix;
        $output_prefix = $self->output_prefix;
    }

    $self->combined_gtf_path(File::Spec->rel2abs($output_prefix .'.combined.gtf'));
    $self->loci_path(File::Spec->rel2abs($output_prefix .'.loci'));
    $self->stats_path(File::Spec->rel2abs($output_prefix .'.stats'));
    $self->tracking_path(File::Spec->rel2abs($output_prefix .'.tracking'));

    push @output_files, $self->combined_gtf_path;
    push @output_files, $self->loci_path;
    push @output_files, $self->stats_path;
    push @output_files, $self->tracking_path;

    my @refmap_files;
    my @tmap_files;
    unless ($self->generate_tracking_files) {
        $cmd .= ' -T';
    } else {
        for my $input_gtf (@input_gtfs) {
            my ($output_basename,$output_dirname) = File::Basename::fileparse($output_prefix);
            my ($input_basename,$input_dirname,$input_suffix) = File::Basename::fileparse($input_gtf,qw/\.gtf/);
            my $map_prefix = $input_dirname . $output_basename .'.'. $input_basename . $input_suffix;
            push @refmap_files, $map_prefix .'.refmap';
            push @tmap_files, $map_prefix .'.tmap';
        }
    }
    if (@refmap_files && @tmap_files) {
        $self->refmap_paths(\@refmap_files);
        $self->tmap_paths(\@tmap_files);
        push @output_files, @refmap_files;
        push @output_files, @tmap_files;
    }
    if ($self->verbose) {
        $cmd .= ' -V';
    }

    if ($self->include_contained) {
        $cmd .= ' -C';
    }

    if ($self->generic_gtf_input) {
        $cmd .= ' -G';
    }

    if ($self->reference_fasta_path) {
        # TODO: should we test if the fasta is a file or directory?
        push @input_files, $self->reference_fasta_path;
        $cmd .= ' -s '. $self->reference_fasta_path;
    }

    if ($self->transcript_prefix) {
        $cmd .= ' -p '. $self->transcript_prefix;
    }

    if (defined($self->max_distance)) {
        $cmd .= ' -d '. $self->max_distance;
    }

    if ($self->discard_single_exon eq 'transfrags_and_transcripts') {
        $cmd .= ' -M';
    } elsif ($self->discard_single_exon eq 'transcripts') {
        $cmd .= ' -N';
    } elsif ($self->discard_single_exon eq 'none') {
        # Nothing to do
    } else {
        die('Invalid option for discard_single_exon: '. $self->discard_single_exon);
    }

    $cmd .= ' '. $input_gtf_string;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
    );
    return 1;
}

1;
