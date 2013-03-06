
package Genome::Model::Tools::Picard::FastqToSam;

# originally copied from Genome::Model::Tools::Picard::SortSam

use strict;
use warnings FATAL => 'all';

use Genome;

my $DEFAULT_SORT_ORDER = 'queryname';
my $DEFAULT_MAX_RECORDS_IN_RAM = 500000;


class Genome::Model::Tools::Picard::FastqToSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        fastq => {
            is  => 'String',
            doc => 'First input FASTQ file.',
        },
        fastq2 => {
            is          => 'String',
            doc         => 'Second input FASTQ file (if paired end).',
            is_optional => 1,
        },
        output => {
            is  => 'String',
            doc => 'The resulting sorted SAM/BAM file.  File type is determined by suffix.',
        },
        platform_unit => {
            is  => 'String',
            doc => 'The platform unit (often run_barcode.lane) to insert into the read group header.',
            is_optional => 1,
        },
        quality_format => {
            is  => 'String',
            doc => 'A value describing how the quality values are encoded in the fastq files.',
            valid_values => [qw(Solexa Illumina Standard)],
        },
        platform => {
            is  => 'String',
            doc => 'The platform type (e.g. illumina, solid) to insert into the read group header.',
            is_optional => 1,
        },
        sort_order => {
            is  => 'String',
            doc => 'The sort order for the output sam/bam file.  default_value=' . $DEFAULT_SORT_ORDER,
            valid_values  => [ 'coordinate', 'unsorted', 'queryname' ],
            default_value => $DEFAULT_SORT_ORDER,
            is_optional   => 1,
        },
        sample_name => {
            is  => 'String',
            doc => 'Sample name to insert into the read group header.',
        },
        library_name => {
            is  => 'String',
            doc => 'The library name to place into the LB attribute in the read group header.',
            is_optional => 1,
        },
        read_group_name => {
            is  => 'String',
            doc => 'Read group name Default value: A.',
            is_optional => 1,
        },
        sequencing_center => {
            is  => 'String',
            doc => 'The sequencing center from which the data originated.',
            is_optional => 1,
        },
        predicted_insert_size => {
            is  => 'Integer',
            doc => 'Predicted median insert size, to insert into the read group header.',
            is_optional => 1,
        },
        comment => {
            is  => 'String',
            #is_many => 1,
            doc => 'Comment(s) to inlcude in the merged output file\'s header.',
            is_optional => 1,
        },
        description => {
            is  => 'String',
            doc => 'Inserted into the read group header',
            is_optional => 1,
        },
        run_date => {
            is  => 'iso8601Date',
            doc => 'Date the run was produced, to insert into the read group header.',
            is_optional => 1,
        },
        min_q => {
            is  => 'Integer',
            doc => 'Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.',
            is_optional => 1,
        },
        max_q => {
            is  => 'Integer',
            doc => 'Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to create SAM/BAM file from FASTQ using Picard';
}

sub help_detail {
    return <<EOS
    Tool to create SAM/BAM file from FASTQ using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#FastqToSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path . '/FastqToSam.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }

    my $args = join( ' ',
        map {
            my $value = $self->$_;
            defined($value) ? ( uc($_) . "='$value'" ) : ()
        } sort qw(
            fastq fastq2 output platform_unit quality_format platform
            sort_order sample_name library_name read_group_name max_records_in_ram
        )
    );

    my $cmd = $jar_path . " net.sf.picard.sam.FastqToSam $args";
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [ $self->fastq, ( $self->fastq2 ? $self->fastq2 : () ) ],
        output_files => [ $self->output ],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
__END__

