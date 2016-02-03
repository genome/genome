package Genome::Model::Tools::Picard::FilterSamReads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::FilterSamReads {
    is => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
          is  => 'String',
          doc => 'The SAM or BAM file that will be filtered.',
          picard_param_name => 'INPUT',
      },
      output_file => {
          is  => 'String',
          doc => 'The resulting SAM/BAM file',
          picard_param_name => 'OUTPUT',
      },
      filter => {
          is           => 'String',
          doc          => 'Filter to pass to picard, if selected filter uses a Read List, you must specify one',
          valid_values => ['includeAligned', 'excludeAligned', 'includeReadList', 'excludeReadList'],
          picard_param_name => 'FILTER',
      },
      read_list_file => {
          is          => 'String',
          doc         => 'File containing reads that will be included or excluded from the output',
          is_optional => 1,
          picard_param_name => 'READ_LIST_FILE',
      },
      sort_order => {
          is           => 'String',
          doc          => 'Specify only if you want a different sort order than the input file',
          is_optional  => 1,
          valid_values => ['unsorted', 'queryname', 'coordinate'],
          picard_param_name => 'SORT_ORDER',
      },
      write_reads_files => {
          is            => 'Boolean',
          doc           => 'Create .reads files during run',
          is_optional   => 1,
          default_value => '1',
          picard_param_name => 'WRITE_READS_FILES',
      },
    ],
};

sub help_brief {
    'Tool to filter reads out of a SAM/BAM file';
}

sub help_detail {
   return <<EOS
    Tool to filter reads out of a SAM/BAM file. For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#FilterSamReads
EOS
}

sub minimum_version_required { '1.77'; }
sub _jar_name { 'FilterSamReads.jar'; }
sub _java_class {
    return qw(picard sam FilterSamReads);
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files               => [$self->input_file],
        output_files              => [$self->output_file],
        skip_if_output_is_present => 0,
    );
}


1;
