package Genome::Model::Tools::Picard::FilterSamReads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::FilterSamReads {
    is => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
          is  => 'String',
          doc => 'The SAM or BAM file that will be filtered.',
      },
      output_file => {
          is  => 'String',
          doc => 'The resulting SAM/BAM file'
      },
      filter => {
          is           => 'String',
          doc          => 'Filter to pass to picard, if selected filter uses a Read List, you must specify one',
          valid_values => ['includeAligned', 'excludeAligned', 'includeReadList', 'excludeReadList'],
      },
      read_list_file => {
          is          => 'String',
          doc         => 'File containing reads that will be included or excluded from the output',
          is_optional => 1,
      },
      sort_order => {
          is           => 'String',
          doc          => 'Specify only if you want a different sort order than the input file',
          is_optional  => 1,
          valid_values => ['unsorted', 'queryname', 'coordinate'],
      },
      write_reads_file => {
          is            => 'String',
          doc           => 'Create .reads file during run',
          is_optional   => 1,
          default_value => 'true',
          valid_values  => ['true', 'false']
      },
    ],
};

sub help_brief {
    'Tool to filter reads out of a SAM/BAM file';
}

sub help_detail {
   return <<EOS
    Tool to filter reads out of a SAM/BAM file. For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#FilterSamReads
EOS
}

sub execute {
    my $self = shift;

    my $input_file = $self->input_file;
    my $output_file = $self->output_file;

    my $jar_path = $self->picard_path . '/FilterSamReads.jar';
    unless ( -e $jar_path ) {
        die( "Failed to find $jar_path!" .
          "This command may not be available in version " . $self->use_version );
    }

    my $cmd = $self->_generate_cmd_string( $jar_path );

    $self->run_java_vm(
        cmd                       => $cmd,
        input_files               => [$input_file],
        output_files              => [$output_file],
        skip_if_output_is_present => 0,
    );

    return 1;
}

sub _generate_cmd_string {
    my $self = shift;
    my $jar_path = shift;

    my $input_file = $self->input_file;
    my $output_file = $self->output_file;

    #build command, including required args by default
    my $cmd = $jar_path . " net.sf.picard.sam.FilterSamReads O=$output_file" .
      " I=$input_file FILTER=" . $self->filter .
      " WRITE_READS_FILES=" . $self->write_reads_file;

    #optional params
    if ( defined( $self->read_list_file ) ) {
        $cmd .= " READ_LIST_FILES=" . $self->read_list_file;
    }
    if ( defined( $self->sort_order ) ) {
        $cmd .= " SORT_ORDER=" . $self->sort_order;
    }

    return $cmd;
}

1;
