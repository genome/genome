package Genome::Model::Tools::SmrtAnalysis::CsvMerge;

class Genome::Model::Tools::SmrtAnalysis::CsvMerge {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        input_csv_files => {
        },
        output_csv_file => {
            is_output => 1,
        },
    ],
    has_optional_input => [
        remove_originals => {
            default_value => 1,
        },
        separator => {
            default_value => ',',
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $input_files = delete($params{input_csv_files});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->input_csv_files($input_files);
    return $self;
}

sub execute {
    my $self = shift;
    my @input_csv_files;
    if (ref($self->input_csv_files) eq 'ARRAY') {
        @input_csv_files = @{$self->input_csv_files};
    } else {
        @input_csv_files = split(',',$self->input_csv_files);
    }
    my $writer;
    for my $file (@input_csv_files) {
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            input => $file,
            separator => $self->separator,
        );
        unless ($writer) {
            $writer = Genome::Utility::IO::SeparatedValueWriter->create(
                output => $self->output_csv_file,
                headers => $reader->headers,
            );
        }
        while (my $data = $reader->next) {
            $writer->write_one($data);
        }
    }
    if ($self->remove_originals) {
        for my $file (@input_csv_files) {
            unlink($file) || die('Failed to remove original file '. $file);
        }
    }
    return 1;
}


1;
