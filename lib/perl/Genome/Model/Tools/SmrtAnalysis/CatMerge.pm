package Genome::Model::Tools::SmrtAnalysis::CatMerge;

class Genome::Model::Tools::SmrtAnalysis::CatMerge {
    is  => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        input_files => {
        },
        output_file => {
            is_output => 1,
        },
    ],
    has_optional_input => [
        remove_originals => {
            default_value => 1,
        },
    ],

};

sub create {
    my $class = shift;
    my %params = @_;
    my $input_files = delete($params{input_files});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->input_files($input_files);
    return $self;
}

sub execute {
    my $self = shift;
    my @input_files;
    if (ref($self->input_files) eq 'ARRAY') {
        @input_files = @{$self->input_files};
    } else {
        @input_files = split(',',$self->input_files);
    }
    # Some of the PacBio Fofn files are missing newline characters so regular cat won't work
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    for my $input_file (@input_files) {
        my $input_fh = Genome::Sys->open_file_for_reading($input_file);
        while ( my $line = $input_fh->getline) {
            chomp($line);
            print $output_fh $line ."\n";
        }
        $input_fh->close;
    }
    $output_fh->close;
    if ($self->remove_originals) {
        for my $file (@input_files) {
            unlink($file) || die('Failed to remove original file '. $file);
        }
    }
    return 1;
}


1;
