package Genome::Model::Tools::454::IsolatePrimerTag;

use strict;
use warnings;

use Genome;

use File::Basename;

class Genome::Model::Tools::454::IsolatePrimerTag {
    is => ['Genome::Model::Tools::454'],
    has => [
            in_sff_file => {
                            doc => 'The sff file to operate',
                            is => 'String',
                        },
            out_sff_file => {
                             is => 'String',
                             doc => 'The output file path',
                         },
            primer_length => {
                              is => 'Integer',
                              default_value => 20,
                              doc => 'The length of the expected primer to isolate(default=20)',
                          },
        ],
};

sub help_brief {
    "isolate the the sequence primer from a set of reads"
}

sub help_detail {
    return <<EOS
create a new sff file with first n(default=20) base pair were the expected primer should be found
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('This genome-model tool '. $self->command_name .' will only run on 64-bit');
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    unless (-e $self->in_sff_file) {
        die ('Failed to find file '. $self->in_sff_file ."\n");
    }

    my $basename = basename($self->in_sff_file);
    my $tmp_trim_file = $self->_tmp_dir . '/'. $basename .'.trim';

    #for testing
    #my $tmp_trim_file = $self->in_sff_file .'.trim';

    my $writer = Genome::Utility::454TrimFile::Writer->create(file => $tmp_trim_file);
    unless ($writer) {
        die ('Failed to create writer to '. $self->out_sff_file ."\n");
    }
    my $cmd = $self->bin_path .'/sffinfo -a '. $self->in_sff_file;
    open(ACC,"$cmd |");
    while (my $line = <ACC>) {
        chomp($line);
        my %record = (
                      accession => $line,
                      start => 1,
                      end => $self->primer_length,
                  );
        $writer->write_record(\%record);
    }
    $writer->close;

    my $read_trimmer = Genome::Model::Tools::454::Sfffile->create(
                                                                  in_sff_files => [$self->in_sff_file],
                                                                  out_sff_file => $self->out_sff_file,
                                                                  params => '-tr '. $tmp_trim_file,
								  version => $self->version,
								  version_subdirectory => $self->version_subdirectory,
                                                              );
    unless ($read_trimmer->execute) {
        die ('Failed to trim file '. $self->in_sff_file);
    }
    return 1;
}

1;


