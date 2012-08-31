package Genome::Model::Tools::Blat::PslToLayers;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Blat::PslToLayers {
    is => ['Genome::Model::Tools::Blat','Genome::Sys'],
    has => [
            psl_file => {
                         is => 'Text',
                     },
            layers_file => {
                            is => 'Text',
                        }
        ],
    has_optional => [
                     randomize => {
                                   is => 'Boolean',
                                   doc => 'run shuf to randomize the reads',
                                   default_value => 0,
                               },
                 ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless ($self->validate_file_for_reading($self->psl_file)) {
        $self->error_message('Failed to validata map file for reading '. $self->psl_file);
        return;
    }
    unless ($self->validate_file_for_writing($self->layers_file)) {
        $self->error_message('Failed to validate layers file for writing '. $self->layers_file);
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::PSL::Reader->create( file => $self->psl_file);
    unless ($reader) {
        $self->error_message('Could not create a reader for file '. $self->psl_file);
        return;
    }
    my $layers_fh = $self->open_file_for_writing($self->layers_file);
    unless ($layers_fh) {
        $self->error_message('Failed to open layers file for writing '. $self->layers_file);
        return;
    }
    while (my $record = $reader->next) {
        #blat output is zero based coordinate system so output start and stop with +1
        my $read_name = $$record{qName};
        my $ref_name = $$record{tName};
        my $block_count = $$record{blockCount};
        my @block_sizes = split(",",$$record{blockSizes});
        my @start_sites = split(",",$$record{tStarts});
        my @seq_blocks = split(",",$$record{tSeq});
        unless ($block_count == scalar(@block_sizes)) {
            die ('block count '. $block_count .' does not match number of block sizes '. scalar(@block_sizes));
        }
        for (my $i = 0; $i < $block_count; $i++) {
            my $start = $start_sites[$i] + 1;
            my $size = $block_sizes[$i];
            # This seems logical but the tEnd does not add up to the $stop for the last block???
            # but then the blockSize plus tStart overlaps other alignments?? shouldn't there be gaps??
            my $stop = ($start + $size) - 1;
            my $seq = $seq_blocks[$i];
            print $layers_fh $read_name ."\t". $start ."\t". $stop  ."\t". $ref_name ."\t". uc($seq) ."\n";
        }
    }
    $reader->close;
    $layers_fh->close;
    if ($self->randomize) {
        my $basename = File::Basename::basename($self->psl_file);
        my $random_file = $self->create_temp_file_path($basename .'.randomized');
        my $cmd = '/gsc/pkg/coreutils/coreutils-6.10-64/shuf -o '. $random_file .' '. $self->layers_file;
        $self->shellcmd(
                        cmd => $cmd,
                        input_files => [$self->layers_file],
                        output_files => [$random_file],
                    );
        unless (unlink $self->layers_file) {
            $self->error_message('Failed to remove un-randomized layers file'. $self->layers_file .": $!");
            die($self->error_message);
        }
        unless ($self->copy_file($random_file,$self->layers_file)) {
            $self->error_message('Failed to copy randomized file '. $random_file .' to output layers file '. $self->layers_file .":  $!");
            die($self->error_message);
        }
    }
    return 1;
}



1;

