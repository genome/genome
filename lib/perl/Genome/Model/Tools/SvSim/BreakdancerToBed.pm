package Genome::Model::Tools::SvSim::BreakdancerToBed;

use Genome;
use strict;
use warnings;

class Genome::Model::Tools::SvSim::BreakdancerToBed {
    is => "Command::V2",
    has_input => [
        input_file => {
            is => "Text",
            doc => "SV calls in breakdancer format",
        },

        insertion_output_bed => {
            is => "Text",
            doc => "Output bed file of insertions",
            is_optional => 1,
        },

        deletion_output_bed => {
            is => "Text",
            doc => "Output bed file of deletions",
            is_optional => 1,
        },

        inversion_output_bed => {
            is => "Text",
            doc => "Output bed file of inversions",
            is_optional => 1,
        },

        ctx_output_bedpe => {
            is => "Text",
            doc => "Output bedpe file of CTX",
            is_optional => 1,
        },

        include_score => {
            is => "Boolean",
            is_optional => 1,
            default_value => 0,
        }
    ],
};

sub _process_bed {
    my ($self, $fh, @fields) = @_;

    my $seq = $fields[0];
    my $size = abs($fields[7]);
    my $start = $fields[1] - 1;
    my $stop = $fields[4];
    my $score = $fields[8];

    my @out = ($seq, $start, $stop, $size);
    push @out, $score if $self->include_score;
    $fh->print(join("\t", @out) . "\n");
}

sub _process_bedpe {
    my ($self, $fh, @fields) = @_;

    my $size = abs($fields[7]);

    my $seq1 = $fields[0];
    my $start1 = $fields[1] - 1;
    my $stop1 = $start1 + $size;

    my $seq2 = $fields[3];
    my $start2 = $fields[4] - 1;
    my $stop2 = $start2 + $size;


    $fh->print(join("\t",
        $seq1, $start1, $stop1,
        $seq2, $start2, $stop2,
        $size) . "\n");

}

sub execute {
    my $self = shift;


    my %paths = (
        INS => $self->insertion_output_bed,
        INV => $self->inversion_output_bed,
        DEL => $self->deletion_output_bed,
        CTX => $self->ctx_output_bedpe,
        );

    my %fhs = map {
        $_ => Genome::Sys->open_file_for_writing($paths{$_}) if defined $paths{$_}
        } keys %paths;

    die "No output files specified" unless scalar(keys(%fhs)) > 0;

    my $fh = Genome::Sys->open_file_for_reading($self->input_file);
    while (my $line = $fh->getline) {
        next if ($line =~ /^#/);
        chomp $line;
        my @fields = split(/\t/, $line);
        my $type = $fields[6];
        if ($type eq "CTX" and exists $fhs{CTX}) {
            $self->_process_bedpe($fhs{$type}, @fields);
        }
        else {
            $self->_process_bed($fhs{$type}, @fields) if exists $fhs{$type};
        }
    }

    print Data::Dumper::Dumper(\%paths);
    print Data::Dumper::Dumper(\%fhs);
    map {$_->close()} values %fhs;

    for my $p (grep {-s $_} values(%paths)) {
        print "Sorting $p\n";
        my $tmp = Genome::Sys->create_temp_file_path;
        my $cmd = Genome::Model::Tools::Joinx::Sort->create(
            input_files => [$p],
            output_file => $tmp
            );
        $cmd->execute();
        unlink($p);
        Genome::Sys->move_file($tmp, $p);
    }

    return 1;
}

1;
