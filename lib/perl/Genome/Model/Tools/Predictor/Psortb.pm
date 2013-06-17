package Genome::Model::Tools::Predictor::Psortb;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Predictor::Psortb {
    is => 'Genome::Model::Tools::Predictor::Base',
    doc => 'execute psortb',
};

sub requires_chunking {
    return 1;
}

sub run_predictor {
    my $self = shift;

    my $gram_stain = $self->gram_stain;
    unless (defined $gram_stain) {
        warn "No value provided for gram stain!";
    } elsif ($gram_stain eq 'positive') {
        $gram_stain = '-p';
    }
    elsif ($gram_stain eq 'negative') {
        $gram_stain = '-n';
    }
    else {
        die "Unexpected value $gram_stain for gram stain!";
    }
    
    # TODO Need to get path for psort-b based on version
    my $cmd = join(' ', 'psort-b', $gram_stain, '-o', 'terse', $self->input_fasta_file,
        ">", $self->raw_output_path, "2>", $self->debug_output_path);
    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    unless ($rv) {
        die "Could not execute psortb!";
    }

    return 1;
}

sub parse_output {
    my $self = shift;

    my $output_fh = IO::File->new($self->raw_output_path, 'r');
    unless ($output_fh) {
        die "Could not get file handle for raw output file " . $self->raw_output_path;
    }

    my @features;
    while (my $line = $output_fh->getline) {
        chomp $line;

        if ($line =~ /^SeqID/) {
            next;
        }
        
        my ($gene, $class, $score) = split(/\t/,$line);
        $gene =~ s/\s$//; # psort-b has been appending a space to this...
        
        my $feature = Bio::SeqFeature::Generic->new(
            -display_name => $gene,
        );
        $feature->add_tag_value('psort_localization', $class);
        $feature->add_tag_value('psort_score', $score);
        push @features, $feature;
    }

    $self->unfiltered_bio_seq_features(\@features);
    return 1;
}

sub filter_results {
    my $self = shift;
    my @features;
    for my $feature ($self->unfiltered_bio_seq_features) {
        my ($class) = $feature->get_tag_values('psort_localization');
        if ($class =~ /unknown/i) {
            next;
        }
        push @features, $feature;
    }
    $self->bio_seq_features(\@features);
    return 1;
}

sub create_ace_file {
    my $self = shift;
    my $raw_output = $self->raw_output_path;
    my $ace_file_path = $self->ace_file_path;

    my $raw_fh = Genome::Sys->open_file_for_reading($raw_output);
    my $ace_fh = Genome::Sys->open_file_for_writing($ace_file_path);

    readline $raw_fh; #skip header line
    while(my $line = <$raw_fh>) {
        #chomp $line;  #keep newline to print out to ace file
        my @f = split("\t", $line);
        next if $f[1] eq 'Unknown';

        #need extra newline at end to separate blocks
        print $ace_fh "Sequence " . $f[0] . "\nPSORT_B " . $f[1] . ' ' . $f[2] . "\n";
    }

    $raw_fh->close;
    $ace_fh->close;

    return 1;
}

# TODO Need to implement this instead of just relying on 'psort-b' working on command line
sub tool_path_for_version {
    my $self = shift;
    return '1.0';
}

sub raw_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'psortb.output');
}

sub debug_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'psortb.debug');
}

sub dump_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'psortb.dump');
}

sub ace_file_path {
    my $self = shift;
    return join('/', $self->output_directory, 'psortb.ace');
}

1;

