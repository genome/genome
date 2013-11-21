package Genome::Model::GenePrediction::Command::Pap::PsortB;

use strict;
use warnings;

use PAP;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use English;
use IO::File;
use IPC::Run;

class Genome::Model::GenePrediction::Command::Pap::PsortB {
    is  => 'Command::V1',
    has => [
        psortb_archive_dir => {
            is => 'Path',
            doc => 'Raw psortb output is placed in this directory',
            is_input => 1,
        },
        fasta_file => { 
            is  => 'Path',
            doc => 'Path to fasta file',
            is_input => 1,
        },
        gram_stain => {
            is  => 'Text',
            doc => 'gram stain (positive/negative)',
            valid_values => ['positive', 'negative'],
            is_input => 1,
        },
    ],
    has_optional => [
        bio_seq_feature => { 
            is => 'ARRAY', 
            doc => 'array of Bio::Seq::Feature', 
            is_output => 1,
        },
    ],
    has_param => [
        lsf_resource => { 
            default_value => "-q $ENV{GENOME_LSF_QUEUE_SHORT} -R 'select[type==LINUX64] rusage[tmp=100]'", 
        },
    ],
};

sub help_brief {
    "Run psort-b";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;
    
    my $fasta_file  = $self->fasta_file();
    my $gram_stain  = $self->gram_stain();

    if ($gram_stain eq 'positive') {
        $gram_stain = '-p';
    }
    elsif ($gram_stain eq 'negative') {
        $gram_stain = '-n';
    }

    my $temp_fh = File::Temp->new(
        TEMPLATE => 'psortb_raw_output_XXXXXX',
        SUFFIX => '.out',
        DIR => $self->psortb_archive_dir,
        CLEANUP => 0,
        UNLINK => 0,
    );
    my $output_file = $temp_fh->filename();
    $temp_fh->close();
    
    my @psortb_command = (
                          'psort-b',
                          $gram_stain, 
                          '-o',
                          'terse',
                          $fasta_file,
                      );
    
    my ($psortb_err);

    IPC::Run::run(
                  \@psortb_command,
                  \undef,
                  '>',
                  $output_file,
                  '2>',
                  \$psortb_err,
              ) || die "psort-b failed: $CHILD_ERROR";
    
    my $feature_ref = $self->parse_psortb_terse($output_file);
    
    $self->bio_seq_feature($feature_ref);
    return 1;
}

sub parse_psortb_terse {
    
    my ($self, $psort_fn) = (@_);
    
    
    my @features = ( );

    my $psort_fh = IO::File->new();
    $psort_fh->open("$psort_fn") or die "Can't open '$psort_fn': $OS_ERROR";

    LINE: while (my $line = <$psort_fh>) {

        chomp $line;
        
        if ($line =~ /^SeqID/) {
            next;
        }
        
        my ($gene, $class, $score) = split(/\t/,$line);
        $gene =~ s/\s$//; # psort-b has been appending a space to this...

        if ($class =~ /unknown/i) { next LINE; }
        
        my $feature = Bio::SeqFeature::Generic->new(
                                                    -display_name => $gene,
                                                );

        $feature->add_tag_value('psort_localization', $class);
        $feature->add_tag_value('psort_score', $score);
        
        push @features, $feature;
        
    }

    return \@features;
    
}

1;
