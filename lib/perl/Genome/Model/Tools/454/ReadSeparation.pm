package Genome::Model::Tools::454::ReadSeparation;

use strict;
use warnings;

use Genome;

use File::Basename;


class Genome::Model::Tools::454::ReadSeparation {
    is => ['Genome::Model::Tools::454'],
    has => [
            sff_file => {
                         is => 'String',
                         is_input => 1,
                         doc => 'The sff format sequence file to separate reads by primer',
                     },
            _primer_fasta => {
                              is => 'String',
                              doc => 'A private variable to store the tmp path of the created primer db',
                          },
        ],
    has_many => [
                 primers => {
                             is => 'String',
                             is_input => 1,
                             doc => 'From the command line, a comma separated list of primer names(default_value=M13,MID1).  As an object, an array ref of primer names.',
                             is_optional => 1,
                         },
             ],
};

sub help_brief {
    "tool to separate reads for 454"
}

sub help_detail {
    return <<EOS
Given an sff file this tool will separate the reads based on their alignment to a primer sequence.
Currently, the primer fasta must be found here:
/gscmnt/sata180/info/medseq/biodb/shared/Vector_sequence/\$PRIMER-primers.fasta
First, this tool isolates the first 20 bp of the reads.
Second, performs a cross_match alignment of the first 20 bp to the primer fasta.
Finally, separates the whole reads based on their alignment to the primer fasta.
The output is an sff format file with the same root name as the input sff and the primer name such as, \$INPUT.\$PRIMER.sff
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    my @default_primers = qw(M13 MID1);
    unless ($self->primers) {
        $self->warning_message('Using default primers: '. join(' ', @default_primers));
        $self->primers(\@default_primers);
    }

    my $primer_dir = '/gscmnt/sata180/info/medseq/biodb/shared/Vector_sequence';

    my @primer_fastas;
    my $db_name = $self->_tmp_dir .'/';

    for my $primer ($self->primers) {
        $db_name .= $primer .'-';
        my $primer_fasta = $primer_dir .'/'. $primer .'-primers.fasta';
        unless (-s $primer_fasta) {
            die("Primer fasta file '$primer_fasta' is missing or zero size!");
        }
        push @primer_fastas, $primer_fasta;
    }
    $db_name .= 'primers.fasta';
    my $cmd = 'cat ';
    for my $primer_fasta (@primer_fastas) {
        $cmd .= $primer_fasta .' ';
    }
    $cmd .= '> '. $db_name;
    my $rv = system($cmd);
    unless ($rv == 0) {
        die("Failed to execute system command: '$cmd'");
    }

    $self->_primer_fasta($db_name);

    return $self;
}


sub execute {
    my $self = shift;
$DB::single = 1;
    my $sff_file_dirname = dirname($self->sff_file);
    my $sff_file_basename = basename($self->sff_file);
    $sff_file_basename =~ s/\.sff$//;

    my $out_sff_file = $self->_tmp_dir .'/'. $sff_file_basename .'_20bp.sff';
    my $cross_match_file = $self->_tmp_dir .'/'. $sff_file_basename .'_20bp.cm';

    my $isolate_primer = Genome::Model::Tools::454::IsolatePrimerTag->create(
                                                                             in_sff_file => $self->sff_file,
                                                                             out_sff_file => $out_sff_file,
									     version => $self->version,
									     version_subdirectory => $self->version_subdirectory,
                                                                         );
    unless ($isolate_primer->execute) {
        $self->error_message('Failed to execute '. $isolate_primer->command_name);
        return;
    }
    my $cross_match = Genome::Model::Tools::454::CrossMatchPrimerTag->create(
                                                                             cross_match_file => $cross_match_file,
                                                                             sff_file => $out_sff_file,
                                                                             primer_fasta => $self->_primer_fasta,
									     version => $self->version,
									     version_subdirectory => $self->version_subdirectory,
                                                                         );
    unless ($cross_match->execute) {
        $self->error_message('Failed to execute '. $cross_match->command_name);
        return;
    }

    my $separate_reads = Genome::Model::Tools::454::SeparateReadsWithCrossMatchAlignment->create(
                                                                                                 cross_match_file => $cross_match_file,
                                                                                                 sff_file => $self->sff_file,
												 version => $self->version,
												 version_subdirectory => $self->version_subdirectory,
                                                                                             );
    unless ($separate_reads->execute) {
        $self->error_message('Failed to execute '. $separate_reads->command_name);
        return;
    }
    return 1;
}


1;

