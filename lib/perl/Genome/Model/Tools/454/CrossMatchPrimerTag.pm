package Genome::Model::Tools::454::CrossMatchPrimerTag;

use strict;
use warnings;

use Genome;

use File::Basename;

class Genome::Model::Tools::454::CrossMatchPrimerTag {
    is => ['Genome::Model::Tools::454'],
    has => [
            sff_file => {
                         is_input => 1,
                         is => 'String',
                         doc => 'The sff file to operate',
                     },
            primer_fasta => {
                             is_input => 1,
                             is => 'String',
                             doc => 'A fasta file of primer sequences',
                         },
            cross_match_file => {
                                 is_output => 1,
                                 is => 'String',
                                 doc => 'The output file path',
                         },
        ],
};

sub help_brief {
    "identify the sequence primer in a set of reads"
}

sub help_detail {
    return <<EOS
runs cross_match on the read set with a given input sff file of reads and a primer/adaptor fasta file
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('This genome-model tool '. $self->command_name .' will only run on 64-bit');
        return;
    }
    if (-s $self->cross_match_file) {
        $self->error_message('cross_match output file '. $self->cross_match_file .' already exists with data.');
        return;
    }
    unless (-s $self->sff_file) {
        $self->error_message('sff file '. $self->sff_file .' does not exist or contains no data.');
        return;
    }
    unless (-s $self->primer_fasta) {
        $self->error_message('primer fasta file '. $self->primer_fasta .' does not exist or contains no data.');
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    unless (-e $self->sff_file) {
        die ('Failed to find file '. $self->sff_file ."\n");
    }

    my $basename = basename($self->sff_file);
    #my $tmp_fasta_file = $self->_tmp_dir . '/'. $basename .'.fasta';

    #for testing
    my $tmp_fasta_file = $self->sff_file .'.fasta';
    my $tmp_cm_file = $tmp_fasta_file .'.cm';
    my $fasta_sffinfo = Genome::Model::Tools::454::Sffinfo->create(
                                                                   sff_file => $self->sff_file,
                                                                   output_file => $tmp_fasta_file,
                                                                   params => '-s',
								   version => $self->version,
								   version_subdirectory => $self->version_subdirectory,
                                                               );
    unless ($fasta_sffinfo->execute) {
        die('Failed to dump fasta data from '. $self->sff_file);
    }

    my $tmp_qual_file = $tmp_fasta_file .'.qual';
    my $qual_sffinfo = Genome::Model::Tools::454::Sffinfo->create(
                                                                  sff_file => $self->sff_file,
                                                                  output_file => $tmp_fasta_file .'.qual',
                                                                  params => '-q',
								  version => $self->version,
								  version_subdirectory => $self->version_subdirectory,
                                                              );
    unless ($qual_sffinfo->execute) {
        die('Failed to dump fasta data from '. $self->sff_file);
    }

    my $cross_match_tool = Genome::Model::Tools::Crossmatch::Run->create(
                                                                         query_file => $self->primer_fasta,
                                                                         subject_file => $tmp_fasta_file,
                                                                         cm_params => '-tags -masklevel 0 -gap1_only -minmatch 6 -minscore 6 -minmargin 1',
                                                                         alignment_file => $self->cross_match_file
                                                                     );
    unless ($cross_match_tool) {
        $self->error_message('Failed to create cross_match tool');
        return;
    }
    unless ($cross_match_tool->execute) {
        $self->error_message('Failed to execute command '. $cross_match_tool->command_name);
        return;
    }
    return 1;
}

1;


