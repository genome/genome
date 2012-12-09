package Genome::Model::Event::Build::AmpliconAssembly::PrepareInstrumentData;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Copy;

class Genome::Model::Event::Build::AmpliconAssembly::PrepareInstrumentData {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    my $amplicons = $self->build->get_amplicons
        or return;
    
    my $method = '_prepare_instrument_data_for_phred_phrap';
    
    for my $amplicon ( @$amplicons ) {
        $self->$method($amplicon)
            or return;
    }

    return 1;
}

sub _prepare_instrument_data_for_phred_phrap {
    my ($self, $amplicon) = @_;

    $amplicon->create_scfs_file
        or return;

    my $scf2phd = Genome::Model::Tools::PhredPhrap::ScfToPhd->create(
        chromat_dir => $self->build->chromat_dir,
        phd_dir => $self->build->phd_dir,
        phd_file => $amplicon->phds_file,
        scf_file => $amplicon->scfs_file,
    );
    unless ( $scf2phd ) {
        $self->error_message("Can't create scf to phd command");
        return;
    }
    unless ( $scf2phd->execute ) {
        $self->error_message("Can't execute scf to phd command");
        return;
    } 
    
    my $phd2fasta = Genome::Model::Tools::PhredPhrap::PhdToFasta->create(
        fasta_file => $amplicon->fasta_file,
        phd_dir => $self->build->phd_dir,
        phd_file => $amplicon->phds_file,
    );
    unless ( $phd2fasta ) {
        $self->error_message("Can't create phd to fasta command");
        return;
    }
    unless ( $phd2fasta->execute ) {
        $self->error_message("Can't execute phd to fasta command");
        return;
    }
    
    # Create raw reads fasta and qual file
    my $reads_fasta = $amplicon->reads_fasta_file;
    unlink $reads_fasta if -e $reads_fasta;
    File::Copy::copy(
        $amplicon->fasta_file,
        $reads_fasta,
    );
    my $reads_qual = $amplicon->reads_qual_file;
    unlink $reads_qual if -e $reads_qual;
    File::Copy::copy(
        $amplicon->qual_file,
        $reads_qual,
    );
    
    return 1;
}

1;

#$HeadURL$
#$Id$
