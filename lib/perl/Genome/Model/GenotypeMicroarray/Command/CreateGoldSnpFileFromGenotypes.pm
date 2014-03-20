package Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpFileFromGenotypes;

use strict;
use warnings;
use Genome;

class Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpFileFromGenotypes {
    is => 'Command',
    has => [
        output_file => {
            is => 'FilePath',
            doc => 'Gold snp output file',
        },
    ],
    has_optional => [
        genotype_file => {
            is => 'FilePath',
            doc => 'input file of genotypes from one platform',
        },
        reference_sequence_build_id => {
            is => 'Number',
            doc => 'Build id of reference sequence build',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_sequence_build_id',
            doc => 'Reference sequence build',
        },
    ],
};

sub help_brief { return 'Creates a gold snp file from two genotype files' };
sub help_synopsis { return help_brief() };
sub help_detail { return help_brief() };

sub execute {
    my $self = shift;

    # Make sure we got a reference
    unless ($self->reference_sequence_build) {
        Carp::confess 'Could not resolve reference sequence build!';
    }

    # Check and open filehandles    
    my $genotype_fh = IO::File->new($self->genotype_file);
    unless($genotype_fh) {
        Carp::confess "Failed to open filehandle for: " . $self->genotype_file;
        return;
    }
    my $output_fh = IO::File->new($self->output_file,"w");
    unless($output_fh) {
        Carp::confess "Failed to open filehandle for: " . $self->output_file;
        return;
    }

    my %chromosomes;
    $chromosomes{$_} = 1 for @{$self->reference_sequence_build->chromosome_array_ref};
    while(my $line = $genotype_fh->getline) {
        chomp $line;

        my ($chr, $pos, $genotype) = split /\s+/, $line;

        #intersecting position
        #check genotypes
        my @alleles = split //, uc($genotype);
        unless(exists $chromosomes{$chr}) {
            $self->warning_message('Chromosome "%s" is not in the reference "%s". Skipping position.', $chr, $self->reference_sequence_build->name);
            next;
        }
        my $ref = $self->reference_sequence_build->sequence($chr, $pos, $pos);

        if($genotype ne '--' && $genotype =~ /[ACTGN]/ ) {

            #print genotypes with call
            my $type1 = ($alleles[0] eq $ref) ? 'ref' : 'SNP';
            my $type2 = ($alleles[1] eq $ref) ? 'ref' : 'SNP';
            print $output_fh "$chr\t$pos\t$pos\t$alleles[0]\t$alleles[1]\t$type1\t$type2\t$type1\t$type2\n";
        }

    }

    return 1;
}

1;

