package Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;

class Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed {
    is => 'Genome::Command::Base',
    doc => 'Convert the genotype microarray file to the .bed format used by the GoldSnp view',
    has_input => [
        input_file => {
            is => 'File',
            is_optional => 1,
            doc => 'The input file Genotype Microarray file to convert.',
            shell_args_position => 1,
        },
        output_file => {
            is => 'File',
            is_optional => 1,
            doc => 'The output BED file to write.',
            shell_args_position => 2,
        },
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_optional => 1,
            id_by => 'reference_id',
            doc => 'The reference sequence build id to use',
            shell_args_position => 3,
        },
        reference_id => {
            is => 'Integer',
            is_optional => 1,
        },
        build_id => {
            is => 'Number',
            is_optional => 1,
            doc => 'ID of genotype microarray build',
        },
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            is_optional => 1,
            id_by => 'build_id',
            doc => 'GenotypeMicroarray build',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->{bases_files} = {};

    my ($input_file, $output_file, $reference) = ($self->input_file, $self->output_file, $self->reference);
    unless (defined $input_file and defined $output_file and defined $reference) {
        unless ($self->build) {
            Carp::confess 'Must give either a build id or an input file, output file, and reference id!';
        }

        $input_file = $self->build->formatted_genotype_file_path;
        $output_file = $self->build->snvs_bed;
        $reference = $self->build->model->reference_sequence_build;
    }

    # Unfortunately, some of these genotype files seem to be unsorted.
    # It would be better sort first, so the reference sequence can be
    # accessed sequentially. However, joinx sort can't deal with -- in the
    # stop field of the input, which seems to happen sometimes, so we let 
    # _convert_file filter those out then sort after.

    my $tmpfile = Genome::Sys->create_temp_file_path;
    my $infile = $input_file;
    my $ifh = new IO::File("<$infile") || die "Failed to open input file $infile";
    my $ofh = new IO::File(">$tmpfile") || die "Failed to open output file $tmpfile";
    my $rv = $self->_convert_file($reference, $ifh, $ofh);
    $ifh->close();
    $ofh->close();

    my $sort_cmd = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [$tmpfile],
        output_file => $output_file,
    );
    if (!$sort_cmd->execute) {
        $self->error_message("Failed to sort input genotype file " . $input_file);
        return;
    }

    return $rv;
}

sub _convert_file {
    my ($self, $reference, $ifh, $ofh) = @_;

    while (my $line = <$ifh>) {
        chomp $line;
        my ($chr, $start, $end, $allele1, $allele2, $a1t1, $a2t1, $a1t2, $a2t2) = split("\t", $line);

        if ($end =~ /^-/) {
            $self->warning_message("Invalid end position '$end', skipping");
            next;
        }

        my $ref = $reference->sequence($chr, $end, $end); 
        if (!$ref) {
            $self->warning_message("No reference for position $chr $end $end, skipping $line\n");
            next;
        }
        $ref = uc($ref);

        if ($allele1 eq $allele2) {
            if ($a1t1 ne $a2t1 or $a1t2 ne $a2t2) {
                print STDERR "Inconsistent types within a platform $line\n";
                next;
            }

            if ($a1t1 eq 'ref' and $allele1 ne $ref) {
                print STDERR "Reference mismatch at $chr $end: expected $ref, got $allele1, line: $line\n";
                next;
            }

            if ($a1t1 ne 'ref' and $allele1 eq $ref) {
                print STDERR "Gold SNP is listed as reference in B36 sequence: $line\n";
                next;
            }
        } else {
            if ($a1t1 eq $a2t1) {
                print STDERR "Het snp where both alleles are non-ref: $line\n";
                next;
            }

            if ($a1t1 ne 'SNP' and $a2t1 ne 'SNP') {
                print STDERR "Supposedly het snp not labeled as such: $line\n";
                next;
            }
        }

        my $iub = Genome::Info::IUB::iub_for_alleles($allele1, $allele2);

        my $new_call = "$ref/$iub";
        --$start;
        my $score = 0;
        my $depth = 0;
        $ofh->print(join("\t", $chr, $start, $end, $new_call, $score, $depth, $a1t1, $a2t1, $a1t2, $a2t2) . "\n");

    }
    $ifh->close();
    $ofh->close();
    return 1;
}

1;
