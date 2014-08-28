package Genome::Model::ClinSeq::Command::GenerateClonalityPlots::Readcounts;

use strict;
use warnings;

use Genome;

class Genome::Model::ClinSeq::Command::GenerateClonalityPlots::Readcounts {
    is => 'Command::V2',
    has_input => [
        sites_file => {
            is => 'Text',
            doc => 'List of chromosome sites for which to get BAM readcounts (format: "1	16949959	16949959	G	R")',
        },
        bam_files => {
            is => 'Text',
            doc => 'The BAM files to interrogate, with labels (e.g. "Tumor:/path/to/tumor.bam")',
            is_many => 1,
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference used to align the data in the BAMs',
        },
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'Where to save the results',
        },
        bam_readcount_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'The version of bam-readcounts to use',
        },
    ],
    doc => 'Get readcounts from a series of BAMs for the positions specified in the sites file',
};

use Cwd 'abs_path';
sub execute {
    my $self = shift;

    my $reference_fasta = $self->reference_build->full_consensus_path('fa');
    my @bams = $self->bam_files;
    my @types;
    my $out_fh = Genome::Sys->open_file_for_writing($self->output_file);

    #headers
    print $out_fh join("\t", '#Chr', 'Start', 'Stop');
    for my $bam (@bams) {
        my ($type,$bam) = split /:/,$bam;
        push @types,$type;
        print $out_fh map { "\t${type}_$_" } qw(Ref Var Ref_Count Var_Count Var_Freq);
    }
    print $out_fh "\n";

    my ($site_list_for_readcount_fh, $site_list_for_readcount) = Genome::Sys->create_temp_file();
    my $sites_fh = Genome::Sys->open_file_for_reading($self->sites_file);

    #begin to load output hash
    my %output;
    while (my $line = $sites_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop) = split /\t/,$line;
        $output{$chr}{$start} = {};

        $site_list_for_readcount_fh->say(join("\t", $chr, $start, $stop));
    }
    $sites_fh->close;

    #TODO Bring script into a module
    my $script_dir = Cwd::abs_path(File::Basename::dirname(__FILE__) . '/../../OriginalScripts/') . '/';
    unless (-d $script_dir) {
        die $self->error_message("failed to find script dir $script_dir!")
    }

    #cycle through bams to get readcounts
    for my $bam (@bams) {
        my $temp_readcount_file = Genome::Sys->create_temp_file_path();
        my ($type,$bam) = split /:/,$bam;

        my $rv = Genome::Model::Tools::Sam::Readcount->execute(
            bam_file => $bam,
            minimum_mapping_quality => 1,
            output_file => $temp_readcount_file,
            reference_fasta => $reference_fasta,
            region_list => $site_list_for_readcount,
            use_version => $self->bam_readcount_version,
        );

        my $sites = $self->sites_file;
        my $statscmd = "$script_dir/borrowed/ndees/miller-bam-readcount-to-stats.noheader.pl $sites $temp_readcount_file |";
        open(STATS,$statscmd) or die "Couldn't open stats command: $!";
        while (my $line = <STATS>) {
            chomp $line;
            my ($chr,$pos,$stats) = split /\t/,$line,3;
            $output{$chr}{$pos}{$type} = $stats;
        }
    }

    #print output
    for my $chr (sort keys %output) {
        for my $pos (sort keys %{$output{$chr}}) {
            print $out_fh join("\t", $chr, $pos, $pos);
            for my $type (@types) {
                if (exists $output{$chr}{$pos}{$type}) {
                    print $out_fh "\t".$output{$chr}{$pos}{$type};
                } else {
                    print $out_fh ("\tNA" x 5);
                }
            }
            print $out_fh "\n";
        }
    }

    return 1;
}

1;
