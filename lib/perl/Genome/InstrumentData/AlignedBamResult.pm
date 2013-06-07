package Genome::InstrumentData::AlignedBamResult;

use Genome;

use warnings;
use strict;

class Genome::InstrumentData::AlignedBamResult {
    is => 'Genome::SoftwareResult::Stageable',
    is_abstract => 1,
    attributes_have => [
        is_output => { is => 'Boolean', is_optional => 1, },
    ],
    has_constant => [
        bam_file => { 
            is_output => 1,
            calculate => q| return $self->output_dir.'/'.$self->id.'.bam'; |, 
        },
    ],
};

sub run_flagstat_on_output_bam_file {
    my $self = shift;
    $self->status_message('Run flagstat on output bam file...');

    my $bam_file = $self->bam_file;
    if ( not $bam_file or not -s $bam_file ) {
        $self->error_message('Bam file not set or does not exist!');
        return;
    }

    my $flagstat_file = $bam_file.'.flagstat';
    $self->status_message("Flagstat file: $flagstat_file");
    my $cmd = "samtools flagstat $bam_file > $flagstat_file";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_file ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
    $self->status_message('Flagstat output:');
    $self->status_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );
    if ( not $flagstat->{total_reads} > 0 ) {
        $self->error_message('Flagstat determined that there are no reads in bam! '.$bam_file);
        return;
    }

    $self->status_message('Run flagstat on output bam file...done');
    return $flagstat;
}

1;

