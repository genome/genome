package Genome::Model::Tools::Bed::Convert::Snv::RtgSnpToBed;

use strict;
use warnings;

use Genome;
use POSIX;
use Genome::Info::IUB;

class Genome::Model::Tools::Bed::Convert::Snv::RtgSnpToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert samtools SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv rtg-snp-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in samtools format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;
    my $input_fh = $self->_input_fh;

    #dump header lines
    for (1..5){
        $input_fh->getline;
    }

    while(my $line = <$input_fh>) {
        my ($chromosome, $position, undef, $reference, $consensus, $quality, $depth, @extra) = split("\t", $line);
        next if $consensus =~ m/^x$/;
        my $cons;
        my @cons = split ":", $consensus;
        $cons = (scalar(@cons)>1) ? Genome::Info::IUB->iub_for_alleles(@cons) : $consensus ;

        if((length $cons)> 1){
            next;
        }
        if((length $reference) > 1){
            next;
        }
        $quality = floor($quality);
        if(defined($cons)){
            $self->write_bed_line($chromosome, $position - 1, $position, $reference, $cons, $quality, $depth);
        } else {
            $self->error_message("Could not get an IUB code for ". $line);
        }
    }
    
    return 1;
}

sub convert_bed_to_detector {
    my $self = shift;
    my $detector_file = $self->detector_style_input;    
    my $bed_file = $self->source;
    my $output = $self->output;

    my $ofh = Genome::Sys->open_file_for_writing($output);
    my $detector_fh = Genome::Sys->open_file_for_reading($detector_file);
    my $bed_fh = Genome::Sys->open_file_for_reading($bed_file);
    OUTER: while(my $line = $bed_fh->getline){
        chomp $line;
        my ($chr,$start,$stop,$refvar,@data) = split "\t", $line;
        my ($ref,$var) = split "/", $refvar;
        my $scan=undef;
        while(my $dline = $detector_fh->getline){
            chomp $dline;
            my ($dchr,$dpos,$dref,$dvar) = split "\t", $dline;
            if(($chr eq $dchr)&&($stop == $dpos)&&($ref eq $dref)&&($var eq $dvar)){
                print $ofh $dline."\n";
                next OUTER;
            }
        }
    }
    $bed_fh->close;
    $ofh->close;
    $detector_fh->close;
    return 1;
}

1;
