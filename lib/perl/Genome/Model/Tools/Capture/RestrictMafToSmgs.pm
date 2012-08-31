package Genome::Model::Tools::Capture::RestrictMafToSmgs;     # rename this when you give the module file a different name <--

use strict;
use warnings;
use FileHandle;
use Genome; # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::RestrictMafToSmgs {
  is => 'Command',

  has => [ # specify the command's single-value properties (parameters) <---
    maf_file  => { is => 'Text', doc => "Maf File", is_optional => 0 },
    smg_file  => { is => 'Text', doc => "Output from SMG Test", is_optional => 0 },
    output_file  => { is => 'Text', doc => "Output MAF restricted to SMGs", is_optional => 0 },
    output_bed_smgs  => { is => 'Text', doc => "Output bed file of MAF variants restricted to SMGs", is_optional => 1 },
  ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
  "Input a MAF and output a MAF restricted to SMGs"
}

sub help_synopsis {
  return <<EOS
Input a MAF and output a MAF restricted to SMGs
EXAMPLE:   gmt capture restrict-maf-to-smgs ...
EOS
}

sub help_detail {
  return <<HELP
Input a MAF and output a MAF restricted to SMGs
HELP
}

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;
    ## Get required parameters ##
    my $maf_file = $self->maf_file;
    my $smg_file = $self->smg_file;
    my $output_file = $self->output_file;
    my $fh = Genome::Sys->open_file_for_writing($output_file);
    my $fh2;
    if($self->output_bed_smgs) {
        $fh2 = Genome::Sys->open_file_for_writing($self->output_bed_smgs);
    }

    my %gene_hash;
    my $inFh_positions = IO::File->new( $smg_file ) || die "can't open $smg_file\n";
    while(my $line = $inFh_positions->getline ) {
        chomp($line);
        my ($gene_name, @rest) = split(/\t/, $line);
        if ($gene_name eq 'Gene') {
            next;
        }
        $gene_hash{$gene_name}++;
    }


    my $inFh_maf = IO::File->new( $maf_file ) || die "can't open $maf_file\n";
    while(my $line = $inFh_maf->getline ) {
        chomp($line);
        my ($gene_name, @rest) = split(/\t/, $line);
        if ($line =~ m/^#/ || $gene_name eq 'Hugo' || defined($gene_hash{$gene_name})) {
            print $fh "$line\n";
            next;
        }
        if($self->output_bed_smgs) {
            my $chr = $rest[3];
            my $start = $rest[4];
            my $stop = $rest[5];
            my $type = $rest[8];
            my $ref = $rest[9];
            my $var;
            if ($rest[9] eq $rest[10]) {
                $var = $rest[11];
            }
            else {
                $var = $rest[10];
            }
            if ($type eq "INS") {
                --$stop;
            }
            elsif ($type eq "DEL" || $type eq "SNV") {
                --$start;
            }
            print $fh2 "$chr\t$start\t$stop\t$ref\t$var\n";
        }
    }
    return 1; # No error

}


1;




