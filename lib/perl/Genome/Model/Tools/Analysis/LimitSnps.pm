package Genome::Model::Tools::Analysis::LimitSnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LimiSnps.pm -    Limit SNPs by various factors, like ROI positions.
#               
#   AUTHOR:      Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#   CREATED:   03/20/2009 by D.K.
#   MODIFIED:   03/20/2009 by D.K.
#
#   NOTES:   
#         
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::LimitSnps {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
    variants_file   => { is => 'Text', doc => "File containing SNPs" },
    regions_file   => { is => 'Text', doc => "File containing ROI chromosome positions (chr1\t939339)"},
    reference_list => {is => 'String', doc => 'File containing the lengths of each reference sequence', is_optional => 1, example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai']},
    output_file   => { is => 'Text', doc => "Output file to receive limited SNPs", is_optional => 1 },
    not_file   => { is => 'Text', doc => "Output file to receive excluded SNPs", is_optional => 1 },
    verbose   => { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Limit a SNPs file to set regions"                 
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $variants_file = $self->variants_file;
    my $regions_file = $self->regions_file;
    my $reference_list = $self->reference_list;
    my $notfile = $self->not_file if(defined($self->not_file));
    my $outfile = $self->output_file if(defined($self->output_file));
    my $verbose = 1 if(defined($self->verbose));

    unless(-e $regions_file) {
        $self->error_message("Positions file does not exist!");
        return;
    }
    unless(-e $variants_file) {
        $self->error_message("Variants file does not exist!");
        return;
    }
    unless(-e $reference_list) {
        $self->error_message("Reference list does not exist!");
        return;
    }

    my %stats = ();
    $stats{'total_snps'} = $stats{'limit_snps'} = 0;
    my $input, my $lineCounter;

    my %genome;
    $input = new FileHandle ($reference_list);
    unless($input) {
        $self->error_message("Couldn't open $reference_list");
        return;
    }

    #set up bit vector representative of genome
    while(my $line = $input->getline) {
        chomp $line;
        #parse the reference name and length
        my ($chr, $length) = split /\s+/, $line;
        $chr = uc $chr;
        $genome{$chr} = Bit::Vector->new($length);
    }
    $input->close; #done creating the bit vector representation of the genome

    $self->debug_message("Parsing regions...");

    if($regions_file) {
        ## Parse the alignment blocks ##

        $input = new FileHandle ($regions_file);
        $lineCounter = 0;

        while (<$input>) {
            chomp;
            my $line = $_;
            my ($chrom, $start, $stop, $name) = split(/\t/, $line);
            $chrom = uc $chrom;
            unless(exists($genome{$chrom})) {
                $self->error_message("$chrom not found in reference list. Skipping region $name");
                next;
            }

            #assuming 1-based regions file
            unless($start >= 1 && $stop <= ($genome{$chrom}->Size() - 1)) {
                $self->error_message("Coordinates for region $name fall outside of the chromosome coordinates indicated in the reference list. Skipping region $name");
                next;
            }

            #otherwise add region to bitvector
            $genome{$chrom}->Interval_Fill($start-1, $stop-1);
        }
        $lineCounter = $input->input_line_number;
        $self->debug_message("$lineCounter regions loaded\n");
        $input->close;
    }

    print "Parsing SNPs...\n";
    ## Parse the combined SNPs ##

    ## Open the outfile ##
    if($outfile)
    {
        open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
    }

    if($notfile)
    {
        open(NOTFILE, ">$notfile") or die "Can't open outfile: $!\n";
    }

    $input = new FileHandle ($variants_file);
    $lineCounter = 0;
    while (<$input>)
    {
        chomp;
        my $line = $_;

        $stats{'total_snps'}++;
        (my $chrom, my $position) = split(/\t/, $line);
        $chrom = uc $chrom;

        if(!$regions_file || (exists($genome{$chrom}) && $genome{$chrom}->contains($position-1)))
        {
            $stats{'limit_snps'}++;
            print OUTFILE "$line\n" if($outfile);
        }
        else {
            print NOTFILE "$line\n" if($notfile);
        }
    }

    close($input);

    close(OUTFILE) if($outfile);
    close(NOTFILE) if($notfile);

    print "$stats{'total_snps'} total SNPs\n";

    ## Get percentage ##

    if($stats{'total_snps'})
    {
        $stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%';
        print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after filter\n";
    }




    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

