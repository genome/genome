package Genome::Model::Tools::SnpArray::LohCalls;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#
#       AUTHOR:         Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#       CREATED:        04/01/2009 by D.K.
#       MODIFIED:       04/01/2009 by D.K.
#
#       NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();

class Genome::Model::Tools::SnpArray::LohCalls {
    is => 'Command',

    has => [ # specify the command's single-value properties (parameters) <---
             tumor_genotype_file => {
                 is => 'Text',
                 doc => "Three-column file of genotype calls chrom, pos, genotype",
                 is_optional => 0,
                 is_input => 1 },

             normal_genotype_file => {
                 is => 'Text',
                 doc => "Three-column file of genotype calls chrom, pos, genotype",
                 is_optional => 0,
                 is_input => 1 },

             output_file => {
                 is => 'Text',
                 doc => "Output file for QC result",
                 is_optional => 1,
                 is_input => 1},

             include_all_mismatches => {
                 is => 'Text',
                 doc => "If set to 1, include any site where genotypes differ",
                 is_optional => 1,
                 is_input => 1},

             cbs_output => {
                 is => 'Text',
                 doc => "If set to 1, output all sites in format appropriate for feeding to CBS",
                 is_optional => 1,
                 is_input => 1},

        ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Makes LOH calls from SNP array genotypes from tumor-normal pairs"
}

sub help_synopsis {
    return <<EOS
        This command compares SAMtools variant calls to array genotypes
EXAMPLE:        gmt snp-array loh-calls --tumor-genotype-file tumor.affy --normal-genotype-file normal.affy --output-file tumor-normal.loh
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
    my $tumor_genotype_file = $self->tumor_genotype_file;
    my $normal_genotype_file = $self->normal_genotype_file;
    my $output_file = $self->output_file;

    open(OUTFILE, ">$output_file") or die "Can't open outfile $output_file: $!\n";

    print "Loading tumor genotypes...\n";
    my %tumor_genotypes = load_genotypes($tumor_genotype_file);

    print "Loading normal genotypes...\n";
    my %normal_genotypes = load_genotypes($normal_genotype_file);


    ## Go through all normal genotypes ##

    foreach my $snp_key (keys %normal_genotypes)
    {
        $stats{'num_normal_genotypes'}++;

        if($tumor_genotypes{$snp_key})
        {
            $stats{'num_compared'}++;

            my $normal_gt = $normal_genotypes{$snp_key};
            my $tumor_gt = $tumor_genotypes{$snp_key};
            
            if($normal_gt eq $tumor_gt)
            {
                $stats{'num_matched'}++;
                if($self->cbs_output){
                    print OUTFILE join("\t", ($snp_key, "0")) . "\n";
                }
            }
            elsif(is_heterozygous($normal_gt) && is_homozygous($tumor_gt))
            {
                $stats{'num_LOH'}++;
                if($self->cbs_output){
                    print OUTFILE join("\t", ($snp_key, "1")) . "\n";
                } else {
                    print OUTFILE join("\t", $snp_key, $normal_gt, $tumor_gt) . "\n";
                }    
                
            } elsif(is_homozygous($normal_gt) && is_heterozygous($tumor_gt)) {
                $stats{'num_GOH'}++;
                if($self->include_all_mismatches)
                {
                    if($self->cbs_output){
                        print OUTFILE join("\t", ($snp_key, "NA")) . "\n";
                    } else {
                        print OUTFILE join("\t", $snp_key, $normal_gt, $tumor_gt) . "\n";
                    }
                }
            }
            else
            {
                $stats{'num_other_mismatch'}++;
                if($self->include_all_mismatches)
                {
                    print OUTFILE join("\t", $snp_key, $normal_gt, $tumor_gt) . "\n";
                }
            }
        }
    }

    print $stats{'num_normal_genotypes'} . " normal genotypes\n";
    print $stats{'num_compared'} . " had genotype call in tumor\n";
    print $stats{'num_matched'} . " matched the normal genotype\n";
    print $stats{'num_LOH'} . " were LOH\n";
    print $stats{'num_GOH'} . " were GOH\n";
    print $stats{'num_other_mismatch'} . " were another kind of mismatch\n";

    close(OUTFILE);
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_genotypes
{                               # replace with real execution logic.
    my $genotype_file = shift(@_);
    my %genotypes = ();

    my $input = new FileHandle ($genotype_file);
    my $lineCounter = 0;
    my $gtCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        (my $chrom, my $position, my $genotype) = split(/\t/, $line);

        if($chrom && $position)
        {
            my $key = "$chrom\t$position";

            if($genotype && $genotype ne "--" && $genotype ne "NN")
            {
                $genotypes{$key} = $genotype;
                $gtCounter++;
            }
        }

    }
    close($input);

    print "$gtCounter genotypes loaded\n";

#       print "$gtCounter genotypes loaded\n";

    return(%genotypes);                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub is_heterozygous
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);
    return(1) if($a1 ne $a2);
    return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub is_homozygous
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);
    return(1) if($a1 eq $a2);
    return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub flip_genotype
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);

    if($a1 eq "A")
    {
        $a1 = "T";
    }
    elsif($a1 eq "C")
    {
        $a1 = "G";
    }
    elsif($a1 eq "G")
    {
        $a1 = "C";
    }
    elsif($a1 eq "T")
    {
        $a1 = "A";
    }

    if($a2 eq "A")
    {
        $a2 = "T";
    }
    elsif($a2 eq "C")
    {
        $a2 = "G";
    }
    elsif($a2 eq "G")
    {
        $a2 = "C";
    }
    elsif($a2 eq "T")
    {
        $a2 = "A";
    }

    $gt = $a1 . $a2;
    $gt = sort_genotype($gt);
    return($gt);
}

################################################################################################
# Load Genotypes
#
################################################################################################

sub sort_genotype
{
    my $gt = shift(@_);
    (my $a1, my $a2) = split(//, $gt);

    my @unsorted = ($a1, $a2);
    my @sorted = sort @unsorted;
    $a1 = $sorted[0];
    $a2 = $sorted[1];
    return($a1 . $a2);
}



sub code_to_genotype
{
    my $code = shift(@_);

    return("AA") if($code eq "A");
    return("CC") if($code eq "C");
    return("GG") if($code eq "G");
    return("TT") if($code eq "T");

    return("AC") if($code eq "M");
    return("AG") if($code eq "R");
    return("AT") if($code eq "W");
    return("CG") if($code eq "S");
    return("CT") if($code eq "Y");
    return("GT") if($code eq "K");

#       warn "Unrecognized ambiguity code $code!\n";

    return("NN");
}



sub commify
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}


1;

