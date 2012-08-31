package Genome::Model::Tools::Annotate::CosmicUpdate;

########################################################################################################################
# CosmicUpdate.pm - A module for updating the local COSMIC database to reflect the current COSMIC database.
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	3/09/2010 by W.S.
#	MODIFIED:	3/09/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

   use warnings;
   use strict;
   use Net::FTP;
   use Cwd;
   use Getopt::Long;

class Genome::Model::Tools::Annotate::CosmicUpdate {

	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		cosmic_folder	=> { is => 'Text', doc => "Path to the current local COSMIC files", is_optional => 1 , default => '/gscuser/wschierd/git/genome/genome-db-cosmic/' },
		cosmic_url	=> { is => 'Text', doc => "URL to the online COSMIC repository", is_optional => 1 , default => 'ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/genes/' },
		output_file	=> { is => 'Text', doc => "Output file name for flatfile of amino acid changes" , is_optional => 1 , default => 'Cosmic_Database.tsv' },
	],
};

sub help_brief {                            # keep this to just a few words <---
    "Update Local COSMIC database"
}

sub help_synopsis {
    return <<EOS
A module for updating the downloaded COSMIC database to reflect the current COSMIC database.
EXAMPLE:	gmt annotate cosmic-update ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
A module for updating the downloaded COSMIC database to reflect the current COSMIC database.
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.

my $self = shift;
my $URL = $self->cosmic_url;
my $cosmicdir = $self->cosmic_folder;
my $cosmicdb = $self->output_file;

chdir ($cosmicdir);
my $dir = getcwd;
print "Working Directory: $dir";

print "Retrieving File(s) from COSMIC\n";

#directory on COSMIC server with all GENE folders
my $wget = 'wget';
system ( $wget, "\-r" , "\-A.csv" , "\-nd" , "\-N" , "\-l2" , $URL);

print "Finished Downloading COSMIC File(s), Moving on to Writing COSMIC Flatfile\n";

opendir(IMD, $cosmicdir) || die("Cannot open directory"); 

my @cosmicdirfiles= readdir(IMD);
closedir(IMD);
my @cosmiccsvfiles;

foreach my $file (@cosmicdirfiles) {
  unless ( ($file eq ".") || ($file eq "..") ) {
    if ( $file =~ m/.csv/ ) {
       push (@cosmiccsvfiles,$file);
    }
  }
}


unless (open(COSMIC_DB,">$cosmicdb")) {
    die "Could not open output file '$cosmicdb' for writing";
}
#Sample Name	COSMIC Sample ID	Amino Acid	Nucleotide	Primary Tissue	Tissue subtype 1	Tissue subtype 2	Histology	Histology subtype 1	Histology subtype 2	Pubmed ID	studies	Mutation ID	Somatic Status	Sample Source	Zygosity	Chromosome NCBI36	Genome Start NCBI36	Genome Stop NCBI36	Chromosome GRCh37	Genome Start GRCh37	Genome Stop GRCh37
#Sample Name       COSMIC Sample ID        Amino Acid      Nucleotide      Primary Tissue  Tissue subtype 1        Tissue subtype 2        Histology       Histology subtype 1     Histology subtype 2     Pubmed ID       studies Mutation ID     Somatic Status  Sample Source Zygosity        Chromosome      Genome Start    Genome Stop

print COSMIC_DB "Gene\tChromosome\tGenome Start\tGenome Stop\tAmino Acid\tNucleotide\tSomatic Status\tPrimary_Tissue\tTissue_subtype_1\tTissue_subtype_2\tHistology\tHistology_subtype_1\tHistology_subtype_2\tChromosome Build37\tGenome Start Build37\tGenome Stop Build37\n";

my $i = 0;
foreach my $genefile (@cosmiccsvfiles) {
	print ".";
	my $gene = $genefile;
	$gene =~ s/.csv//g;
	chomp ($gene);
	unless (open(GENE_FILE,"<$genefile")) {
	    die "Could not open output file '$genefile' for writing";
	}
	my $firstline = 0;
	my $chr_col;
	my $start_col;
	my $stop_col;
	my $chr_col_37;
	my $start_col_37;
	my $stop_col_37;
	my $amino_col;
	my $nucleo_col;
	my $somatic_col;
	my $primary_tissue_col;
	my $tissue_sub_1_col;
	my $tissue_sub_2_col;
	my $histology_col;
	my $histology_sub_1_col;
	my $histology_sub_2_col;
	my $chr;
	my $start;
	my $stop;
	my $chr_37;
	my $start_37;
	my $stop_37;
	my $amino;
	my $nucleo;
	my $somatic;
	my $primary_tissue;
	my $tissue_sub_1;
	my $tissue_sub_2;
	my $histology;
	my $histology_sub_1;
	my $histology_sub_2;
	while (my $line = <GENE_FILE>) {
		if ($line =~ m/(Amino Acid)/ ) {
			$firstline = 1;
			my @parser = split(/\t/, $line);
			my $parsecount = 0;
			my %parsehash;
			foreach my $item (@parser) {
				$parsehash{$item} = $parsecount;
				$parsecount++;
			}
			if ($line =~ m/(GRCh37)/) {
				$chr_col = $parsehash{'Chromosome NCBI36'};
				$start_col = $parsehash{'Genome Start NCBI36'};
				$stop_col = $parsehash{'Genome Stop NCBI36'};
				$chr_col_37 = $parsehash{'Chromosome GRCh37'};
				$start_col_37 = $parsehash{'Genome Start GRCh37'};
				$stop_col_37 = $parsehash{'Genome Stop GRCh37'};
			}
			else {
				$chr_col = $parsehash{'Chromosome'};
				$start_col = $parsehash{'Genome Start'};
				$stop_col = $parsehash{'Genome Stop'};
				$chr_col_37 = 'NotListed';
				$start_col_37 = 'NotListed';
				$stop_col_37 = 'NotListed';
			}
			$amino_col = $parsehash{'Amino Acid'};
			$nucleo_col = $parsehash{'Nucleotide'};
			$somatic_col = $parsehash{'Somatic Status'};
			$primary_tissue_col = $parsehash{'Primary Tissue'};
			$tissue_sub_1_col = $parsehash{'Tissue subtype 1'};
			$tissue_sub_2_col = $parsehash{'Tissue subtype 2'};
			$histology_col = $parsehash{'Histology'};
			$histology_sub_1_col = $parsehash{'Histology subtype 1'};
			$histology_sub_2_col = $parsehash{'Histology subtype 2'};
	
			unless (defined($chr_col) && defined($start_col) && defined($stop_col) && defined($amino_col) && defined($nucleo_col) && defined($somatic_col)) {
			    die "Line: $line\nAbove line could not be parsed for gene: $gene";
			}
			next;
		}
		unless ($firstline == 1) {
			next;
		}
		my @parser = split(/\t/, $line);
		chomp($parser[$chr_col],$parser[$start_col],$parser[$stop_col],$parser[$amino_col],$parser[$nucleo_col],$parser[$somatic_col]);
		unless ($parser[$chr_col] ne ' ' || $parser[$start_col] ne ' ' || $parser[$stop_col] ne ' ' || $parser[$amino_col] ne ' ' || $parser[$nucleo_col] ne ' ' || $parser[$somatic_col] ne ' ') {
		    next;
		}
		$chr = $parser[$chr_col];
		$start = $parser[$start_col];
		$stop = $parser[$stop_col];
		$amino = $parser[$amino_col];
		$nucleo = $parser[$nucleo_col];
		$somatic = $parser[$somatic_col];
		$primary_tissue = $parser[$primary_tissue_col];
		$tissue_sub_1 = $parser[$tissue_sub_1_col];
		$tissue_sub_2 = $parser[$tissue_sub_2_col];
		$histology = $parser[$histology_col];
		$histology_sub_1 = $parser[$histology_sub_1_col];
		$histology_sub_2 = $parser[$histology_sub_2_col];
		if ($chr_col_37 eq 'NotListed') {
			$chr_37 = $chr_col_37;
			$start_37 = $start_col_37;
			$stop_37 = $stop_col_37;
		}
		else {
			$chr_37 = $parser[$chr_col_37];
			$start_37 = $parser[$start_col_37];
			$stop_37 = $parser[$stop_col_37];
		}
		chomp($chr,$start,$stop,$amino,$nucleo,$somatic,$primary_tissue,$tissue_sub_1,$tissue_sub_2,$histology,$histology_sub_1,$histology_sub_2,$chr_37,$start_37,$stop_37);
		print COSMIC_DB "$gene\t$chr\t$start\t$stop\t$amino\t$nucleo\t$somatic\t$primary_tissue\t$tissue_sub_1\t$tissue_sub_2\t$histology\t$histology_sub_1\t$histology_sub_2\t$chr_37\t$start_37\t$stop_37\n";
	}
}


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}
