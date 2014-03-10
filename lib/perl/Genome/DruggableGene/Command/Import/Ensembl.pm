package Genome::DruggableGene::Command::Import::Ensembl;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::Import::Ensembl {
    is => 'Genome::DruggableGene::Command::Import::Base',
    has => [
        genes_outfile => {
            is => 'Path',
            is_input => 1,
            example_values => ['/gscmnt/sata132/techd/mgriffit/DruggableGenes/TSV/Ensembl_WashU_TARGETS.tsv'],
            doc => 'PATH.  Path to .tsv file for genes (targets)',
        },
        input_gtf_url => {
            is => 'Text',
            is_input => 1,
            doc => 'URL PATH.  URL to an Ensembl transcript GTF file on the Ensembl FTP site (http://useast.ensembl.org/info/data/ftp/index.html)',
        },
        citation_base_url => {
            default => 'http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=',
        },
        citation_site_url => {
            default => 'http://ensembl.org/index.html',
        },
        citation_text => {
            default => "Ensembl 2011. Flicek P, Amode MR, ..., Vogel J, Searle SM. Nucleic Acids Res. 2011 Jan;39(Database issue)800-6. Epub 2010 Nov 2. PMID: 21045057.",
        },
    ],
    doc => 'Import genes from Ensembl transcript GTF file',
};


sub _doc_copyright_years {
    (2011);
}

sub _doc_license {
    my $self = shift;
    my (@y) = $self->_doc_copyright_years;  
    return <<EOS
Copyright (C) $y[0] Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the 
associated LICENSE file in this distribution.
EOS
}

sub _doc_authors {
    return <<EOS
 Malachi Griffith, Ph.D.
EOS
}

#Citations and related things
=cut
sub _doc_credits {
    return ('','None at this time.');
}
=cut

sub _doc_see_also {
    return <<EOS
B<genome>(1)
EOS
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
}

sub help_synopsis {
    return <<HELP
genome druggable-gene import ensembl --version=68_37 --input-gtf-url=ftp://ftp.ensembl.org/pub/release-68/gtf/homo_sapiens/Homo_sapiens.GRCh37.68.gtf.gz
HELP
}

sub help_detail {
#TODO: make this accurate
    my $summary = <<HELP
This importer downloads genes from a current Ensembl Human GTF file.
HELP
}

sub execute {
    my $self = shift;
    $self->input_to_tsv(); #Get the data from source files and write to a temp file
    $self->status_message(".tsv file created for use in 'rake dgidb:import:ensembl <.tsv_file_path> <source_db_version>'");
    return 1;
}

sub input_to_tsv {
    my $self = shift;
    my $genes_outfile_path = $self->genes_outfile;
    my $input_gtf_url = $self->input_gtf_url;

    #TODO: Take in the input ensembl file, make a tsv file at $genes_outfile_path  
    #The user will specify a URL path to a GTF file on the Ensembl FTP site and a version name
    #e.g. ftp://ftp.ensembl.org/pub/release-64/gtf/homo_sapiens/Homo_sapiens.GRCh37.64.gtf.gz
    my ($fh, $path) = Genome::Sys->create_temp_file('temp_gtf_file');
    my $wget_cmd = "wget $input_gtf_url -O $path";
    my $retval = Genome::Sys->shellcmd(cmd=>$wget_cmd);

    unless ($retval == 1){
      self->error_message('Failed to wget the specified URL');
      return;
    }

    my %ensembl_map;
    open (ENSG, "zcat $path |") || die "\n\nCould not open zcat pipe to file: $path\n\n";
    while(<ENSG>){
      chomp($_);
      my @line = split("\t", $_);
      # gene_id "ENSG00000237375"; transcript_id "ENST00000327822"; exon_number "1"; gene_name_report "BX072566.1"; gene_biotype "protein_coding"; transcript_name "BX072566.1-201";
      $DB::single=1;
      my @gene_info = split("; ", $line[8]);
      my ($gene_id_string) = grep ($_ =~ /gene_id/i, @gene_info);
      $gene_id_string =~ s/\"//g; #Kill the quotes
      $gene_id_string =~ s/ gene_id //i; #Kill the gene_id part leaving the actual ENSG id
      my $gene_id = uc($gene_id_string);

      my ($gene_name_report_string) = grep ($_ =~ /gene_name/i, @gene_info);
      $gene_name_report_string =~ s/\"//g; #Kill the quotes
      $gene_name_report_string =~ s/gene_name //i; #Kill the gene_id part leaving the actual ENSG id
      my $gene_name_report = uc($gene_name_report_string);

      my ($gene_biotype_string) = grep ($_ =~ /gene_biotype/i, @gene_info);
      my $gene_biotype = '';
      if (defined($gene_biotype_string)){
        $gene_biotype_string =~ s/\"//g; #Kill the quotes
        $gene_biotype_string =~ s/gene_biotype //i; #Kill the gene_id part leaving the actual ENSG id
        $gene_biotype = uc($gene_biotype_string);
      }else{
        $gene_biotype = "N/A";
      }
      $ensembl_map{$gene_id}{name} = $gene_name_report;
      $ensembl_map{$gene_id}{biotype} = $gene_biotype;
    }
    close(ENSG);
    $fh->close;

    open(TARGETS, ">$genes_outfile_path") || die "\n\nCould not open outfile: $genes_outfile_path\n\n";
    my $targets_header = join("\t", 'ensembl_id', 'ensembl_gene_symbol', 'ensembl_gene_biotype');
    print TARGETS "$targets_header\n";
    for my $ensembl_id (keys %ensembl_map){
      print TARGETS join("\t", $ensembl_id, $ensembl_map{$ensembl_id}{name}, $ensembl_map{$ensembl_id}{biotype}), "\n";
    }
    close(TARGETS);

    return 1;
}

1;
