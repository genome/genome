package Genome::Model::Tools::Predictor::Ber::Setup;

use strict;
use warnings;
use Carp;

use Genome;

class Genome::Model::Tools::Predictor::Ber::Setup {
    is  => ['Command'],
    has => [
        output_directory => {
            is => 'DirectoryPath',
            is_input => 1,
            predictor_specific => 1,
            doc => 'Directory in which raw and parsed output from this predictor should go',
        },
        input_fasta_file => {
            is => 'FilePath',
            is_input => 1,
            doc => 'File containing assembly sequence (typically fasta) to be used as input to predictor',
        },
    ],
};

sub execute { 

    my $self = shift;

    return $self->setup;
}

sub _get_locus_from_fasta_header {
    my $self = shift;
    
    #get headers
    my $output_fh = Genome::Sys->open_file_for_reading($self->input_fasta_file) 
        or die "Could not get file handle for " . $self->input_fasta_file;
    
    my $fasta_header = $output_fh->getline;
    chomp $fasta_header;
    
    #check header format 
    my ($locus_id) = $fasta_header =~ /^\>\w*?\-?(\w+)_Contig/;

    #get locus id and ensure only one exists
    return $locus_id;
}


sub setup {
    my $self = shift;
    my $locusid = $self->_get_locus_from_fasta_header;
    my $config_dir =$self->output_directory.'/BER/';
    make_path($config_dir) or croak "failed to make path: $config_dir\n";
    
    # link in source fasta for easy reference
    my $fasta_link = $config_dir.$self->id.'_input.fasta';
    Genome::Sys->create_symlink($self->input_fasta_file, $fasta_link) 
            or croak "failed to create link $fasta_link";

    my $fasta_dir = $config_dir.'fasta/';
    make_path($fasta_dir) or croak "failed to make path: $fasta_dir";
    my $hmm_dir   = $config_dir.'hmm/';
    make_path($hmm_dir) or croak "failed to make path: $hmm_dir";
    my $ber_dir   = $config_dir.'ber/';
    make_path($ber_dir) or croak "failed to make path: $ber_dir";

    #to satisfy the Ber software it expects a directory in its data path
    my $ber_genomes_dir = $self->tool_path_for_version.'/data/genomes/'.$locusid;
    Genome->symlink_directory($config_dir, $ber_genomes_dir) or 
        croak "failed to symlink dir $config_dir to $ber_genomes_dir";
    my $asm_file = $locusid.'_asm_feature';
    my $asm_fh = Genome::Sys->open_file_for_writing($config_dir.$asm_file) 
        or croak "cann't open the file:$!.";
    my $asmbl_file = $locusid.'_asmbl_data';
    my $asmbl_fh = Genome::Sys->open_file_for_writing($config_dir.$asmbl_file) 
        or croak "cann't open the file:$!.";
    my $indent_file = $locusid.'_indent_fh';
    my $indent_fh = Genome::Sys->open_file_for_writing($config_dir.$indent_file) 
        or croak "cann't open the file:$!.";
    my $stan_file = $locusid.'_stan';
    my $stan_fh = Genome::Sys->open_file_for_writing($config_dir.$stan_file) 
        or croak "cann't open the file:$!.";

    $asm_fh->print(join("\t", qw(asmbl_id end3 end5 feat_name feat_type))."\r\n");
    $asmbl_fh->print(join("\t", qw(id name type))."\r\n");
    $indent_fh->print(join("\t", qw(complete feat_name locus))."\r\n");
    $stan_fh->print(join("\t", qw(asmbl_data_id asmbl_id iscurrent))."\r\n");

    my $count=0;
    my $asmblid=0;
    my %fasta_files;
    my $last_header;
    my $last_gene;
    my $fh = Genome::Sys->open_file_for_reading($self->input_fasta_file) or croak "failed to open: ".$self->input_fasta_file;
    
    while (my $line = $fh->getline)  {
        if($line =~ $locusid) {
            my ($contig, $start, $stop) = $line =~ /\_(Contig.*)\s(\d+)\s(\d+)$/;
            my $locus_count = sprintf("%05d", $count);
            $asm_fh->print(join("\t", ($count, $start, $stop, $contig, 'ORF'))."\r\n");
            $asmbl_fh->print(join("\t", ($count, 'Contig','contig'))."\r\n");
            $indent_fh->print(join("\t", (' ', $contig, $locusid.$locus_count))."\r\n");
            $stan_fh->print(join("\t", ($count, $count, 1))."\r\n");
            $count++;
            $last_header = $line;
            $last_gene   = $contig;
        }
        else {
            my $gene_file = $fasta_dir.$last_gene.'.fasta';
            my $gene_fh = Genome::Sys->open_file_for_writing($gene_file)
                or croak "cann't open the file:$!.";
            $gene_fh->print($last_header);
            $gene_fh->print($line);
            $gene_fh->close;
            $fasta_files{$last_gene} = $gene_file;
        }
    }

    $asm_fh->close;
    $asmbl_fh->close;
    $indent_fh->close;
    $stan_fh->close;
    
    #link these files into the BER data directories also, BLAH!
    my $ber_db_csv_dir = $self->tool_path_for_version.'/data/db/CSV/';
    Genome->create_symlink($config_dir.$asm_file, $ber_db_csv_dir.$asm_file) or 
        croak "failed to symlink dir $config_dir to $ber_genomes_dir";
    Genome->create_symlink($config_dir.$asmbl_file, $ber_db_csv_dir.$asmbl_file) or 
        croak "failed to symlink dir $config_dir to $ber_genomes_dir";
    Genome->create_symlink($config_dir.$indent_file, $ber_db_csv_dir.$indent_file) or 
        croak "failed to symlink dir $config_dir to $ber_genomes_dir";
    Genome->create_symlink($config_dir.$stan_file, $ber_db_csv_dir.$stan_file) or 
        croak "failed to symlink dir $config_dir to $ber_genomes_dir";
    

    return 1;
}

1;
