package Genome::Model::Tools::Dbsnp::Import::Flatfile;

use strict;
use warnings;

use Genome;
use LWP::Simple;
use List::AllUtils qw( :all );

class Genome::Model::Tools::Dbsnp::Import::Flatfile {
    is => 'Genome::Model::Tools::Dbsnp',
    has => 
        [
            flatfile => {
                is => 'Text',
                is_input => 1,
                doc => 'URL of the dbsnp flat file',
            },
            output_file => {
                is => 'Path',
                doc => 'File that dbsnps should be written/appended to',
            },
            reference_coordinates => {
                is => 'Text',
                is_input => 1,
                default => 'GRCh37\.p[0-9]+',
                doc => 'reference_coordinates whose coordinates will be used, regex syntax accepted for matching multiple patch levels',
            },
            use_contig => {
                is => 'Boolean',
                is_input => 1,
                default => 0,
                doc => 'True if the contig coords should be used.  Otherwise the chromosome coords will be used',
            },
        ],
        has_optional => [
            contig_name_translation_file => {
                is => 'Path',
                doc => 'File path that contains translations of contig names',
            },
            from_names_column => {
                is => 'Number',
                doc => '0-based column number containing names you want to translate from',
            },
            to_names_column => {
                is => 'Number',
                doc => '0-based column number containing names you want to translate to',
            },
        ],
        has_optional_transient => [
            _translate => {
                is => 'HASH',
            },
        ],
};

sub help_brief {
    'Create formatted tsv from DbSnp flat file, with one output file per contig'
}

sub help_synopsis {
    return <<EOS
gmt dbsnp import flatfile --flatfile ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/ASN1_flat/ds_flat_ch1.flat.gz --output-file /path/to/dir/ds.flat.bed
EOS
}

sub help_detail {
    return <<EOS
This command is used for importing flat file based DbSnp files.  It creates a .tsv that can be
merged with other such files and created as an ImportedVariationList build.
EOS
}

my %ds_type_conv = ( 'in-del' => 'INDEL',
                  'microsatellite' => 'MICROSATELLITE',
                  'mixed' => 'MIXED',
                  'multinucleotide-polymorphism' => 'MNP',
                  'named-locus' => 'NAMEDLOCUS',
                  'snp' => 'SNP',
                  'heterozygous' => 'HETEROZYGOUS',
                );
my %val_type_conv = ( 'by2Hit2Allele' => 'is_validated_by_allele',
                      'byCluster' => 'is_validated_by_cluster',
                      'byFrequency' => 'is_validated_by_frequency',
                      'byHapMap' => 'is_validated_by_hap_map',
                      'byOtherPop' => 'is_validated_by_other_pop',
                    );
        
my $csv_delimiter = "\t";
my @fd_order = qw(
    ds_chr
    ds_start
    ds_stop
    ds_allele
    ds_id
    ds_type
    submitter
    rs_id
    strain
    is_validated
    is_validated_by_allele
    is_validated_by_cluster
    is_validated_by_frequency
    is_validated_by_hap_map
    is_validated_by_other_pop
);


sub execute {
    my $self = shift;
    $self->debug_message("Processing file ".$self->flatfile);
    $self->debug_message("Using reference_coordinates $self->reference_coordinates");

    if ($self->contig_name_translation_file) {
        my %translate;
        my $translate_fh = Genome::Sys->open_file_for_reading($self->contig_name_translation_file);
        while (my $line = <$translate_fh>) {
            chomp $line;
            my @fields = split /\t/, $line;
            $translate{$fields[$self->from_names_column]} = $fields[$self->to_names_column];
        }
        $translate_fh->close;
        $self->_translate(\%translate);
    }
    my $flatfile_fh; 

    my $flatfile_path = Genome::Sys->create_temp_file_path;

    my $response = getstore($self->flatfile, $flatfile_path);

    die($self->error_message("Unable to download the flat files at: " . $self->flatfile)) unless $response ==     
    RC_OK;

    unless (-s $flatfile_path) {
        $self->warning_message("File ".$flatfile_path." is empty, continuing");
        return 1;
    }

    if (Genome::Sys->file_is_gzipped($flatfile_path)) {
        $flatfile_fh = Genome::Sys->open_gzip_file_for_reading($flatfile_path);
    }
    else {
        $flatfile_fh = Genome::Sys->open_file_for_reading($flatfile_path);
    }
    my $output_fh;
    if (-s $self->output_file) {
        $output_fh = new IO::File;
        my $open_string = ">> ".$self->output_file;
        $output_fh->open($open_string);
    }
    else {
        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    }

    my @block = ();
    while (<$flatfile_fh>) {
        chomp;
        next if ($. <= 3); # each file has a 3-line header
            my @split_line = split(/\s*\|\s*/, $_);
        if (@split_line == 0) { # blank line
            $self->process_block($output_fh, @block);
            @block = ();
        } else {
            push @block, \@split_line;
        }
    }
    $flatfile_fh->close;
    $output_fh->close;
    return 1;
}

sub process_block {
    # ensure the ss_pick=YES field is first
    my $self = shift;
    my $output_fh = shift;
    my $reference_coordinates = $self->reference_coordinates;
    my @ss = sort { $b->[-1] cmp $a->[-1] } (grep { $_->[0] =~ /^ss/ } @_);
    my @submitters = uniq map { $_->[1] } @ss;
    my ($snp) = grep { $_->[0] eq 'SNP' } @_;
    my ($val) = grep { $_->[0] eq 'VAL' } @_;
    my @ctgs = grep { $_->[0] eq 'CTG' && $_->[1] =~ /assembly=$reference_coordinates/ } @_;
    
    my $converted_type = $ds_type_conv{$_[0][3]};
    unless ($converted_type) {
        $self->warning_message("Could not convert type ".$_[0][3]);
    }
    my %record = ('ds_id'        => 0,
                  'rs_id'        => $_[0][0],
                  'ds_type'      => $converted_type,
                  'is_validated' => ($val->[1] eq 'validated=YES') || 0,
                  'is_validated_by_allele'    => 0,
                  'is_validated_by_cluster'   => 0,
                  'is_validated_by_frequency' => 0,
                  'is_validated_by_hap_map'   => 0,
                  'is_validated_by_other_pop' => 0,
    );


    my ($alleles) = ($snp->[1] =~ /alleles=\'(.*)\'/);
    my @ref_var = split('/', $alleles);
    
    ($record{'ds_allele'}) = ($snp->[1] =~ /alleles=\'(.*)\'/);

    if ($record{'is_validated'} and $val->[5]){
        my @vals = split /,/, $val->[5];
        foreach my $val_type (@vals){
            #TODO: this is janky, but it quiets the warnings for now
            next if $val_type eq "suspect";
            unless($val_type_conv{$val_type}) {
                $self->warning_message("$val_type not found in table");
            }
            $record{$val_type_conv{$val_type}} = 1;
        }
    }
    for my $ctg (@ctgs){
        my $contig_name;
        if ($self->use_contig) {
            $contig_name = $ctg->[4] or next;
            ($record{'ds_start'}) = ($ctg->[5] =~ /ctg-start=(\d+)/) or next;
            ($record{'ds_stop'})   = ($ctg->[6] =~ /ctg-end=(\d+)/) or next;
        }
        else {
            ($contig_name) = ($ctg->[2] =~ /chr=(.*)/) or next;
            my ($chr_pos) = ($ctg->[3] =~ /chr-pos=(\d+)/) or next;
            my ($ctg_start) = ($ctg->[5] =~ /ctg-start=(\d+)/) or next;
            my ($ctg_end)   = ($ctg->[6] =~ /ctg-end=(\d+)/) or next;

            $record{'ds_start'} = $chr_pos-1;
            $record{'ds_stop'}  = $chr_pos + ($ctg_end - $ctg_start);
        }
        if ($self->contig_name_translation_file) {
            if ($self->_translate->{$contig_name}) {
                $contig_name = $self->_translate->{$contig_name};
            }
            else {
                next;
            }
        }
        ($record{'strain'}) = ($ctg->[8] =~ /orient=([-\+])/) or next;
        $record{'ds_chr'} = $contig_name;
        my ($loctype)   = ($ctg->[7] =~ /loctype=(\d)/) or next;

        for my $sub (@submitters){
            $record{'submitter'} = $sub;
            my @vals = map { $record{$_} } @fd_order;
            $output_fh->print(join($csv_delimiter,@vals), "\n"); 
        }
    }
}

1;
