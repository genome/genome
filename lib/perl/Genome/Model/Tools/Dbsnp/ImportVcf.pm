package Genome::Model::Tools::Dbsnp::ImportVcf;

use warnings;
use strict;

use Genome;
use LWP::Simple;
use List::MoreUtils 'uniq';

class Genome::Model::Tools::Dbsnp::ImportVcf {
    is => 'Command::V2',
    has_input => [
        vcf_file_url => {
            is => 'Text',
            doc => 'Path to the full VCF file on the ftp site'
        },
        flat_file_pattern => {
            is => 'Text',
            default => 'ds_flat_ch[X].flat.gz',
            is_optional => 1,
            doc => 'String representing the pattern that the flat file filenames follow with [X] substituted in for the chromosome number'
        },
        output_file_path => {
            is => 'Path',
            is_output => 1,
            doc => 'Path to the final amended VCF file',
        },
        flat_url => {
            is_optional => 1,
            is => 'Text',
            doc => "URL of directory of flat files on ftp site",
        },
    ],
    has_transient_optional => [
        _submitter_map => {
            is => 'HASH',
            default => {},
        },
        _current_downloaded_chromosome => {
            is => 'Text',
            default => '0',
        }
    ],
};


sub help_brief {
    'Import a Dbsnp VCF file from the web - merging in submitter data from the flatfile format'
}

sub help_synopsis {
    return <<EOS
gmt dbsnp import-vcf --vcf_file_url ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz --output_file out
EOS
}

sub help_detail {
    return <<EOS
This command takes the path for a VCF file on the web, downloads it and its corresponding flatfiles, and merges submitter data in from the flatfiles.
EOS
}

sub execute {
    my $self = shift;
    #have we started the ##INFO block or are we ending it?
    my $start_info = 0;
    my $end_info = 0;

    my $vcf_download_location = Genome::Sys->create_temp_file_path();

    #check for successful status, else die
    my $response = getstore($self->vcf_file_url, $vcf_download_location);

    die($self->error_message("Unable to download the VCF file at: " . $self->vcf_file_url)) unless $response == RC_OK;

    my $vcf_input_fh  = Genome::Sys->open_gzip_file_for_reading($vcf_download_location);
    my ($vcf_output_fh, $vcf_temp_output) = Genome::Sys->create_temp_file();

    my @vcf_row = ();
    while (my $line = <$vcf_input_fh>) {                                                                                                                                                                    
        chomp $line;                                                                                                                                                                                  

         if ($line =~ /^#/ ){
            if($start_info && !$end_info){
                if($line !~ /^##INFO/){
                    #write the new metadata for submitter
                    print $vcf_output_fh qq(##INFO=<ID=SUB,Number=.,Type=String,Description="Name of DBSnp Submitter">\n);
                    $end_info = 1;
                }
            }elsif(!$end_info){
                $start_info = $line =~ /^##INFO/;
            }
            print $vcf_output_fh  $line . "\n"; 
            # this is a normal data line
        } else {
            my @line = split '\s+', $line;
            my $submitter_name = $self->_get_submitter_for_chromosome_and_snp_id($line[0], $line[2]);
            if($submitter_name){
                print $vcf_output_fh "$line;SUB=$submitter_name\n";
            }else{
                print $vcf_output_fh  $line . "\n";
            }
        }
    }

    $vcf_input_fh->close;
    $vcf_output_fh->close;
    my $rv = Genome::Model::Tools::Bed::ChromSort->execute(input => $vcf_temp_output, output => $self->output_file_path);
    unless ($rv) {
        $self->error_message("Failed to sort dbsnp VCF file");
        return;
    }
    return 1;
}

sub _get_submitter_for_chromosome_and_snp_id {
    my $self = shift;
    my $chromosome = shift;
    my $id = shift;

    unless($chromosome eq $self->_current_downloaded_chromosome){
        $self->_fetch_flat_file_for_chromosome($chromosome);
    }

    return $self->_submitter_map->{$id};
}

sub _fetch_flat_file_for_chromosome {
    my $self = shift;
    my $chromosome = shift;

    (my $flat_file_path = $self->flat_file_pattern) =~ s/\[X\]/$chromosome/;

    my $flat_download_location  = Genome::Sys->create_temp_file_path();

    my $response = getstore(join('/', $self->flat_url, $flat_file_path), $flat_download_location);

    die("Unable to download flat file for chromosome $chromosome at " . join('/', $self->flat_url, $flat_file_path)) unless $response == RC_OK;

    $self->_build_hash_table_for_flat_file($flat_download_location);

    $self->_current_downloaded_chromosome($chromosome);
    return 1;
}

sub _build_hash_table_for_flat_file {
    my $self = shift;
    my $flat_file_path = shift;

    #new hash table
    $self->_submitter_map({});

    my $flat_input_fh = Genome::Sys->open_gzip_file_for_reading($flat_file_path);

    my @block = ();
    while(<$flat_input_fh>){
        chomp;
        next if ($. <= 3); # each file has a 3-line header                                                                                                                                      

        my @split_line = split(/\s*\|\s*/, $_);                                                                                                                                             
        if (@split_line == 0) { # blank line                                                                                                                                                    
            my @submitters = uniq(map {$_->[1]} (sort { $b->[-1] cmp $a->[-1] } (grep { $_->[0] =~ /^ss/ } @block)));
            #add hash table entry
            $self->_submitter_map->{$block[0][0]} = join(',', @submitters);
            @block = ();                                                                                                                                                                        
        } else {                                                                                                                                                                                
            push @block, \@split_line;                                                                                                                                                          
        }                                                                                                                                                                                       
    }

    $flat_input_fh->close;
    return 1;
}

1;
