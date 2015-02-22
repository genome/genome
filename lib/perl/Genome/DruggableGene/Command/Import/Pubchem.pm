package Genome::DruggableGene::Command::Import::Pubchem;

use strict;
use warnings;

use Genome;
use IO::File;
use Data::UUID;
use POSIX qw/strftime/;
use DBI;

class Genome::DruggableGene::Command::Import::Pubchem {
    is => 'Command::V2',
    has_optional => [
        postgres_host => {
            is => 'text',
            doc => 'The hostname of the postgres box to which data will be imported',
            default => 'vmpool50',
        },
        synonyms_file => {
            is => 'text',
            doc => 'Path to synonyms file, optionally unzipped, downloaded from ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz',
            default => '',
        },
        scratch_dir => {
            is => 'path',
            doc => 'Path to scratch dir where we have room to split synonyms file into parts',
            default => '/tmp/pubchem_import/',
        },
        cleanup => {
            is => 'boolean',
            doc => 'delete generated files',
            default => 1,
        },
        skip_db_copy => {
            is => 'boolean',
            doc => 'Dont touch db, just generate tsv',
            default => 0,
        },
        generate_uuids_locally => {
            is => 'boolean',
            doc => 'Dont rely on postgres to generate UUIDs with its uuid-ossp extension',
            default => 0,
        },
        citation_base_url => {
            default => 'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=',
        },
        citation_site_url => {
            default => 'http://pubchem.ncbi.nlm.nih.gov/',
        },
        citation_text => {
            default => "PubChem's BioAssay Database. Wang Y, Xiao J, ..., Gindulyte A, Bryant SH. Nucleic Acids Res. 2012 Jan;40(Database issue):D400-12. Epub 2011 Dec 2. PMID: 22140110",
        },
    ],
    doc => 'Import synonyms from pubchem into a postgres db',
};

sub help_brief { 'Import synonyms from pubchem into a postgres db' }

sub help_detail { help_brief() }

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);
  my $test_psql = `which psql`;
  chomp($test_psql);

  unless ($test_psql) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
                                            properties => ['psql'],
	                                          desc => "psql command not found on this machine!",
                                          );
  }
  return @errors;
}

sub execute {
    my $self = shift;

    use_scratch($self->scratch_dir);
    download_unzip_split($self->synonyms_file, $self->cleanup);
    $self->write_tsvs($self->cleanup, $self->generate_uuids_locally, get_hash_of_existing_drugs($self->postgres_host));
    copy_into_postgres($self->postgres_host, $self->cleanup) unless $self->skip_db_copy;
    print_timestamp('Done!!!');

    1;
}

sub use_scratch {
    my $scratch_dir = shift;
    `mkdir -p $scratch_dir`;
    die "Unable to cd to $scratch_dir" unless chdir $scratch_dir;
    print_timestamp("Now working within $scratch_dir");
}

sub download_unzip_split {
    #Synonyms file can be the gzipped download or an unzipped file, in any directory
    #It will only be deleted (cleaned up) if it is in the scratch_dir
    my $synonyms_file = shift;
    my $cleanup = shift;

    if (-e 'pubchem_a' or -e 'drug_name_report_pubchem.tsv' or -e 'drug_name_report_association_pubchem.tsv' or -e 'citation_pubchem.tsv'){
        print_timestamp('Skipping download_unzip_split, found file splits or tsvs');
        return;
    }

    unless ($synonyms_file or -e 'CID-Synonym-filtered' or -e 'CID-Synonym-filtered.gz'){ #Do we already have it downloaded?
        print_timestamp('Downloading ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz');
        `wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz`;
    }
    $synonyms_file = 'CID-Synonym-filtered.gz' if not $synonyms_file and -e 'CID-Synonym-filtered.gz' and not -e 'CID-Synonym-filtered';
    $synonyms_file = 'CID-Synonym-filtered' unless $synonyms_file;


    if ($synonyms_file =~ /gz$/ or `file -- $synonyms_file` =~ /gzip/){ #Is it already unzipped?
        print_timestamp('Unzipping file');
        `gunzip -c $synonyms_file > CID-Synonym-filtered`;

        print_timestamp('Removing gzip file if it exists in scratch_dir') if $cleanup;
        `rm -f CID-Synonym-filtered.gz` if $cleanup;
    }
    $synonyms_file = 'CID-Synonym-filtered';

    print_timestamp("Splitting $synonyms_file into parts of 10 million lines each");
    `split -a1 -l10000000 $synonyms_file pubchem_`;

    print_timestamp('Removing unzipped synonyms file if it exists in scratch_dir') if $cleanup;
    `rm -f CID-Synonym-filtered` if $cleanup;
}

sub get_hash_of_existing_drugs {
    my $postgres_host = shift;
    my %drugs;

    my $dbh = DBI->connect("dbi:Pg:dbname=genome;host=$postgres_host;",'genome','TGI_pg_1');
    map{$drugs{pop @$_}++} @{$dbh->selectall_arrayref("select alternate_name from drug_name_report_association")};

    return %drugs;
}

sub write_tsvs {
    my $self = shift;
    my $cleanup = shift;
    my $generate_uuids_locally = shift;
    my %existing_drugs = @_;

    my $synonym_file_lines = 72000000; #hard coded to save time calculating for its only used in progress counter and won't affect generated data
    my $uuid_gen = new Data::UUID;

    if (-e 'drug_name_report_pubchem.tsv' and -e 'drug_name_report_association_pubchem.tsv' and -e 'citation_pubchem.tsv'){
        print_timestamp('Skipping tsv generation, tsvs are already created');
        return;
    }
    if (-e 'drug_name_report_pubchem.tsv' or -e 'drug_name_report_association_pubchem.tsv' or -e 'citation_pubchem.tsv'){
        print_timestamp('One or more tsvs are missing, deleting the left-overs');
        `rm -f drug_name_report_pubchem.tsv drug_name_report_association_pubchem.tsv citation_pubchem.tsv`;
    }

    my $citation_tsv = IO::File->new('citation_pubchem.tsv', 'w');
    my $drug_name_report_tsv = IO::File->new('drug_name_report_pubchem.tsv', 'w');
    my $drug_name_report_association_tsv = IO::File->new('drug_name_report_association_pubchem.tsv', 'w');

    my ($current_line, $cur_pubchem_id) = (0,0);
    my @cur_cid_drugs;

    print_timestamp('Creating citation sql');
    my $citation_id = uuid($uuid_gen);
    my $version = strftime('%d-%b-%Y', localtime);
    my $citation_text = $self->citation_text;
    my $base_url = $self->citation_base_url;
    my $site_url = $self->citation_site_url;
    print $citation_tsv "$citation_id\tPubChem\t$version\t$citation_text\t$base_url\t$site_url\tPubChem\n";

    for my $file (glob 'pubchem_*'){
        print_timestamp("Slurping file $file");
        my $fh = IO::File->new($file, 'r');
        my @lines = <$fh>;
        print_timestamp('Generating sql');
        for my $line (@lines){
            $current_line++;

            my ($id, $name) = split("\t",$line);
            chomp $name;
            $name = uc $name;

            if ($cur_pubchem_id != $id){
                for my $pubchem_drug (@cur_cid_drugs) {
                    if ($existing_drugs{$pubchem_drug}) {

                        my $cur_drug_id = uuid($uuid_gen);
                        print $drug_name_report_tsv "$cur_drug_id\t$cur_pubchem_id\t\\N\tpubchem\t$citation_id\n";

                        my $pubchem_primary_name = shift @cur_cid_drugs;
                        my $alt_name_id = $generate_uuids_locally ? uuid($uuid_gen) : '\\N';
                        print $drug_name_report_association_tsv "$alt_name_id\t$cur_drug_id\t$pubchem_primary_name\t1\tpubchem_primary_name\n";

                        my $popularity = 2;
                        for my $pubchem_alt_name (@cur_cid_drugs){
                            my $alt_name_id = $generate_uuids_locally ? uuid($uuid_gen) : '\\N';
                            print $drug_name_report_association_tsv "$alt_name_id\t$cur_drug_id\t$pubchem_alt_name\t$popularity\tpubchem_alt_name\n";
                            $popularity++;
                        }

                        last;
                    }
                }

                $cur_pubchem_id = $id;
                @cur_cid_drugs = ();
            }

            push @cur_cid_drugs, $name;

            if($current_line % 1000000 == 0){
                print_timestamp(sprintf "$current_line out of $synonym_file_lines , %.2f%%", $current_line / $synonym_file_lines * 100);
            }
        }
    }

    $citation_tsv->close;
    $drug_name_report_tsv->close;
    $drug_name_report_association_tsv->close;


    print_timestamp('Removing synonyms file parts') if $cleanup;
    `rm -f pubchem_*` if $cleanup;
}

sub copy_into_postgres {
    my $postgres_host = shift;
    my $cleanup = shift;

    print_timestamp("Setting user/pass for postgres on $postgres_host in ~/.pgpass");
    `echo $postgres_host:5432:genome:genome:TGI_pg_1 >> ~/.pgpass`;
    `chmod 600 ~/.pgpass`;
    print_timestamp('Copying citation_pubchem.tsv into citation table');
    `psql -h $postgres_host -U genome -c '\\copy citation from citation_pubchem.tsv'`;
    print_timestamp('Copying drug_name_report_pubchem.tsv into drug_name_report table');
    `psql -h $postgres_host -U genome -c '\\copy drug_name_report from drug_name_report_pubchem.tsv'`;
    print_timestamp('Copying drug_name_report_association_pubchem.tsv into drug_name_report_association table');
    `psql -h $postgres_host -U genome -c '\\copy drug_name_report_association from drug_name_report_association_pubchem.tsv'`;

    print_timestamp('Removing tsvs') if $cleanup;
    `rm -f drug_name_report_pubchem.tsv drug_name_report_association_pubchem.tsv citation_pubchem.tsv` if $cleanup;
}

sub uuid {
    my $uuid_gen = shift;
    my $a = $uuid_gen->to_string($uuid_gen->create());
    return $a;
}

sub print_timestamp {
    my $data = shift;
    my ($sec,$min,$hour) = localtime time;
    printf "$hour:$min:%02d", $sec;
    print "  $data\n";
}

1;
