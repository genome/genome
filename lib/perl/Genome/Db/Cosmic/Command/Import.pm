package Genome::Db::Cosmic::Command::Import;
use strict;
use warnings;
use Genome;

class  Genome::Db::Cosmic::Command::Import {
    is => 'Command::V2',
    has_input => [
        source_name         => { is => 'Text', is_constant => 1, value => 'cosmic' },
        database_name       => { is => 'Text', is_constant => 1, value => '' },           
        external_version    => { is => 'Text', },
        import_iteration    => { is => 'Text', },             
    ],
    has_param => [
        result_version      => { is => 'Integer', default_value => 1 },
    ],
    has_output => [
        output_dir  => {
            is => 'FilesystemPath',
            doc => 'the output directory',
        },
    ],
    #    source_directory    => { is => 'Text', is_constant => 1, value => undef },                       
    #    data_directory      => { is => 'Text', is_constant => 1, value => undef },           
    doc => 'download primary COSMIC data'
};

sub help_detail {
    return 'Import COSMIC primary data.';
}

sub import_source_url { 
    # TODO: how might we extend the API to take inputs like this?
    'ftp://ftp.sanger.ac.uk/pub/CGP/cosmic' 
}

sub repository_url {
    'https://github.com/genome-vendor/genome-db-cosmic-data.git'
}

sub execute {
    my $self = shift;
    my $external_version = $self->external_version;
    my $import_source_url = $self->import_source_url;
    my $output_dir = $self->output_dir;
    
    print "output_dir $output_dir, external_version $external_version, import_source_url $import_source_url\n";
   
    chdir($output_dir) || die "failed to chdir to $output_dir";

    my @files;
    for my $pattern (
        "CosmicCompleteExport_v${external_version}_*.tsv.gz", 
        "CosmicInsMutExport_v${external_version}_*.tsv.gz", 
    ) {
        #Genome::Sys->shellcmd(cmd => "ncftpget $import_source_url/data_export/$pattern");
        `touch CosmicCompleteExport_v65_280513.tsv CosmicInsMutExport_v65_280513.tsv; gzip -f CosmicCompleteExport*.tsv; gzip -f CosmicInsMutExport*.tsv;`;
        my @matches = glob($pattern);
        my @downloaded = grep { -e $_ } @matches;
        unless (@matches == 1 and @downloaded == 1) {
            die "no matches for @matches!";
        }
        push @files, @downloaded;
    }

    #for my $file (@files) {
    #    Genome::Sys->shellcmd(cmd => "gunzip $file");
    #}

    my $complete = $files[0];
    my $insmut = $files[1];

    my @symlink_targets = split(':',$ENV{GENOME_DB});
    my $symlink_dir = $symlink_targets[0];

    
    Genome::Sys->shellcmd(cmd => "find $output_dir");
    die "TODO: Add symlinks to the corret path!\n";
    
    return 1;
}

sub _sync_repository {

}

sub _symlink_results {

}

1;

__END__
#Originally from Cyriac Kandoth

#=================================
#Parse and cleanup COSMIC variants
ad
#=================================

cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/Cosmic_v65/

#Download the massive tab-delim file that contains curated mutation data from COSMIC, except for fusions:
wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_v65_280513.tsv.gz
gunzip CosmicCompleteExport_v65_280513.tsv.gz

#Download the tab-delim file that contains the inserted sequence for each insertion (COSMIC stores this separately for whatever reason):
wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicInsMutExport_v65_280513.tsv.gz
gunzip CosmicInsMutExport_v65_280513.tsv.gz

#There appears to be a column in CosmicCompleteExport file named "Genome-wide screen" which likely indicates whether a mutation is from a WGS sample. 
#Pull out all such sample IDs into a file:
cut -f 5,11 CosmicCompleteExport_v65_280513.tsv | grep -w y | cut -f 1 | sort -u > wgs_samples_in_cosmic_per_cosmic

#Sample IDs in wgs_samples_in_cosmic_per_cosmic appear to be more reliable:
cat wgs_samples_in_cosmic_per_cosmic | perl -ne 'chomp; @t=`grep -w $_ CosmicCompleteExport_v65_280513.tsv`; print "$_\t".scalar(@t)."\n";' > entries_per_sample.tsv

#Parse out mutations of WGS cases, and prep per-sample Build37 variant lists for WU annotation:
cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/
perl bin/parse_cosmic_for_wgs_muts.pl Cosmic_v65/CosmicCompleteExport_v65_280513.tsv Cosmic_v65/CosmicInsMutExport_v65_280513.tsv Cosmic_v65

#Run WU annotator on the variant lists of each sample:
#Remove certain categories of variants.  e.g. 'silent', '3_prime_untranslated_region', '5_prime_untranslated_region', etc.
#cat anno_files_v65/cosmic.NCBI-human.ensembl.67_37l_v2.anno | cut -f 14 | sort | uniq -c | sort -r

#NCBI-human.ensembl/67_37l_v2
cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/Cosmic_v65/
gmt annotate transcript-variants --variant-file='cosmic.var' --annotation-filter='none' --output-file='cosmic.NCBI-human.ensembl.67_37l_v2.anno'  --reference-transcripts='NCBI-human.ensembl/67_37l_v2'  --use-version='3' --nono-headers
grep -v -P "silent|3_prime_flanking_region|5_prime_flanking_region|intronic|3_prime_untranslated_region|5_prime_untranslated_region" cosmic.NCBI-human.ensembl.67_37l_v2.anno >  cosmic.NCBI-human.ensembl.67_37l_v2.anno.filt

#NCBI-human.ensembl/68_37m -> Annotator is not working with this annotation build...
cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/Cosmic_v65/
gmt annotate transcript-variants --variant-file='cosmic.var' --annotation-filter='none' --output-file='cosmic.NCBI-human.ensembl.68_37m.anno'  --reference-transcripts='NCBI-human.ensembl/68_37m'  --use-version='3' --nono-headers
grep -v -P "silent|3_prime_flanking_region|5_prime_flanking_region|intronic|3_prime_untranslated_region|5_prime_untranslated_region" cosmic.NCBI-human.ensembl.68_37m.anno > cosmic.NCBI-human.ensembl.68_37m.anno.filt

#NCBI-human.ensembl/69_37n_v3
cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/Cosmic_v65/
gmt annotate transcript-variants --variant-file='cosmic.var' --annotation-filter='none' --output-file='cosmic.NCBI-human.ensembl.69_37n_v3.anno'  --reference-transcripts='NCBI-human.ensembl/69_37n_v3'  --use-version='3' --nono-headers
grep -v -P "silent|3_prime_flanking_region|5_prime_flanking_region|intronic|3_prime_untranslated_region|5_prime_untranslated_region" cosmic.NCBI-human.ensembl.69_37n_v3.anno > cosmic.NCBI-human.ensembl.69_37n_v3.anno.filt

#NCBI-human.ensembl/70_37_v5
cd /gscmnt/sata132/techd/mgriffit/reference_annotations/Cosmic/Cosmic_v65/
gmt annotate transcript-variants --variant-file='cosmic.var' --annotation-filter='none' --output-file='cosmic.NCBI-human.ensembl.70_37_v5.anno'  --reference-transcripts='NCBI-human.ensembl/70_37_v5'  --use-version='3' --nono-headers
grep -v -P "silent|3_prime_flanking_region|5_prime_flanking_region|intronic|3_prime_untranslated_region|5_prime_untranslated_region" cosmic.NCBI-human.ensembl.70_37_v5.anno > cosmic.NCBI-human.ensembl.70_37_v5.anno.filt

