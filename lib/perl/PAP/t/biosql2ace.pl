#!/usr/bin/env genome-perl


use strict;
use warnings;

use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;
use Bio::Seq;

use Data::Dumper;
use English;
use Getopt::Long;
use IO::File;


my $dev_flag = 0;

GetOptions(
           "--dev" => \$dev_flag,
       );


my $product_fh = IO::File->new();
my $product_fn = 'productID.ace';

$product_fh->open(">$product_fn") or die "Can't open '$product_fn': $OS_ERROR";

my $psortb_fh = IO::File->new();
my $psortb_fn = 'psort.ace';

$psortb_fh->open(">$psortb_fn") or die "Can't open '$psortb_fn': $OS_ERROR";

my $blastp_fh = IO::File->new();
my $blastp_fn = 'briefID.fof.ace';

$blastp_fh->open(">$blastp_fn") or die "Can't open '$blastp_fn': $OS_ERROR";

my $interpro_fh = IO::File->new();
my $interpro_fn = 'merged.raw.sorted.ace';

$interpro_fh->open(">$interpro_fn") or die "Can't open '$interpro_fn': $OS_ERROR";

my $interpro_ipr_fh = IO::File->new();
my $interpro_ipr_fn = 'merged.raw.sorted.ipr.ace';

$interpro_ipr_fh->open(">$interpro_ipr_fn") or die "Can't open '$interpro_ipr_fn': $OS_ERROR";

my $kegg_fh = IO::File->new();
my $kegg_fn = 'REPORT-top_new.ks.ace';

$kegg_fh->open(">$kegg_fn") or die "Can't open '$kegg_fn': $OS_ERROR";


## Need a unique list of IPR numbers per gene prediction
my %ipr = ( );


my $dbadp;

if ($dev_flag) {
    
    $dbadp = Bio::DB::BioDB->new(
                                 -database => 'biosql',
                                 -user     => 'sg_user',
                                 -pass     => 'sgus3r',
                                 -dbname   => 'DWDEV',
                                 -driver   => 'Oracle'
                             );
    
}
else {
    
    $dbadp = Bio::DB::BioDB->new(
                                 -database => 'biosql',
                                 -user     => 'sg_user',
                                 -pass     => 'sg_us3r',
                                 -dbname   => 'DWRAC',
                                 -driver   => 'Oracle'
                             );
    
}

my $adp = $dbadp->get_object_adaptor("Bio::SeqI");

my $query = Bio::DB::Query::BioQuery->new();

$query->datacollections([
                         "Bio::PrimarySeqI s",
                        ]);

$query->where(["s.display_id like '$ARGV[0]%'"]);

my $result = $adp->find_by_query($query);

while (my $seq = $result->next_object()) {

    my @features = $seq->get_SeqFeatures();

    FEATURE: foreach my $feature (@features) {

        ## Only the protein coding_genes have protein annotation to dump.
        ## The 'InterPro' features are special in that they have their 
        ## source_tag set to something other than 'gene', but we want them,
        ## too.  Perhaps PAP::Command::InterProScan should just set the 
        ## source to gene...
        unless (
                ($feature->primary_tag() eq 'gene')
                ||
                ($feature->source_tag() eq 'InterPro')
               ) { 
            next FEATURE; 
        }
        
        my $display_name = $feature->display_name();
        
        my $new_display_name;

        if (
            ($display_name =~ /^(\w+\d+)\.(\w+)\.(\d+)$/) 
            ||
            ($display_name =~ /^(\w+\d+)\.(\w+)\.(\d+)\.InterPro\.\d+$/) 
            ||
            ($display_name =~ /^(\w+\d+\.\d+)\.(\w+)\.(\d+)$/)
            ||
            ($display_name =~ /^(\w+\d+\.\d+)\.(\w+)\.(\d+)\.InterPro\.\d+$/)
           ) {
            my ($seq_id, $source, $number) = ($1, $2, $3);
            $new_display_name = join('.', $seq_id, $source, 'p5_hybrid', $number);
        }
        else {
            die "failed to parse '$display_name':\n - does not match expected format (seqid.predictor.sequence)";
        }

        
        ## InterPro
        ##
        ## This thing is not like the others.  Whereas all
        ## the other annotation got shoved into the features
        ## for the predicted genes, due to the need to have
        ## multiple sets of tags and multiple dblinks and
        ## have them properly linked, InterPro gets its own
        ## features. 
        if ($feature->source_tag eq 'InterPro') {

            $new_display_name =~ s/\.InterPro\.\d+//g;
            
            my ($dblink) = grep { $_->database() eq 'InterPro' } $feature->annotation->get_Annotations();

            if (defined($dblink)) {

                my $ipr_number = $dblink->primary_id();

                $ipr{$new_display_name}{$ipr_number} = 1;
                
                my ($analysis) = $feature->each_tag_value('interpro_analysis');
                my ($evalue)   = $feature->each_tag_value('interpro_evalue');
                my ($desc)     = $feature->each_tag_value('interpro_description');

                ## Bug for bug replication...
                if ($evalue == 1e10) { $evalue = ''; }
                
                print $interpro_fh qq{Sequence $new_display_name}, "\n";
                print $interpro_fh qq{Interpro   "$analysis : $ipr_number $desc : pval $evalue"}, "\n\n";

            }

            next FEATURE;

        }
        
        ## Psort-B
        if (
            $feature->has_tag('psort_localization') &&
            $feature->has_tag('psort_score')
        ) {
            
            my ($psort_localization) = $feature->each_tag_value('psort_localization');
            my ($psort_score)        = $feature->each_tag_value('psort_score');

            print $psortb_fh qq{Sequence $new_display_name}, "\n";
            print $psortb_fh qq{PSORT_B $psort_localization $psort_score}, "\n\n";

        }

        
        ## BlastP Product Categorization
        if ($feature->has_tag('blastp_category')) {
        
            my ($product) = $feature->each_tag_value('blastp_category');

            print $product_fh qq{Sequence $new_display_name}, "\n";
            print $product_fh qq{Product_ID "$product"}, "\n\n";
            
        }
        
        my $annotation_collection = $feature->annotation();
        
        my @annotations = $annotation_collection->get_Annotations();
        my @dblinks     = grep { $_->isa('Bio::Annotation::DBLink') } @annotations;

        
        ## BlastP
        if (@dblinks > 0) {
            
            my ($top_hit) = grep { $_->database() eq 'GenBank' } @dblinks;

            if (defined($top_hit)) {
                
                my $hit_accession = $top_hit->primary_id();
                
                my ($bit_score)  = $feature->each_tag_value('blastp_bit_score');
                my ($evalue)     = $feature->each_tag_value('blastp_evalue');
                my ($percent_id) = $feature->each_tag_value('blastp_percent_identical');
                my ($q_start)    = $feature->each_tag_value('blastp_query_start');
                my ($q_end)      = $feature->each_tag_value('blastp_query_end');
                my ($s_start)    = $feature->each_tag_value('blastp_subject_start');
                my ($s_end)      = $feature->each_tag_value('blastp_subject_end');
                my ($hit_name)   = $feature->each_tag_value('blastp_hit_name');
                my ($hit_desc)   = $feature->each_tag_value('blastp_hit_description');
                
                my $brief_id = join(
                                    "\t",
                                    qq{Brief_identification "Stats: $bit_score},
                                    $evalue,
                                    $percent_id,
                                    $q_start,
                                    $q_end,
                                    $s_start,
                                    $s_end,
                                    qq{ TopHit: $hit_name $hit_desc "},
                                );
                
                print $blastp_fh qq{Sequence $new_display_name}, "\n";
                print $blastp_fh qq{$brief_id}, "\n\n";
                
            }
            
        }


        ## KEGG
        if (@dblinks > 0) {

            my @kegg_dblinks = grep { $_->database() eq 'KEGG' } @dblinks;
            
            my ($gene_dblink)      = grep { $_->primary_id() =~ /^\w{3}\:\w+/ } @kegg_dblinks;
            my ($orthology_dblink) = grep { $_->primary_id() =~ /^K\d+(\:K\d+)*$/   } @kegg_dblinks;

            if (
                $feature->has_tag('kegg_evalue') &&
                $feature->has_tag('kegg_description') &&
                defined($gene_dblink)
            ) {
                
                my $gene_id            = $gene_dblink->primary_id();
               
                my ($kegg_evalue)      = $feature->each_tag_value('kegg_evalue');
                my ($kegg_description) = $feature->each_tag_value('kegg_description');

                my $kegg_ace = qq{KEGG   "$gene_id $kegg_evalue $kegg_description};

                if (defined($orthology_dblink)) {
                    my $orthology_id       = $orthology_dblink->primary_id();    
                    $kegg_ace = join('  ', $kegg_ace, $orthology_id);
                }

                $kegg_ace = join('', $kegg_ace, qq{"});
                
                print $kegg_fh qq{Sequence $new_display_name}, "\n";
                print $kegg_fh $kegg_ace, "\n\n";
                
            }
            
        }
        
    }
    
}

foreach my $gene (sort keys %ipr) {

    my $ipr_string = join " ", (sort { $a cmp $b } keys %{$ipr{$gene}});

    print $interpro_ipr_fh qq{Sequence $gene}, "\n";
    print $interpro_ipr_fh qq{IPR_ID "$ipr_string"}, "\n\n";

}
