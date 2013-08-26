package Genome::Model::Tools::Dgidb::QueryGene;

use strict;
use warnings;

use JSON;
use Genome;
use LWP::UserAgent;
use HTTP::Request::Common;

my $DOMAIN   = 'http://dgidb.genome.wustl.edu/';
my $API_PATH = '/api/v1/interactions.json';

my %OPTIONAL_PROPERTIES = (
    interaction_sources => {
        is  => 'Text',
        doc => 'Limit results to those from particular data sources. e.g. DrugBank, PharmGKB, TALC, TEND, TTD, MyCancerGenome',
    },
    interaction_types => {
        is  => 'Text',
        doc => 'Limit results to interactions with drugs that have a particular mechanism of action. e.g. inhibitor, antibody, etc.',
    },
    gene_categories => {
        is  => 'Text',
        doc => 'Limit results to genes with a particular druggable gene type. e.g. KINASE, ION CHANNEL, etc.',
    },
    source_trust_levels => {
        is  => 'Text',
        doc => 'Limit results based on trust level of the interaction source. e.g. \'Expert curated\' or Non-curated',
    },
);

class Genome::Model::Tools::Dgidb::QueryGene {
    is  => 'Command::V2',
    has => [
        genes => {
            is  => 'Text',
            doc => 'List of gene symbols. Use offical entrez symbols for best results. Separated by comma'
        },
    ],
    has_optional => [
        %OPTIONAL_PROPERTIES,
        output_file => {
            is  => 'Text',
            doc => 'A file path to store the output. Default is to STDOUT',
        },
        antineoplastic_only => {
            is  => 'Boolean',
            doc => 'Limit results to anti-cancer drugs only',
            default_value => 0,
        },
    ],
    has_optional_output => [
        output_hash_ref => {
            is  => 'HASH',
        },
    ],
};


sub help_brief {
    'Tool to query genes from DGIDB database';
}

sub help_detail {
    return <<EOS
    Tool to query genes from DGIDB database with different criteria. The example commands are:
    gmt dgidb query-gene --genes=FLT3 --output-file=myDGIDB.out
    gmt dgidb query-gene --genes=FLT3,EGFR,KRAS
    gmt dgidb query-gene --genes=FLT3,EGFR --interaction-sources=TALC,TEND
    gmt dgidb query-gene --genes=FLT3,EGFR --gene-categories=KINASE
    gmt dgidb query-gene --genes=FLT3,EGFR --interaction-types=inhibitor
    gmt dgidb query-gene --genes=FLT3,EGFR --source-trust-levels=\'Expert curated\'
    gmt dgidb query-gene --genes=FLT3,EGFR --antineoplastic-only
    gmt dgidb query-gene --genes=FLT3,EGFR,KRAS --interaction-sources=TALC,TEND,MyCancerGenome --gene-categories=KINASE --interaction-types=inhibitor --antineoplastic-only --output-file=myTest.out
EOS
}

sub execute {
    my $self = shift;

    my %params = (genes => $self->genes);
    $params{drug_types} = 'antineoplastic' if $self->antineoplastic_only;

    for my $property (keys %OPTIONAL_PROPERTIES) {
        $params{$property} = $self->$property if $self->$property;
    }

    my $ua   = LWP::UserAgent->new;
    my $resp = $ua->request(POST $DOMAIN . $API_PATH, \%params);

    if ($resp->is_success) {
        $self->write_output(decode_json($resp->content));
    } 
    else {
        die $self->error_message("Something went wrong! Did you specify any genes?\n");
    }

    return 1;
}


sub write_output {
    my ($self, $json_ref) = @_;
    my $out_file = $self->output_file;

    my $writer;
    my @headers = qw(gene_name drug_name interaction_type source gene_categories);

    my %params = (
        headers   => \@headers,
        separator => "\t",
    );
    $params{output} = $out_file if $out_file; #if no outfile, default to STDOUT

    $writer = Genome::Utility::IO::SeparatedValueWriter->create(%params);
    die $self->error_message("Failed to create IO SeparatedValueWriter for $out_file") unless $writer;
    
    my $output = {};

    #loop over each search term that was definitely matched
    for my $matched (@{$json_ref->{matchedTerms}}) {
        my $gene_categories = join ",", sort @{$matched->{geneCategories}};
        $gene_categories    = $gene_categories ? lc $gene_categories : 'n/a';

        #loop over any interactions for this gene
        for my $interaction (@{$matched->{interactions}}) {
            my %content;
            @content{@headers} = (
                $matched->{geneName}, 
                $interaction->{drugName}, 
                $interaction->{interactionType},
                $interaction->{source},
                $gene_categories,
            );

            $writer->write_one(\%content);
            push @{$output->{matchedTerms}}, \%content;
        }
    }

    #loop over each search term that wasn't matched definitely
    for my $unmatched (@{$json_ref->{unmatchedTerms}}) {
        my %content = (
            searchTerm  => $unmatched->{searchTerm},
            suggestions => $unmatched->{suggestions},
        );
        push @{$output->{unmatchedTerms}}, \%content;
    }

    $self->output_hash_ref($output);
    return 1;
}

1;
