package Genome::Model::Tools::Vcf::FilterVcfInfo;

use Genome;
use strict;
use warnings;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;

class Genome::Model::Tools::Vcf::FilterVcfInfo {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "filtered VCF file",
            is_optional => 0,
        },
        vcf_file => {
            is => 'Text',
            is_input => 1,
            doc => "mutations in Vcf format",
            is_optional => 0,
        },
        filters => {
            is => 'Text',
            is_input => 1,
            doc => "comma-separated string of filters to apply. See help for syntax",
            is_optional => 0,
        },
        filter_descriptions => {
            is => 'Text',
            is_input => 1,
            doc => "comma-separated string of filters descriptions to add to the header",
            example_values => ['has GMAF>0.05,is non-nonsense mutation'],
            is_optional => 0,
        },
        non_existent_fields_are_filtered => {
            is => 'Boolean',
            is_input => 1,
            doc => "How to handle entries that do not have the parameter being filtered. If true, a filter of GMAF>0.05 will also filter entries with no GMAF defined. If false, entries with no GMAF defined will not be filtered.",
            is_optional => 1,
            default => 1,
        },

    ],
};


sub help_brief { 
    "apply filters to a VCF file based on items in the INFO field"
}


sub help_synopsis {
    <<'HELP';
    apply filters to a VCF file based on items in the INFO field
HELP
}

sub help_detail {
    <<'HELP';
    Applies filters to a VCF file based on items in the INFO field. The filters param takes a comma separated list of filters to apply. Each filter looks like the following:
    
  filtername:<infoname><operation><value>

  Examples:
    GMAFFILTER:GMAF>0.05
    TIERFILTER:TIER=2,NONSENSEFILTER:CSQ!~NONSENSE

  Supported operations are currently:
    =   equal
    !=  not equal
    >   greater than
    >=  greater than or equal to
    <   less than
    <=  less than or equal to
    ~   contains
    !~  does not contain
    

  If a line does not contain the value being searched for 

HELP
}


sub make_relational_info_filter{
    my ($relational_op, $field, $rhs, $filter_name, $filter_undef) = @_;
    
        return sub {
        my $entry = shift;
        my @common = grep {
            my $value = $entry->info_for_allele($_, $field);
            my $result;
            if($filter_undef){
                $result = defined $value && $relational_op->($value,$rhs);
            } else {
                $result = defined $value || $relational_op->($value,$rhs);
            }
            $result;

            } @{$entry->{alternate_alleles}};
        push(@common, $entry->{reference_allele});        
        $entry->filter_calls_involving_only(filter_name => $filter_name, alleles => \@common);
        return 1; # don't completely skip the entry
    }
}

sub getFilter{
    my ($filter_name, $filter_string, $filter_undef)= @_;


    my %relational_operators = (
        '>' => sub { $_[0] > $_[1]; },
        '<' => sub { $_[0] < $_[1]; },
        '<=' => sub { $_[0] <= $_[1]; },
        '>=' => sub { $_[0] >= $_[1]; },
        '=' => sub { $_[0] eq $_[1]; },
        '!=' => sub { $_[0] ne $_[1]; },
        '~' => sub { $_[0] =~ $_[1]; },
        '!~' => sub { $_[0] !~ $_[1]; },    
        );

    if($filter_string =~/([^><=!~]+)([>|<|!|=|~][=~]?)(.+)/){
        return(make_relational_info_filter($relational_operators{$2}, $1, $3, $filter_name, $filter_undef));
    } else {
        die("filter string unparsable: $filter_string\n")
    }
}


# sub wildtype_filter {
#     my $entry = shift;
#     $entry->filter_calls_involving_only(filter_name => "WILDTYPE", alleles => [$entry->{reference_allele}]);

#     return 1;
# }



################################################################################################

sub execute {
    my $self=shift;
    
    my $filter_undef = $self->non_existent_fields_are_filtered;    
    
    my @filters = split(",",$self->filters);
    my @filter_descs = split(",",$self->filter_descriptions);
    
    #sanity check
    if(@filters != @filter_descs){
        die "number of filter descriptions does not equal the number of filters\n";
    }

    # Open input vcf
    my $reader = new Genome::File::Vcf::Reader($self->vcf_file);


    # Make sure vcf header contains FT format tag and the filter names we will apply
    my $header = $reader->header;
    
    $header->add_format_type(
        id => "FT",
        type => "String",
        number => ".",
        description => "Filters",
        skip_if_exists => 1,
        );

    #loop through and apply the filters
    my $i=0;
    for($i=0;$i<@filters;$i++){
        my $filterstr = $filters[$i];
        my $filter_desc = $filter_descs[$i];
        my ($filter_name,$filter_expr) = split(":",$filterstr);

        $header->add_filter(id => $filter_name, description => $filter_desc);
        
        #create filter
        my $filter = getFilter($filter_name, $filter_expr, $filter_undef);        
        $reader->add_filter($filter);
    }


    # Create output writer
    my $writer = new Genome::File::Vcf::Writer($self->output_file, $header);

    # Copy entries from input to output
    while (my $entry = $reader->next) {
        $writer->write($entry);
    }

    return 1;
}

