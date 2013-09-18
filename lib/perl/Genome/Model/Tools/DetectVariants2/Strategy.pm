package Genome::Model::Tools::DetectVariants2::Strategy;

use strict;
use warnings;

use Genome;
use Parse::RecDescent;

# grammar for parsing strategy rules
my $grammar = q{
    startrule: combination end
        { $item[1]; }
    | <error>

    end: /^\Z/

    combination: intersection
                { $item[1]; }
    | union
                { $item[1]; }
    | unique_union 
                { $item[1]; }
    | single
                { $item[1]; }
    | <error>
    
    parenthetical: "(" combination ")"
                { $item[2]; }
                
    intersection: single "intersect" combination
                { $return = { intersect => [$item[1], $item[3] ] }; }
    
    union: single "union" combination
                { $return = { union => [ $item[1], $item[3] ] }; }

    unique_union: single "unique union" combination
                { $return = { unionunique => [ $item[1], $item[3] ] }; }
    
    single: parenthetical
                { $item[1]; }
    | strategy
                { $item[1]; }
    
    strategy: program_spec "filtered by" filter_list
                { $return = { detector => {%{$item[1]}, filters => $item[3]} }; }
    | program_spec 
                { $return = { detector => {%{$item[1]}, filters => []} }; }
    | <error>

    filter_list: program_spec "then" filter_list
                { $return = [$item[1], @{$item[3]}]; }
    | program_spec
                { $return = [$item[1]]; }

    word: /([\w\.:-]|\\\\)+/ { $return = $item[1]; }

    valid_subpackage: "somatic "
                { $return = $item[1]; }

    name: valid_subpackage word
                { $return = "$item[1] $item[2]"; }
    | word
                { $return = $item[1]; }
    | <error>

    version: word { $return = $item[1]; }
    | <error>

    params: {
                my $txt = extract_bracketed($text, '[]');
                $txt =~ s/^\[(.*)\]$/$1/;
                $txt=~ s/\\\\([\[\]])/$1/g;
                $return = $txt;
            } 

    program_spec: name version params
                { $return = {
                    name => $item[1],
                    version => $item[2],
                    params => $item[3],
                    };
                }
    | name version
                { $return = {
                    name => $item[1],
                    version => $item[2],
                    params => '',
                    };
                }
};


class Genome::Model::Tools::DetectVariants2::Strategy {
    is => ['UR::Value'],
    id_by => 'id',
    has => [
        id => { is => 'Text' },
    ],
    doc => 'This class represents a variant detection strategy. It specifies a variant detector, its version and parameters, as well as any filtering to be done.'
};

sub __errors__ {
    my $self = shift;
    my @tags = $self->SUPER::__errors__(@_);

    my $tree = $self->tree;
    push @tags, $self->{_last_error} if exists $self->{_last_error};
    unless ($tree) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['id'],
            desc => "Failed to create detector strategy from id string " . $self->id,
            );
    }

    eval {
        $self->validate($tree);
    };
    if ($@) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['id'],
            desc => $@
            );
    }

    return @tags;
}

sub create {
    my $class = shift;
    return $class->SUPER::create(@_);
}

sub parse {
    my ($self, $str) = @_;
    my $parser = Parse::RecDescent->new($grammar)
        or die "Failed to create parser from grammar";
    my $tree = $parser->startrule($str);
    $self->_add_class_info($tree) if $tree;
    return $tree;
}

sub _set_last_error {
    my ($self, $msg) = @_;
    $self->{_last_error} = UR::Object::Tag->create(
        type => 'error',
        properties => [],
        desc => $msg,
    );
}

sub tree {
    my $self = shift;

    # already computed it
    return $self->{_tree} if defined $self->{_tree};

    eval { $self->{_tree} = $self->parse($self->id); };
    if ($@) {
        $self->_set_last_error($@);
        return;
    }

    return $self->{_tree};
}

sub get_detectors {
    my $self = shift;

    my $tree = $self->parse($self->id);
    return $self->_get_detectors($tree);
}

sub _get_detectors {
    my $self = shift;
    my $tree = shift;

    my @detectors;
    if (ref($tree) eq 'ARRAY') {
        for my $sub_tree (@$tree) {
            my @other_detectors = $self->_get_detectors($sub_tree);
            if (@other_detectors) {
                push @detectors, @other_detectors;
            }
        }
    } elsif (ref($tree) eq 'HASH') {
        for my $key (keys(%$tree)) {
            if ($key eq 'detector') {
                push @detectors, $tree->{$key};
            } else {
                my @other_detectors = $self->_get_detectors($tree->{$key});
                if (@other_detectors) {
                    push @detectors, @other_detectors;
                }
            }
        }
    }
    return @detectors;
}

sub _add_class_info {
    my ($self, $tree) = @_;
    my @keys = keys %$tree;
    for my $key (@keys) {
        if ($key eq 'detector') {
            my $name = $tree->{$key}->{name};
            my $class = $self->detector_class($name);
            $tree->{$key}->{class} = $class;
            # FIXME combine this code with the below
            my $number_of_filters = scalar @{$tree->{$key}->{filters}};
            next unless $number_of_filters > 0;
            for my $index (0..$number_of_filters - 1) {
                my $name = @{$tree->{$key}->{filters}}[$index]->{name};
                my $class = $self->filter_class($name);
                @{$tree->{$key}->{filters}}[$index]->{class} = $class;
            }
        } elsif (ref $tree->{$key} eq 'ARRAY') {
            ($self->_add_class_info($_)) for (@{$tree->{$key}});
        }
    }
}

sub detector_class {
    my $self = shift;
    my $detector = shift;
    
    # Convert things like "hi foo-bar" to "Hi::FooBar"
    $detector = join("::", 
        map { join('', map { ucfirst(lc($_)) } split(/-/, $_))
            } split(' ', $detector));
    
    my $detector_class_base = 'Genome::Model::Tools::DetectVariants2';
    my $detector_class = join('::', ($detector_class_base, $detector));
    
    return $detector_class;
}

# TODO merge this with the above method, probably. The only difference is the ::Filter part.
sub filter_class {
    my $self = shift;
    my $filter = shift;
    
    # Convert things like "hi foo-bar" to "Hi::FooBar"
    $filter = join("::", 
        map { join('', map { ucfirst(lc($_)) } split(/-/, $_))
            } split(' ', $filter));
    
    my $filter_class_base = 'Genome::Model::Tools::DetectVariants2::Filter';
    my $filter_class = join('::', ($filter_class_base, $filter));
    
    return $filter_class;
}

sub validate {
    my ($self, $tree) = @_;
    #return;
    #$DB::single=1;

    my @keys = keys %$tree;
    for my $key (@keys) {
        if ($key eq 'detector') {
            ## Check detector class and params
            my $class = $tree->{$key}->{class};
            my $version = $tree->{$key}->{version};
            unless($self->is_class_valid($class)){
                $self->error_message("could not find class ".$class);
                die $self->error_message;
            }
            unless($self->has_version($class,$version)){
                $self->error_message("could not find version ".$version." for class ".$class);
                die $self->error_message;
            }

            ## Check filters that belong to this detector
            for my $filter (@{$tree->{$key}->{filters}}) {
                my $filter_class = $filter->{class};
                my $filter_version = $filter->{version};
                unless($self->is_class_valid($filter_class)){
                    $self->error_message("could not find class ".$filter_class);
                    die $self->error_message;
                }
                unless($self->has_version($filter_class,$filter_version)){
                    $self->error_message("could not find version ".$filter_version." for class ".$filter_class);
                    die $self->error_message;
                }
            }

        } elsif (ref $tree->{$key} eq 'ARRAY') {
            for my $subtree ( @{$tree->{$key}} ) {
                $self->validate($subtree);
            }
        }
    }
}

sub is_class_valid {
    my $self = shift;
    my $class = shift;
    my $result;
    eval {
        $result = $class->__meta__;
    };
    if($@){
        return 0;
    }
    return 1;
}

sub has_version {
    my $self = shift;
    my ($class,$version) = @_;
    my $answer = 0;
    eval {
        $answer = $class->has_version($version);
    };
    if($@){
        $self->error_message("Could not call has_version on the class ".$class);
        die $self->error_message;
    }
    return $answer;
}

1;
__END__

=pod

=head1 NAME

Genome::Model::Tools::DetectVariants2::Strategy

=head1 SYNOPSIS

    # Detect snvs with sniper version 0.7.3 with the parameters "-q 1 -Q 15".
    'sniper 0.7.3 [-q 1 -Q 15]'

    # Detect snvs with sniper version 0.7.3 with the listed parameters and filter the results by running the "loh" filter version "v1".
    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by loh v1 '

    # Detect snvs: 
    # 1) Run sniper version 0.7.3 with parameters
    # 2) Filter the results by running the loh filter version v1, i
    # 3) Further filter results and then the somatic-score-mapping-quality filter version v1 with parameters.
    # 4) Run samtools version r599 (or steal previous results) 
    # 5) Intersect 3 & 4 
    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by loh v1 , somatic-score-mapping-quality v1 [-min_somatic_quality 40:-min_mapping_quality 40] intersect samtools r599'  
    
    # Detectd indels with: 
    # 1) Run sniper version 0.7.3 with the listed parameters. 
    # 2) Run samtools version r599 
    # 3) Run pindel version v1
    # 4) Intersect 2 and 3
    # 5) Union 1 and 4.
    'sniper 0.7.3 [ -q 1 -Q 15 ] union (samtools r599  intersect pindel v1 )'

    # Detect snvs or indels or both with sniper version 0.7.3 with the listed parameters. 
    # This expression can be set as an snv detection strategy or an indel detection strategy, 
    # and if both are set to the same value sniper will run just once to do both.
    'sniper 0.7.3 [ -q 1 -Q 15 ]' 
    
    # Detect structural variation with breakdancer version 2010_06_24.
    'breakdancer 2010_06_24 ' 

    # Detect snvs: Intersect the results of sniper version 0.7.3 with parameters and samtools version r599.
    'sniper 0.7.3 [ -q 1 -Q 15 ] intersect samtools r599 '
    
    # Detect indels using sniper version 0.7.3 with parameters and filter the results with the library-support filter version v1
    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by library-support v1 ' 
    
    # Detect structural variations using breakdancer version 2010_06_24 and filter the results by applying the tigra-validation filter version v1
    'breakdancer 2010_06_24  filtered by tigra-validation v1 '

    
    # Detect indels using sniper version 0.7.3 with parameters and filter the results with the library-support filter version v1
    'sniper 0.7.3 [ -q 1 -Q 15 ] filtered by library-support v1 ' 
    
    # Detect structural variations using breakdancer version 2010_06_24 and filter the results by applying the tigra-validation filter version v1
    'breakdancer 2010_06_24  filtered by tigra-validation v1 '

=head1 DESCRIPTION

=head2 COMPONENTS

=over 4

    A strategy consists of the following:
    detector-name version [ params ] filtered by filter-name version [ params ],filter-name version [ params ] ...

    * Detector-name is the name of the variant detector as it follows "gmt detect-variants2". For example, "sniper" would reference the tool located at "gmt detect-variants2 sniper".

    * In the same way, filter-name is the name of the filter as it follows "gmt detect-variants2 filter". For example, "loh" would reference the tool located at gmt detect-variants2 filter loh".

    * Version is a version number that pertains to that detector or filter specifically. For sniper this might be "0.7.3". For samtools this might be "r599".
        Many filters are not currently versioned, but may be in the future. In these cases "v1" should be used to denote version 1.

    * The parameter list is a list of all parameters to be passed to the detector or filter and will be specific to that tool. It is passed as a single string and is optional.

    * Filtered by may contain any number of complete filter specifications (separated by commas), including 0. Each filter must be a complete list of name, version, and an optional param list.

=back

=head2 UNIONS AND INTERSECTIONS

=over 4

    * Variant detectors can be intersected or unioned with each other to create variant lists which utilize more than one variant detector. In either case, all variant detectors will be run individually and then processed together.
    * An intersection will run both detectors and then produce a final list of variants that represents where both detectors agree on both the position and the call.
    * A union will run both detectors and then produce a final list of variants that represents every call that both the detectors made, regardless of agreement.
    * Parenthesis may also be utilized around pairs of detectors to specify logical order of operation.

    --- Examples of union and intersection --- 
    --snv-detection-strategy 'sniper 0.7.3 [-q 1 -Q 15] intersect samtools r599
    This represents the desire to run version 0.7.3 of sniper with the above parameter list and version r599 of samtools with no parameters and intersect the results. 
    Both detectors will be run and the final variant list will represent all variants which were called by both detectors.

    --snv-detection-strategy 'sniper 0.7.3 [-q 1 -Q 15] union (samtools r599 intersect pindel v1)
    This represents the desire to run version 0.7.3 of sniper with the above parameters, version r599 of samtools with no parameters, and version v1 of pindel with no parameters.
    Due to the parenthesis, the results of pindel and samtools will first be intersected and then that result will be unioned with the variant calls from sniper.
    In plain language, the resulting set will be any variants that either a) sniper called or b) pindel and samtools both called and agreed on.

=back

=cut

