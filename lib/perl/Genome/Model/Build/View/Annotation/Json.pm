package Genome::Model::Build::View::Annotation::Json;

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Genome;

use JSON;
use IO::File;

class Genome::Model::Build::View::Annotation::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        perspective => {
            value => 'annotation',
        },
    ],
    has => [
        datatables_params => {
            is => 'String',
        },
        request_index => {
            is => 'String',
        },
        request_grep => {
            is => 'String',
        },
        request_grab => {
            is => 'String',
        },
    ],
};

sub check_build {
    my $build = shift;
    my $name = $build->__display_name__;
#    if (!$build->can("snvs_bed") or !defined $build->snvs_bed("v2")) {
#        die "Failed to get snvs for build $name (snvs_bed missing or returned null)\n";
#    }
#
#    if (!$build->can("filtered_snvs_bed") or !defined $build->filtered_snvs_bed("v2")) {
#        die "Failed to get filtered snvs for build $name (filtered_snvs_bed missing or returned null)\n";
#    }
#
#    if (!get_reference($build)) {
#        die "Unable to determine reference sequence for build $name\n";
#    }
}

# maybe when running a sort build a line-based index, use perls sort to do a "slow" but memory efficient sort
# have a function that pulls a given field(s) from a line starting at a given fpos


sub get_lines {
    my ($self, %params) = @_;

    #print "[32m".Data::Dumper::Dumper(\%params)."[0m\n";

    my %defaults = (
        offset => 0,
        limit => 500
    );
    
    for (keys %defaults) {
        # fill in defaults in params unless property is already specified
        $params{$_} = $defaults{$_} unless defined $params{$_};
    }

    my $file      = delete $params{file};
    my $offset    = delete $params{offset};
    my $limit     = delete $params{limit};
    my $delimiter = delete $params{delimiter};
    #my $index  = delete $params{index};  # use an index to limit to certain lines
    #my $query  = delete $params{query};  # what to query
    #my $sort   = delete $params{sort};   # what columns to sort by, + for asc, - for desc
    #my $filter = delete $params{filter}; # what columns to filter by

    die $self->error_message("No file or command specified to get lines from") unless defined $file;
    my $fh = IO::File->new("$file") || die $self->error_message("Can't pipe command '$file'.");
    
    my @lines;

    my $count = 0;
    
    # TODO more robust interface
    # TODO test limits and offsets
    while (my $line = <$fh>) {
        if ( ($count >= $offset) && (scalar(@lines) < $limit) ) {
            if (defined ($delimiter)) {
                push @lines, [split($delimiter, $line)];
            } else {
                push @lines, [$line];
            }
        }
        $count++;
    }

    $fh->close();
    
    my %rv = (
        lines_read => $count,
        lines => \@lines
    );
     
    return \%rv;
}

sub build_gene_index {
    my $self = shift;
    my $file = shift;

    my $fh = IO::File->new($file) || die $self->error_message("Can't open $file.");

    my %hash;

    my $cur_gene;
    my $cur_gene_pos;
    
    my $pos = 0;
    
    while (my $line = <$fh>) {
        # get the byte pos of the 7th field without using split:
        my $index = 1+index($line, "\t", 1+index($line, "\t", 1+index(
            $line, "\t", 1+index($line, "\t", 1+index($line, "\t", 1+index($line, "\t"))))));

        # get the byte pos of the end of the 7th field
        my $end = index($line, "\t", $index);

        my $gene = substr($line, $index, $end - $index);
        if (uc $gene ne $cur_gene) { # if genes differ
            if (defined($cur_gene)) { # and there's a previous gene
                if (defined($hash{$cur_gene})) { # this distinciton may be unecessary...
                    push @{$hash{$cur_gene}}, [$cur_gene_pos, $pos];
                } else {
                    $hash{$cur_gene} = [
                        [$cur_gene_pos, $pos]
                    ];
                }
            }

            $cur_gene = uc $gene;
            $cur_gene_pos = $pos;
        }
        $pos = tell($fh);
    }

    # add the last gene
    if (defined($cur_gene)) {
        if (defined($hash{$cur_gene})) { # this distinciton may be unecessary...
            push @{$hash{$cur_gene}}, [$cur_gene_pos, $pos];
        } else {
            $hash{$cur_gene} = [
                [$cur_gene_pos, $pos]
            ];
        }
    }
    
    $fh->close();
    
    return \%hash;
}

sub get_lines_via_gene_index {
    my $self = shift;
    my $file = shift;
    my $bands = shift;
    my $offset = shift;
    my $limit = shift;

    #print "[36m;".Data::Dumper::Dumper($bands)."[0m\n";

    my $fh = IO::File->new($file) || die $self->error_message("Can't open $file.");
    
    my @lines;
    
    my $count = 0;
    if ($limit == -1) { $limit = 500; }
    
    # TODO more robust interface
    # TODO test limits and offsets
    # TODO may want to explicitly sort bands, but they should be in the order of earlier first...
    for my $band (@{$bands}) {
        my $start = $band->[0];
        my $end = $band->[1];
        seek($fh, $start, 0);
        while (tell($fh) != $end) {
            my $line = $fh->getline();
            if ( ($count >= $offset) && (scalar(@lines) < $limit) ) {
                push @lines, [split("\t", $line)];
            }
            $count++;
            die $self->error_message("Seeked too far.") if (tell($fh) > $end);
        }
    }

    $fh->close();
    
    my %rv = (
        lines_read => $count,
        lines => \@lines
    );
     
    return \%rv;
}

sub _generate_content {
    my $self = shift;

    my $build = $self->subject;

    my $annotations_file = join('/', $build->variants_directory, 'filtered.variants.post_annotation');

    my $cache = Genome::Memcache->server();

    if (defined($self->datatables_params())) {
        my $dtparams;
        my @params = split(",",$self->datatables_params());
        while (scalar(@params)) {
            my ($key, $value) = (shift(@params), shift(@params));
            $dtparams->{$key} = $value;
        }
        #print "[35m". Data::Dumper::Dumper($dtparams)."[0m\n";

        my $offset = $dtparams->{'iDisplayStart'};
        my $limit = $dtparams->{'iDisplayLength'};
        my $sEcho = $dtparams->{'sEcho'};
        my $sSearch = $dtparams->{'sSearch'};
        my $sSearchType = $dtparams->{'sSearchType'};
        
        my $key = sprintf("%s-%s-%s", "Genome::Model::Build::View::Annotation::Html", "build-id", $build->id());

        if (!$sSearch) {
            # TODO should return the default view
            my $lines = $self->get_lines(
                file => $annotations_file,
                offset => $offset,
                limit => $limit,
                delimiter => "\t"
            );
            return $self->datatables_response($lines->{lines}, $lines->{lines_read}, $lines->{lines_read}, $sEcho);

        } elsif ($sSearchType eq 'gene') {
            if ($cache && (my $compressed_str = $cache->get($key))) {
                my $value = from_json($compressed_str, {ascii => 1});

                my $index = $value->{'index'};
                
                if (defined($index->{$sSearch})) {
                    my $lines = $self->get_lines_via_gene_index($annotations_file, $index->{$sSearch}, $offset, $limit);
                    #print "[23m".$value->{'line_count'}."[0m\n";
                    return $self->datatables_response($lines->{lines}, $lines->{lines_read}, $value->{'line_count'}, $sEcho);
                } else {
                    return $self->datatables_response([], 0, $value->{'line_count'}, $sEcho);
                }
            } else {
                my $lines = $self->get_lines(
                    file => qq{ awk -F \"\\t\" '\$7 == "$sSearch"' $annotations_file | },
                    offset => $offset,
                    limit => $limit,
                    delimiter => "\t"
                );
                return $self->datatables_response($lines->{lines}, $lines->{lines_read}, 0, $sEcho);
            }
        } elsif ($sSearchType eq 'grep') { # full text grep
            my $escaped_search = $sSearch;
            $escaped_search =~ s/'/'"'"'/g; # disgusting
            $escaped_search = "'$escaped_search'";
            my $lines = $self->get_lines(
                file => "grep $escaped_search $annotations_file | ",
                offset => $offset,
                limit => $limit,
                delimiter => "\t"
            );

            my $lines_count = 0;
            # get lines read from memcache if we can
            if ($cache && (my $compressed_str = $cache->get($key))) {
                my $value = from_json($compressed_str, {ascii => 1});
                $lines_count = $value->{line_count};
            }

            return $self->datatables_response($lines->{lines}, $lines->{lines_read}, $lines_count, $sEcho);
        }
    } elsif (defined($self->request_index())) {
        my $key = sprintf("%s-%s-%s", "Genome::Model::Build::Set::View::Annotation::Html", "build-id", $build->id());

        my $requested = $self->request_index() > 0 ? 1 : 0;
        my $exists = defined($cache->get($key)) ? 1 : 0;
        my $force = $self->request_index() eq '2' ? 1 : 0;
        my $should_build = $cache && $requested && ($force || !$exists) ? 1 : 0;
        
        if ($should_build) {
            #print "[32mBuilding cache with key '$key'[0m\n";

            my @wc = split(/\s+/,`wc -l $annotations_file`);
            my $index = $self->build_gene_index($annotations_file);

            my $data = {
                index => $index,
                line_count => $wc[0],
            };
            
            $cache->set($key, to_json($data, {ascii => 1}), 600);
        }
        
        return to_json({requested => $requested, built => $should_build, existed => $exists, forced => $force}, {ascii => 1});
    } else {
        # this will eventually dumpout test.html with ids and line count filled in and a minimal set of lines at the head
        my @wc = split(/\s+/,`wc -l $annotations_file`);
        return "<p>lines: $wc[0]</p>";
    }
    
}

sub _hashref_to_json_array {
    my ($h) = @_;

    my @arr = map { [$_, $h->{$_}] } sort {$a <=> $b} keys %$h;

    return to_json(\@arr, {ascii => 1});
}

sub datatables_response {
    my ($self, $lines, $filtered_lines, $total_lines, $sEcho) = @_;

    my $rv = {
        "sEcho" => $sEcho,
        "iTotalRecords" => $total_lines,
        "iTotalDisplayRecords" => $filtered_lines,
        "aaData" => $lines
    };

    return to_json($rv, {ascii => 1});
}

1;

