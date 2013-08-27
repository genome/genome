package Genome::Model::Set::View::Coverage::Xml;

use strict;
use warnings;

class Genome::Model::Set::View::Coverage::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        perspective => { default_value => 'coverage' },
    ],
    has_transient_optional => [
        _genomes_total_bp => {
            is => 'HASH',
            doc => 'store previously computed genome_total_bp to avoid recomputing for each model',
            default_value => {},
        },
        _builds => {
            is_many => 1,
            doc => 'store computed list of builds for models',
        },
    ],
};

sub all_subject_classes {
    my $self = shift;
    my @classes = $self->SUPER::all_subject_classes;

    #If more sophisticated handling is required,
    #can substitute the particular classes of model
    #returned by $self->members.  This is quick and
    #sufficient for now.
    unless(grep($_->isa('Genome::Model'), @classes)) {
        push @classes, 'Genome::Model';
    }

    return @classes;
}

sub members {
    my $self = shift;

    my $set = $self->subject;
    my @members = $set->members;

    return @members;
}

sub _generate_content {
    my $self = shift;
    my $subject = $self->subject();

    my $xml_doc = XML::LibXML->createDocument();
    $self->_xml_doc($xml_doc);

    my $object = $xml_doc->createElement('object');
    $xml_doc->setDocumentElement($object);
    my $time = UR::Context->current->now();
    $object->addChild( $xml_doc->createAttribute('id',$subject->id) );
    $object->addChild( $xml_doc->createAttribute('generated-at',$time) );

    my $name;
    if($subject->can('name') and not $subject->isa('UR::Object::Set')) {
        $name = $subject->name;
    } elsif($subject->can('rule_display')) {
        $name = $subject->rule_display;
        $name =~ s/^UR::BoolExpr=\([\w:]+/Set: /;
        $name =~ s/_/-/g;
        $name =~ s/ =>/:/g;
        $name =~ s/"//g;
        $name =~ s/\)$//;
        $name =~ s/([\w\d]),([\w\d])/$1, $2/g;
    } else {
        $name = $subject->__display_name__;
    }

    $object->addChild( $xml_doc->createAttribute('display_name',$name) );

    my $type = $subject->class;
    $type = 'Genome::Model::Build' if $type->isa('Genome::Model::Build');
    $object->addChild( $xml_doc->createAttribute('type', $type));

    $object->addChild( $self->get_enrichment_factor_node() );
    $object->addChild( $self->get_enrichment_factor_v2_node() );
    $object->addChild( $self->get_alignment_summary_node() );
    $object->addChild( $self->get_alignment_summary_v2_node() );
    $object->addChild( $self->get_coverage_summary_node() );

    return $xml_doc->toString(1);
}

sub builds {
    my $self = shift;

    my @b = $self->_builds;
    return @b if @b;

    #preload data for efficiency
    my @members = $self->members;
    my @model_ids = map($_->id, @members);
    my @builds = Genome::Model::Build->get(model_id => \@model_ids);

    @b = map($self->get_last_succeeded_coverage_stats_build_from_model($_), $self->members);

    $self->_builds(\@b);
    return @b;
}

sub get_last_succeeded_coverage_stats_build_from_model {
    my $self = shift;
    my $model = shift;
    # results from $model->builds are now automatically sorted by date_scheduled.
    my @sorted_builds = $model->builds;
    for my $build (@sorted_builds) {
        my @r = $build->result_users('software_result.subclass_name' => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
        if(scalar @r) {
            return $build;
        }

        #old way of storing coverage report
        if($build->isa('Genome::Model::Build::ReferenceAlignment')){
            my @events = $build->the_events;
            my @coverage_stats_events = grep { $_->class eq 'Genome::Model::Event::Build::ReferenceAlignment::CoverageStats' } @events;
            if (scalar(@coverage_stats_events) == 1) {
                my $coverage_stats_event = $coverage_stats_events[0];
                if ($coverage_stats_event->event_status eq 'Succeeded') {
                    return $build;
                }
            }
        }
    }
    return;
}

sub get_enrichment_factor_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $ef_node = $xml_doc->createElement('enrichment-factor');
    for my $build ($self->builds) {
        my $model = $build->model;

        my @results = $self->_results_for_build($build);

        for my $result (@results) {
            my $model_node = $self->_get_model_node($model, $result);
            $ef_node->addChild( $model_node );

            # get BED file
            my $bedf;
            my $refcovd;

            if($result->isa('Genome::SoftwareResult')) {
                $refcovd = $result->output_dir;
            } else {
                $refcovd = $result->reference_coverage_directory;
            }
            opendir(my $refcovdh, $refcovd) or die "Could not open reference coverage directory at $refcovd";

            while (my $file = readdir($refcovdh)) {
                if ($file =~ /.*.bed/) { $bedf = $refcovd . "/" . $file; }
            }

            # calculate target_total_bp
            my $target_total_bp;

            open(my $bedfh, "<", $bedf) or die "Could not open BED file at $bedf";

            while (<$bedfh>) {
                chomp;
                my @f      = split (/\t/, $_);
                my $start  = $f[1];
                my $stop   = $f[2];
                my $length = ($stop - $start);
                $target_total_bp += $length;
            }

            # calculate genome_total_bp from reference sequence seqdict.sam
            my $genome_total_bp;

            my $refseq_build;
            if($result->isa('Genome::SoftwareResult')) {
                $refseq_build = $result->alignment_result->reference_build;
            } else {
                $refseq_build = $result->reference_sequence_build;
            }

            my $seqdictf = $refseq_build->data_directory . "/seqdict/seqdict.sam";

            my $genomes_total_bp = $self->_genomes_total_bp;
            if($genomes_total_bp and $genomes_total_bp->{$seqdictf}) {
                $genome_total_bp = $genomes_total_bp->{$seqdictf};
            } else {
                open(my $seqdictfh, "<", $seqdictf) or die "Could not open seqdict at $seqdictf";

                while (<$seqdictfh>) {
                    chomp;
                    unless($_ =~ /$@HD/) { # skip the header row
                        my @f = split(/\t/, $_);
                        my $ln = $f[2];
                        $ln =~ s/LN://;
                        $genome_total_bp += $ln;
                    }
                }

                $genomes_total_bp->{$seqdictf} = $genome_total_bp;
                $self->_genomes_total_bp($genomes_total_bp);
            }

            # get wingspan 0 alignment metrics
            # Note: this may need to be changed to alignment_summary_v2_hash_ref()
            my $ws_zero = $result->alignment_summary_hash_ref->{'0'};

            # calculate enrichment factor!
            my $myEF = Genome::Model::Tools::TechD::CaptureEnrichmentFactor->execute(
                capture_unique_bp_on_target    => $ws_zero->{'unique_target_aligned_bp'},
                capture_duplicate_bp_on_target => $ws_zero->{'duplicate_target_aligned_bp'},
                capture_total_bp               => $ws_zero->{'total_aligned_bp'} + $ws_zero->{'total_unaligned_bp'},
                target_total_bp                => $target_total_bp,
                genome_total_bp                => $genome_total_bp
            );

            my $theoretical_max_enrichment_factor = 0;
            my $unique_on_target_enrichment_factor = 0;
            my $total_on_target_enrichment_factor = 0;

            if ($myEF) {
                $theoretical_max_enrichment_factor  = $myEF->theoretical_max_enrichment_factor();
                $unique_on_target_enrichment_factor = $myEF->unique_on_target_enrichment_factor();
                $total_on_target_enrichment_factor  = $myEF->total_on_target_enrichment_factor();
            }

            my $uotef_node = $model_node->addChild( $xml_doc->createElement('unique_on_target_enrichment_factor') );
            $uotef_node->addChild( $xml_doc->createTextNode( $unique_on_target_enrichment_factor ) );

            my $totef_node = $model_node->addChild( $xml_doc->createElement('total_on_target_enrichment_factor') );
            $totef_node->addChild( $xml_doc->createTextNode( $total_on_target_enrichment_factor ) );

            my $tmef_node = $model_node->addChild( $xml_doc->createElement('theoretical_max_enrichment_factor') );
            $tmef_node->addChild( $xml_doc->createTextNode( $theoretical_max_enrichment_factor ) );
        }
    }

    return $ef_node;
}

sub get_enrichment_factor_v2_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $ef_node = $xml_doc->createElement('enrichment-factor-v2');
    for my $build ($self->builds) {
        my $model = $build->model;

        my @results = $self->_results_for_build($build);

        for my $result (@results) {
            my $model_node = $self->_get_model_node($model, $result);

            # get BED file
            my $bedf;
            my $refcovd;

            if($result->isa('Genome::SoftwareResult')) {
                $refcovd = $result->output_dir;
            } else {
                $refcovd = $result->reference_coverage_directory;
            }
            opendir(my $refcovdh, $refcovd) or die "Could not open reference coverage directory at $refcovd";

            while (my $file = readdir($refcovdh)) {
                if ($file =~ /.*.bed/) { $bedf = $refcovd . "/" . $file; }
            }

            # calculate target_total_bp
            my $target_total_bp;

            open(my $bedfh, "<", $bedf) or die "Could not open BED file at $bedf";

            while (<$bedfh>) {
                chomp;
                my @f      = split (/\t/, $_);
                my $start  = $f[1];
                my $stop   = $f[2];
                my $length = ($stop - $start);
                $target_total_bp += $length;
            }

            # calculate genome_total_bp from reference sequence seqdict.sam
            my $genome_total_bp;

            my $refseq_build;
            if($result->isa('Genome::SoftwareResult')) {
                $refseq_build = $result->alignment_result->reference_build;
            } else {
                $refseq_build = $result->reference_sequence_build;
            }

            my $seqdictf = $refseq_build->data_directory . "/seqdict/seqdict.sam";

            my $genomes_total_bp = $self->_genomes_total_bp;
            if($genomes_total_bp and $genomes_total_bp->{$seqdictf}) {
                $genome_total_bp = $genomes_total_bp->{$seqdictf};
            } else {
                open(my $seqdictfh, "<", $seqdictf) or die "Could not open seqdict at $seqdictf";

                while (<$seqdictfh>) {
                    chomp;
                    unless($_ =~ /$@HD/) { # skip the header row
                        my @f = split(/\t/, $_);
                        my $ln = $f[2];
                        $ln =~ s/LN://;
                        $genome_total_bp += $ln;
                    }
                }

                $genomes_total_bp->{$seqdictf} = $genome_total_bp;
                $self->_genomes_total_bp($genomes_total_bp);
            }

            # get wingspan 0 alignment metrics
            # Note: this may need to be changed to alignment_summary_v2_hash_ref()
            my $ws_zero = eval { $result->alignment_summary_v2_hash_ref->{'0'} };
            if ($@) {
                next;
            }

            # calculate enrichment factor!
            my $myEF = Genome::Model::Tools::TechD::CaptureEnrichmentFactor->execute(
                capture_unique_bp_on_target    => $ws_zero->{'unique_target_aligned_bp'},
                capture_duplicate_bp_on_target => $ws_zero->{'duplicate_target_aligned_bp'},
                capture_total_bp               => $ws_zero->{'total_aligned_bp'} + $ws_zero->{'total_unaligned_bp'},
                target_total_bp                => $target_total_bp,
                genome_total_bp                => $genome_total_bp
            );

            my $theoretical_max_enrichment_factor = 0;
            my $unique_on_target_enrichment_factor = 0;
            my $total_on_target_enrichment_factor = 0;

            if ($myEF) {
                $theoretical_max_enrichment_factor  = $myEF->theoretical_max_enrichment_factor();
                $unique_on_target_enrichment_factor = $myEF->unique_on_target_enrichment_factor();
                $total_on_target_enrichment_factor  = $myEF->total_on_target_enrichment_factor();
            }

            $ef_node->addChild( $model_node );

            my $uotef_node = $model_node->addChild( $xml_doc->createElement('unique_on_target_enrichment_factor') );
            $uotef_node->addChild( $xml_doc->createTextNode( $unique_on_target_enrichment_factor ) );

            my $totef_node = $model_node->addChild( $xml_doc->createElement('total_on_target_enrichment_factor') );
            $totef_node->addChild( $xml_doc->createTextNode( $total_on_target_enrichment_factor ) );

            my $tmef_node = $model_node->addChild( $xml_doc->createElement('theoretical_max_enrichment_factor') );
            $tmef_node->addChild( $xml_doc->createTextNode( $theoretical_max_enrichment_factor ) );
        }
    }

    return $ef_node;
}


sub get_alignment_summary_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $as_node = $xml_doc->createElement('alignment-summary');

    for my $build ($self->builds) {
        my $model = $build->model;

        my @results = $self->_results_for_build($build);

        for my $result (@results) {
            my $model_node = $self->_get_model_node($model, $result);
            $as_node->addChild( $model_node );

            my $alignment_summary_hash_ref = $result->alignment_summary_hash_ref;
            for my $ws_key (keys %{$alignment_summary_hash_ref}) {
                my $ws_node = $model_node->addChild( $xml_doc->createElement('wingspan') );
                $ws_node->addChild( $xml_doc->createAttribute('size', $ws_key) );
                for my $param_key (keys %{$alignment_summary_hash_ref->{$ws_key}}) {
                    my $key_node = $ws_node->addChild( $xml_doc->createElement($param_key) );
                    $key_node->addChild( $xml_doc->createTextNode( $alignment_summary_hash_ref->{$ws_key}->{$param_key} ) );
                }
            }
        }
    }
    return $as_node;
}

sub get_alignment_summary_v2_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $as_node = $xml_doc->createElement('alignment-summary-v2');

    for my $build ($self->builds) {
        my $model = $build->model;

        my @results = $self->_results_for_build($build);

        for my $result (@results) {
            my $model_node = $self->_get_model_node($model, $result);

            my $alignment_summary_hash_ref = eval { $result->alignment_summary_v2_hash_ref };
            if ($@) {
                next;
            }

            $as_node->addChild( $model_node );
            for my $ws_key (keys %{$alignment_summary_hash_ref}) {
                my $ws_node = $model_node->addChild( $xml_doc->createElement('wingspan') );
                $ws_node->addChild( $xml_doc->createAttribute('size', $ws_key) );
                for my $param_key (keys %{$alignment_summary_hash_ref->{$ws_key}}) {
                    my $key_node = $ws_node->addChild( $xml_doc->createElement($param_key) );
                    $key_node->addChild( $xml_doc->createTextNode( $alignment_summary_hash_ref->{$ws_key}->{$param_key} ) );
                }
            }
        }
    }
    return $as_node;
}


sub get_result_depths {

    my ($result) = @_;
    my @result_depths;

    if($result->isa('Genome::SoftwareResult')) {
        @result_depths = split(',', $result->minimum_depths);
    } else {
        @result_depths = @{ $result->minimum_depths_array_ref };
    }

    return @result_depths;
}

sub get_coverage_summary_node {
    my $self = shift;
    my $xml_doc = $self->_xml_doc;

    my $cs_node = $xml_doc->createElement('coverage-summary');
    my @builds  = $self->builds();
    my $model_ids_by_depths = {};

    my @min_depths;
    my %builds_with_depth_mismatch;
    BUILD:
    for my $build (@builds) {
        my $model = $build->model;

        my @results = $self->_results_for_build($build);

        for my $result (@results) {
            my @result_depths = get_result_depths($result);
            $model_ids_by_depths->{join(',',@result_depths)}->{$model->id}++;

            unless (@min_depths) {
                @min_depths = sort{ $a <=> $b } @result_depths;
                for my $min_depth (@min_depths) {
                    my $header_node = $cs_node->addChild( $xml_doc->createElement('minimum_depth_header') );
                    $header_node->addChild( $xml_doc->createAttribute('value',$min_depth) );
                }
            } else {
                my @other_min_depths = sort{ $a <=> $b } @result_depths;
                unless (scalar(@min_depths) == scalar(@other_min_depths)) {
                    warn('Model '. $model->name .' has '. scalar(@other_min_depths) .' minimum_depth filters expecting '. scalar(@min_depths) .' minimum_depth filters');
                    next BUILD;
                }
                for (my $i = 0; $i < scalar(@min_depths); $i++) {
                    my $expected_min_depth = $min_depths[$i];
                    my $other_min_depth = $other_min_depths[$i];
                    unless ($expected_min_depth == $other_min_depth) {
                        warn('Model '. $model->name .' has '. $other_min_depth .' minimum_depth filter expecting '. $expected_min_depth .' minimum_depth filter');
                        next BUILD;
                    }
                }
            }

            my $model_node = $self->_get_model_node($model, $result);
            $cs_node->addChild( $model_node );

            my $coverage_stats_summary_hash_ref = $result->coverage_stats_summary_hash_ref;
            for my $min_depth (keys %{$coverage_stats_summary_hash_ref->{0}}) {
                unless (grep {$_ eq $min_depth} @min_depths) {
                    $builds_with_depth_mismatch{$result->build_id} = $min_depth;
                    next BUILD;
                }
                
                my $min_depth_node = $model_node->addChild( $xml_doc->createElement('minimum_depth') );
                $min_depth_node->addChild( $xml_doc->createAttribute('value',$min_depth) );
                for my $key (keys %{$coverage_stats_summary_hash_ref->{0}->{$min_depth}}) {
                    my $key_node = $min_depth_node->addChild( $xml_doc->createElement($key) );
                    $key_node->addChild( $xml_doc->createTextNode( $coverage_stats_summary_hash_ref->{0}->{$min_depth}->{$key} ) );
                }
            }
        }
    }

    # things have gone all wrong so return an error node with links to separate graphs
    # grouped by depth grouping
    if (keys %$model_ids_by_depths > 1) {
        $cs_node = $self->get_error_node('Found conflicting depth groupings (different processing profiles?)');

        for my $depth_grouping (keys %$model_ids_by_depths) {
            my $group_node = $cs_node->addChild( $xml_doc->createElement('depth_group') );
            $group_node->addChild( $xml_doc->createAttribute('description',$depth_grouping) );
            my $url = '/view/genome/model/set/coverage.html?id='
                    .  join('&id=', keys %{$model_ids_by_depths->{$depth_grouping}} );
            $group_node->addChild( $xml_doc->createAttribute('url',$url) );
        }
    }

    # different things have gone wrong: some builds contain data for depths not expected by the processing profile
    if (keys %builds_with_depth_mismatch > 0) {
        my $error = "ERROR: Could not render coverage charts and tables, as some builds in this set provided spurrious data for depths not indicated by their processing profile (expecting depths " . join(', ', @min_depths) . ").";
        $cs_node = $self->get_error_node($error);

        my $builds_node = $cs_node->addChild( $xml_doc->createElement('erroneous_builds') );
        $builds_node->addChild( $xml_doc->createAttribute('expected_depths', join(',', @min_depths)) );
        for my $build_id (keys %builds_with_depth_mismatch) {
            my $build_node = $builds_node->addChild( $xml_doc->createElement('build') );
            $build_node->addChild( $xml_doc->createAttribute('id',$build_id) );
            $build_node->addChild( $xml_doc->createAttribute('unexpected_depth', $builds_with_depth_mismatch{$build_id}) );
        }

    }

    return $cs_node;
}

sub get_error_node {

    my ($self, $msg) = @_;
    my $xml_doc = $self->_xml_doc;
    my $node = $xml_doc->createElement('coverage-summary');
    $node->addChild($xml_doc->createAttribute('error',$msg));

#        my $header_node = $cs_node->addChild( $xml_doc->createElement('error') );
#        $header_node->addChild( $xml_doc->createAttribute('value',$msg) );

    return $node;
}

sub _results_for_build {
    my $self = shift;
    my $build = shift;

    my @result_users = $build->result_users('software_result.subclass_name' => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
    my @results = map($_->software_result, @result_users);

    unless(@results) {
        push @results, $build;
    }

    return @results;
}

sub _get_model_node {
    my $self = shift;
    my $model = shift;
    my $result = shift;

    my $xml_doc = $self->_xml_doc;

    my @idata;
    if($result->isa('Genome::SoftwareResult')) {
        @idata = $result->alignment_result->instrument_data;
    } else { #old build way
        @idata = $result->instrument_data;
    }

    my $subject = $idata[0]->sample;

    my $model_node = $xml_doc->createElement('model');
    $model_node->addChild( $xml_doc->createAttribute('id',$model->id));
    $model_node->addChild( $xml_doc->createAttribute('subject_name',$subject->name));
    $model_node->addChild( $xml_doc->createAttribute('model_name',$model->name));

    $model_node->addChild( $xml_doc->createAttribute('lane_count', scalar(@idata)) );

    if($result->isa('Genome::SoftwareResult')) {
        $model_node->addChild( $xml_doc->createAttribute('result_id', $result->id) );
    }

    return $model_node;
}


1;
