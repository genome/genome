package Genome::Subject::Command::Import;

use strict;
use warnings;

use Genome;
use MIME::Types;
use MIME::Base64;
use Text::CSV;
use IO::Scalar;

class Genome::Subject::Command::Import {
    is => 'Command',
    has => [
       nomenclature_id => { is => 'Text' },
       nomenclature    => { is => 'Genome::Nomenclature', id_by=>'nomenclature_id' },
       nomenclature_name    => { is => 'Text', via=>'nomenclature', to=>'name' },
       subclass_name   => { is => 'Text' },
       project_name    => { is => 'Text' },
       decoded_content => { calculate_from => ['content'], calculate => q|MIME::Base64::decode_base64($content)|}
    ],
    has_optional => [
       content  => { is => 'Text' },
       filename => { is => 'Text' }
    ],
};

sub help_brief {
    return 'Import subjects/attributes from web';
}

sub create {
    my ($class, %p) = @_;

    if (my $f = $p{'filename'}) {

        if ($p{'content'}) {
            warn "You passed a filename and some content- ignoring the filename, using content";
        } else {
            $p{'content'} = do {
                local $/ = undef;
                open(my $fh, $f);
                my $c = <$fh>;
                close($fh); 
                $c;
            };
        }
    } elsif (! $p{'content'}) {
        die "Error: you must pass either --filename or --content as an argument";
    }
 
    return $class->SUPER::create(%p);
}

sub execute {
    my ($self) = @_;

    # one project per upload
    my $project = Genome::Project->get_or_create(name => $self->project_name);

    # Assumes first row contains column names
    # Assumes first col is the name of the object or
    #   blank to create a new object

    my $subclass_name = $self->subclass_name();

    # data can come in from a filename or directly (and encoded)
    my $raw = $self->filename ? $self->content() : $self->decoded_content();

    my $fh = new IO::Scalar \$raw;
    my $csv = Text::CSV->new();
    my @header;
    my $field = {};

    my $i = 0;
    my $changed;
    my $seen = {};
    my $value_count_doesnt_match;

    ROW:
    while (my $row = $csv->getline($fh)) {

        if ( $i++ == 0 ) { 
            @header = @$row; 
            $field = $self->check_types(@header);
            next ROW; 
        }

        my $name = $row->[0];
        if ($seen->{$name}) {
            die "Error: Found non-unique row ($name) in the spreadsheet- should be one sample or individual per row";
        }
        $seen->{$name}++;
        my @values = @$row;  

        if (@header != @values) {
            $value_count_doesnt_match++;
        }

        my (@obj) = $subclass_name->get(name => $name);

        if ( !@obj ) {
            warn "Skipping row- couldnt get a $subclass_name object with name: " . $row->[0];
            next ROW;
        }


        my $j = -1;
        VALUE:
        for my $v (@values) {

            if (++$j == 0) { next VALUE; } # first value should be the obj's "name"

            next VALUE if !$v;

            my $col_name = $header[$j];
            my $f = $field->{$col_name};

            if (!$f) {
                next VALUE; # skipping because the column wasnt named
            }

            if (! $self->validate($f, $v)) {
                warn "Skipping value '$v' because its not a valid " . $f->type();
                next VALUE;
            }

            if (@obj > 1) {
                warn "** Found multiple things with name $name so applying attributes to all of them";
            }

            OBJECT:
            for my $o (@obj) {

                my $sa = Genome::SubjectAttribute->get(
                    subject_id      => $o->id,
                    attribute_label => $col_name,
                    nomenclature    => $f->id
                );

                if ($sa) { 
                    $sa->delete;
                }
                $sa = Genome::SubjectAttribute->create(
                    subject_id      => $o->id,
                    attribute_label => $col_name,
                    attribute_value => $v,
                    nomenclature    => $f->id
                );
                $changed++; 
            }
        }

        # TODO: check for unique: subclass_name, id
        # add each subject obj to the project
        for my $o (@obj) {
            if (! $project->get_part($o)) {
                $project->add_part( entity => $o, role => 'automatic');
            }
        }
    }

    if ($value_count_doesnt_match) {
        warn "** Used $value_count_doesnt_match rows that have different number of values than the header (first row) indicates there should be.";
    }

    return $changed;
}


sub validate {
    my ($self, $field, $value) = @_;

    my $nom = $self->nomenclature();
    my $expected_type = $field->type();

    if ($nom->empty_equivalent &&
        ($value eq $nom->empty_equivalent) ) {
        return 1; 
    }

    if ($expected_type eq 'string') {
        return 1;
    } 

    if ($expected_type eq 'integer'
        or $expected_type eq 'numeric') {
        if ($value =~ /^[+-]?\d+$/) {
            return 1;
        }
    }  

    if ($expected_type eq 'real'
        or $expected_type eq 'numeric') {
        if ($value =~ /^[+-]?\d*\.?\d*$/) {
            return 1;
        }
    }

    if ($expected_type eq 'date') {
        if ($value =~ /^\d{4}-\d{2}-\d{2}$/) {
            return 1;
        }
    }

    if ($expected_type eq 'enumerated') {
        my @acceptable_values = map {$_->value} $field->enumerated_values();
        if (grep /^$value$/, @acceptable_values) {
            return 1;
        }
    }

    return;
}

sub check_types {

    my ($self, @header) = @_;

    my $nom = $self->nomenclature();

    my @fields = $nom->fields();
    my %field = map {$_->name => $_} @fields;

    my $i = 0;
    COLUMN_NAME:
    for my $h (@header) {

        if ($i++ == 0) { next COLUMN_NAME; } # first is the subject name

        if (!$h) {
            warn "! NO COLUMN NAME for field number $i- going to skip all values for that column";
            next COLUMN_NAME;
        }

        if (!defined($field{$h})) {

            if ($nom->accepts_any_field) {
                my $new_field = Genome::Nomenclature::Field->create(
                    name => $h,
                    type => 'string',
                    nomenclature_id => $nom->id,
                );

                $field{$h} = $new_field;
            } else {
                die "Error: column '$h' is not a field in nomenclature '"
                    . $nom->name . "'";
            }
        }
    }

    return \%field;
}


1;


