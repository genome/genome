#!/usr/bin/env perl

# This test is to try prevent developers from mistakenly forgetting or assuming
# that properties on SoftwareSubclasses will be saved to the database. If a
# subclass does not specifiy an alternate table name then UR assumes the
# properties do not need to be saved to the superclasses data source. So this
# test is to try to force properties to be set to either is_input, is_param, or
# is_transient (for non-calculated or non-delegated properties).
#
# We had ~300k objects affected by this, see RT #88011.

use strict;
use warnings;
use Test::More;

use above qw(Genome);

ok(-f '../Genome.pm', 'currently at top of namespace') or die;

my $sr_classname = 'Genome::SoftwareResult';
my @classnames = qx(UR_DBI_NO_COMMIT=0 ur show subclasses --superclass $sr_classname --flat --recalculate);
is($?, 0, 'ur show subclasses command exited 0') or die;
chomp @classnames;
@classnames = grep { defined $_ && $_ && $_ =~ /^[\w:]+$/ } @classnames;
ok(@classnames > 0, "got subclasses of $sr_classname");

my @bad_classnames;
for my $classname (@classnames) {
    my $class = $classname->__meta__;
    ok($class, "loaded $sr_classname");
    ok($classname->isa($sr_classname), "confirmed subclass of $sr_classname");

    my @ambiguous_properties;
    # classes that have table_name will complain if no column exists in DB for a property
    unless ($class->table_name) {
        my @properties = $class->properties;
        # try to define what "ambiguous" means...
        @ambiguous_properties = grep { !(
                $_->property_name eq 'subclass_name'
                || $_->is_param
                || $_->is_input
                || $_->is_calculated
                || $_->is_delegated
                || $_->class_name ne $classname
                || $_->is_transient
            )} @properties;
    }

    ok(@ambiguous_properties == 0, "$classname has no ambiguous properties") or diag sprintf("ambiguous properties: %s\n", join(", ", map { $_->property_name } @ambiguous_properties));
}

done_testing;
