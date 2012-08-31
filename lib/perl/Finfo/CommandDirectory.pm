package Finfo::CommandDirectory;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

require File::Basename;
use Data::Dumper;
require Finfo::CommandLineOptions;
use Finfo::ClassUtils 'use_class';
require Term::ANSIColor;
require Text::Wrap;

my %class_base_name :name(class_base_name:r) :isa(string);
my %fnc :name(_functions_and_classes:p) :ds(hashref);

sub START {
    my $self = shift;

    my $class_base_name = $self->class_base_name; 
    my $class_dir = $class_base_name; 
    $class_dir =~ s#::#/#; 
    my %functions_and_classes;
    for my $inc_dir (@INC){ 
        for my $module ( glob( sprintf('%s/%s/*.pm', $inc_dir, $class_dir) ) ) {
            # Class
            my $class = $module;
            $class =~ s#^$inc_dir/##;
            $class =~ s#\.pm$##;
            $class =~ s#/#::#g;

            # Function
            my $function = $class;
            $function =~ s/$class_base_name\:\://;
            my @words = $function =~ /([A-Z](?:[A-Z]*(?=$|[A-Z][a-z])|[a-z]*))/g;
            $function = join('-', map { lc } @words);

            $functions_and_classes{$function} = $class;
        };
    };

    for my $class ( values %functions_and_classes ) {
        use_class($class);
    }

    #print Dumper(\%functions_and_classes);
    $self->_functions_and_classes(\%functions_and_classes);

    return 1;
}

sub execute {
    my $self = shift;

    my $function = shift @ARGV;
    my $functions_and_classes = $self->_functions_and_classes;
    my @functions = keys %$functions_and_classes;
    unless ( $function and $function !~ /^\-/ ) {
        no warnings 'once';
        $Text::Wrap::columns = 100;
        use warnings;
        no warnings 'reserved';

        my $script_name = File::Basename::basename($0);

        print(
            join(
                "\n",
                Term::ANSIColor::colored("Usage", 'bold'),
                "  $script_name <function> [ <options> OR --help ]",
                Term::ANSIColor::colored("Valid functions", 'bold'),
                '',
            )
        );

        for my $function ( sort { $a cmp $b } @functions ) {
            my $help_brief = ( $functions_and_classes->{$function}->can('help_brief') )
            ? $functions_and_classes->{$function}->help_brief
            : 'No description available'; # or nothing at all??
            
            print(
                sprintf(
                    '  %-40s %s',
                    Term::ANSIColor::colored($function, 'yellow'),
                    $help_brief,
                ),
                "\n",
            );
        }

        exit 1;
    }

    Finfo::Validate->validate(
        attr => 'function',
        value => $function,
        isa => [ 'in_list', @functions ],
        msg => 'fatal',
    );

    my $class = $functions_and_classes->{$function};
    my $clo = Finfo::CommandLineOptions->new (
        classes => [ $class ],
        add_q => 1,
    );
    my $opts = $clo->get_options;
    my $command = $class->new( %{ $opts->{$class} } );

    exit not $command->execute;
}

1;

#$HeadURL$
#$Id$
