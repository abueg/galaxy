# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/fQcgigXnA8/northamerica.  Olson data version 2015e
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::Puerto_Rico;
$DateTime::TimeZone::America::Puerto_Rico::VERSION = '1.92';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::Puerto_Rico::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59902647865, #      utc_end 1899-03-28 16:24:25 (Tue)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59902632000, #    local_end 1899-03-28 12:00:00 (Tue)
-15865,
0,
'LMT',
    ],
    [
59902647865, #    utc_start 1899-03-28 16:24:25 (Tue)
61262625600, #      utc_end 1942-05-03 04:00:00 (Sun)
59902633465, #  local_start 1899-03-28 12:24:25 (Tue)
61262611200, #    local_end 1942-05-03 00:00:00 (Sun)
-14400,
0,
'AST',
    ],
    [
61262625600, #    utc_start 1942-05-03 04:00:00 (Sun)
61366287600, #      utc_end 1945-08-14 23:00:00 (Tue)
61262614800, #  local_start 1942-05-03 01:00:00 (Sun)
61366276800, #    local_end 1945-08-14 20:00:00 (Tue)
-10800,
1,
'AWT',
    ],
    [
61366287600, #    utc_start 1945-08-14 23:00:00 (Tue)
61370283600, #      utc_end 1945-09-30 05:00:00 (Sun)
61366276800, #  local_start 1945-08-14 20:00:00 (Tue)
61370272800, #    local_end 1945-09-30 02:00:00 (Sun)
-10800,
1,
'APT',
    ],
    [
61370283600, #    utc_start 1945-09-30 05:00:00 (Sun)
61378315200, #      utc_end 1946-01-01 04:00:00 (Tue)
61370269200, #  local_start 1945-09-30 01:00:00 (Sun)
61378300800, #    local_end 1946-01-01 00:00:00 (Tue)
-14400,
0,
'AST',
    ],
    [
61378315200, #    utc_start 1946-01-01 04:00:00 (Tue)
DateTime::TimeZone::INFINITY, #      utc_end
61378300800, #  local_start 1946-01-01 00:00:00 (Tue)
DateTime::TimeZone::INFINITY, #    local_end
-14400,
0,
'AST',
    ],
];

sub olson_version {'2015e'}

sub has_dst_changes {2}

sub _max_year {2025}

sub _new_instance {
    return shift->_init( @_, spans => $spans );
}



1;

