# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/_mzyzyR3wa/australasia.  Olson data version 2014g
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Pacific::Norfolk;
$DateTime::TimeZone::Pacific::Norfolk::VERSION = '1.74';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Pacific::Norfolk::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59958190088, #      utc_end 1900-12-31 12:48:08 (Mon)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59958230400, #    local_end 1901-01-01 00:00:00 (Tue)
40312,
0,
'LMT',
    ],
    [
59958190088, #    utc_start 1900-12-31 12:48:08 (Mon)
61536026880, #      utc_end 1950-12-31 12:48:00 (Sun)
59958230408, #  local_start 1901-01-01 00:00:08 (Tue)
61536067200, #    local_end 1951-01-01 00:00:00 (Mon)
40320,
0,
'NMT',
    ],
    [
61536026880, #    utc_start 1950-12-31 12:48:00 (Sun)
DateTime::TimeZone::INFINITY, #      utc_end
61536068280, #  local_start 1951-01-01 00:18:00 (Mon)
DateTime::TimeZone::INFINITY, #    local_end
41400,
0,
'NFT',
    ],
];

sub olson_version { '2014g' }

sub has_dst_changes { 0 }

sub _max_year { 2024 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}



1;

