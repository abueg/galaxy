# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/sdoP4iLTHM/northamerica.  Olson data version 2015a
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::Santa_Isabel;
$DateTime::TimeZone::America::Santa_Isabel::VERSION = '1.85';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::Santa_Isabel::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
60620947200, #      utc_end 1922-01-01 08:00:00 (Sun)
DateTime::TimeZone::NEG_INFINITY, #  local_start
60620919632, #    local_end 1922-01-01 00:20:32 (Sun)
-27568,
0,
'LMT',
    ],
    [
60620947200, #    utc_start 1922-01-01 08:00:00 (Sun)
60684015600, #      utc_end 1924-01-01 07:00:00 (Tue)
60620922000, #  local_start 1922-01-01 01:00:00 (Sun)
60683990400, #    local_end 1924-01-01 00:00:00 (Tue)
-25200,
0,
'MST',
    ],
    [
60684015600, #    utc_start 1924-01-01 07:00:00 (Tue)
60792620400, #      utc_end 1927-06-11 07:00:00 (Sat)
60683986800, #  local_start 1923-12-31 23:00:00 (Mon)
60792591600, #    local_end 1927-06-10 23:00:00 (Fri)
-28800,
0,
'PST',
    ],
    [
60792620400, #    utc_start 1927-06-11 07:00:00 (Sat)
60900879600, #      utc_end 1930-11-15 07:00:00 (Sat)
60792595200, #  local_start 1927-06-11 00:00:00 (Sat)
60900854400, #    local_end 1930-11-15 00:00:00 (Sat)
-25200,
0,
'MST',
    ],
    [
60900879600, #    utc_start 1930-11-15 07:00:00 (Sat)
60912720000, #      utc_end 1931-04-01 08:00:00 (Wed)
60900850800, #  local_start 1930-11-14 23:00:00 (Fri)
60912691200, #    local_end 1931-04-01 00:00:00 (Wed)
-28800,
0,
'PST',
    ],
    [
60912720000, #    utc_start 1931-04-01 08:00:00 (Wed)
60928441200, #      utc_end 1931-09-30 07:00:00 (Wed)
60912694800, #  local_start 1931-04-01 01:00:00 (Wed)
60928416000, #    local_end 1931-09-30 00:00:00 (Wed)
-25200,
1,
'PDT',
    ],
    [
60928441200, #    utc_start 1931-09-30 07:00:00 (Wed)
61261862400, #      utc_end 1942-04-24 08:00:00 (Fri)
60928412400, #  local_start 1931-09-29 23:00:00 (Tue)
61261833600, #    local_end 1942-04-24 00:00:00 (Fri)
-28800,
0,
'PST',
    ],
    [
61261862400, #    utc_start 1942-04-24 08:00:00 (Fri)
61366287600, #      utc_end 1945-08-14 23:00:00 (Tue)
61261837200, #  local_start 1942-04-24 01:00:00 (Fri)
61366262400, #    local_end 1945-08-14 16:00:00 (Tue)
-25200,
1,
'PWT',
    ],
    [
61366287600, #    utc_start 1945-08-14 23:00:00 (Tue)
61374006000, #      utc_end 1945-11-12 07:00:00 (Mon)
61366262400, #  local_start 1945-08-14 16:00:00 (Tue)
61373980800, #    local_end 1945-11-12 00:00:00 (Mon)
-25200,
1,
'PPT',
    ],
    [
61374006000, #    utc_start 1945-11-12 07:00:00 (Mon)
61449609600, #      utc_end 1948-04-05 08:00:00 (Mon)
61373977200, #  local_start 1945-11-11 23:00:00 (Sun)
61449580800, #    local_end 1948-04-05 00:00:00 (Mon)
-28800,
0,
'PST',
    ],
    [
61449609600, #    utc_start 1948-04-05 08:00:00 (Mon)
61474143600, #      utc_end 1949-01-14 07:00:00 (Fri)
61449584400, #  local_start 1948-04-05 01:00:00 (Mon)
61474118400, #    local_end 1949-01-14 00:00:00 (Fri)
-25200,
1,
'PDT',
    ],
    [
61474143600, #    utc_start 1949-01-14 07:00:00 (Fri)
61630790400, #      utc_end 1954-01-01 08:00:00 (Fri)
61474114800, #  local_start 1949-01-13 23:00:00 (Thu)
61630761600, #    local_end 1954-01-01 00:00:00 (Fri)
-28800,
0,
'PST',
    ],
    [
61630790400, #    utc_start 1954-01-01 08:00:00 (Fri)
61640647200, #      utc_end 1954-04-25 10:00:00 (Sun)
61630761600, #  local_start 1954-01-01 00:00:00 (Fri)
61640618400, #    local_end 1954-04-25 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61640647200, #    utc_start 1954-04-25 10:00:00 (Sun)
61653949200, #      utc_end 1954-09-26 09:00:00 (Sun)
61640622000, #  local_start 1954-04-25 03:00:00 (Sun)
61653924000, #    local_end 1954-09-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61653949200, #    utc_start 1954-09-26 09:00:00 (Sun)
61672096800, #      utc_end 1955-04-24 10:00:00 (Sun)
61653920400, #  local_start 1954-09-26 01:00:00 (Sun)
61672068000, #    local_end 1955-04-24 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61672096800, #    utc_start 1955-04-24 10:00:00 (Sun)
61685398800, #      utc_end 1955-09-25 09:00:00 (Sun)
61672071600, #  local_start 1955-04-24 03:00:00 (Sun)
61685373600, #    local_end 1955-09-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61685398800, #    utc_start 1955-09-25 09:00:00 (Sun)
61704151200, #      utc_end 1956-04-29 10:00:00 (Sun)
61685370000, #  local_start 1955-09-25 01:00:00 (Sun)
61704122400, #    local_end 1956-04-29 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61704151200, #    utc_start 1956-04-29 10:00:00 (Sun)
61717453200, #      utc_end 1956-09-30 09:00:00 (Sun)
61704126000, #  local_start 1956-04-29 03:00:00 (Sun)
61717428000, #    local_end 1956-09-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61717453200, #    utc_start 1956-09-30 09:00:00 (Sun)
61735600800, #      utc_end 1957-04-28 10:00:00 (Sun)
61717424400, #  local_start 1956-09-30 01:00:00 (Sun)
61735572000, #    local_end 1957-04-28 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61735600800, #    utc_start 1957-04-28 10:00:00 (Sun)
61748902800, #      utc_end 1957-09-29 09:00:00 (Sun)
61735575600, #  local_start 1957-04-28 03:00:00 (Sun)
61748877600, #    local_end 1957-09-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61748902800, #    utc_start 1957-09-29 09:00:00 (Sun)
61767050400, #      utc_end 1958-04-27 10:00:00 (Sun)
61748874000, #  local_start 1957-09-29 01:00:00 (Sun)
61767021600, #    local_end 1958-04-27 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61767050400, #    utc_start 1958-04-27 10:00:00 (Sun)
61780352400, #      utc_end 1958-09-28 09:00:00 (Sun)
61767025200, #  local_start 1958-04-27 03:00:00 (Sun)
61780327200, #    local_end 1958-09-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61780352400, #    utc_start 1958-09-28 09:00:00 (Sun)
61798500000, #      utc_end 1959-04-26 10:00:00 (Sun)
61780323600, #  local_start 1958-09-28 01:00:00 (Sun)
61798471200, #    local_end 1959-04-26 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61798500000, #    utc_start 1959-04-26 10:00:00 (Sun)
61811802000, #      utc_end 1959-09-27 09:00:00 (Sun)
61798474800, #  local_start 1959-04-26 03:00:00 (Sun)
61811776800, #    local_end 1959-09-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61811802000, #    utc_start 1959-09-27 09:00:00 (Sun)
61829949600, #      utc_end 1960-04-24 10:00:00 (Sun)
61811773200, #  local_start 1959-09-27 01:00:00 (Sun)
61829920800, #    local_end 1960-04-24 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61829949600, #    utc_start 1960-04-24 10:00:00 (Sun)
61843251600, #      utc_end 1960-09-25 09:00:00 (Sun)
61829924400, #  local_start 1960-04-24 03:00:00 (Sun)
61843226400, #    local_end 1960-09-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
61843251600, #    utc_start 1960-09-25 09:00:00 (Sun)
61851715200, #      utc_end 1961-01-01 08:00:00 (Sun)
61843222800, #  local_start 1960-09-25 01:00:00 (Sun)
61851686400, #    local_end 1961-01-01 00:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
61851715200, #    utc_start 1961-01-01 08:00:00 (Sun)
62325014400, #      utc_end 1976-01-01 08:00:00 (Thu)
61851686400, #  local_start 1961-01-01 00:00:00 (Sun)
62324985600, #    local_end 1976-01-01 00:00:00 (Thu)
-28800,
0,
'PST',
    ],
    [
62325014400, #    utc_start 1976-01-01 08:00:00 (Thu)
62334957600, #      utc_end 1976-04-25 10:00:00 (Sun)
62324985600, #  local_start 1976-01-01 00:00:00 (Thu)
62334928800, #    local_end 1976-04-25 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62334957600, #    utc_start 1976-04-25 10:00:00 (Sun)
62351283600, #      utc_end 1976-10-31 09:00:00 (Sun)
62334932400, #  local_start 1976-04-25 03:00:00 (Sun)
62351258400, #    local_end 1976-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62351283600, #    utc_start 1976-10-31 09:00:00 (Sun)
62366407200, #      utc_end 1977-04-24 10:00:00 (Sun)
62351254800, #  local_start 1976-10-31 01:00:00 (Sun)
62366378400, #    local_end 1977-04-24 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62366407200, #    utc_start 1977-04-24 10:00:00 (Sun)
62382733200, #      utc_end 1977-10-30 09:00:00 (Sun)
62366382000, #  local_start 1977-04-24 03:00:00 (Sun)
62382708000, #    local_end 1977-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62382733200, #    utc_start 1977-10-30 09:00:00 (Sun)
62398461600, #      utc_end 1978-04-30 10:00:00 (Sun)
62382704400, #  local_start 1977-10-30 01:00:00 (Sun)
62398432800, #    local_end 1978-04-30 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62398461600, #    utc_start 1978-04-30 10:00:00 (Sun)
62414182800, #      utc_end 1978-10-29 09:00:00 (Sun)
62398436400, #  local_start 1978-04-30 03:00:00 (Sun)
62414157600, #    local_end 1978-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62414182800, #    utc_start 1978-10-29 09:00:00 (Sun)
62429911200, #      utc_end 1979-04-29 10:00:00 (Sun)
62414154000, #  local_start 1978-10-29 01:00:00 (Sun)
62429882400, #    local_end 1979-04-29 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62429911200, #    utc_start 1979-04-29 10:00:00 (Sun)
62445632400, #      utc_end 1979-10-28 09:00:00 (Sun)
62429886000, #  local_start 1979-04-29 03:00:00 (Sun)
62445607200, #    local_end 1979-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62445632400, #    utc_start 1979-10-28 09:00:00 (Sun)
62461360800, #      utc_end 1980-04-27 10:00:00 (Sun)
62445603600, #  local_start 1979-10-28 01:00:00 (Sun)
62461332000, #    local_end 1980-04-27 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62461360800, #    utc_start 1980-04-27 10:00:00 (Sun)
62477082000, #      utc_end 1980-10-26 09:00:00 (Sun)
62461335600, #  local_start 1980-04-27 03:00:00 (Sun)
62477056800, #    local_end 1980-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62477082000, #    utc_start 1980-10-26 09:00:00 (Sun)
62492810400, #      utc_end 1981-04-26 10:00:00 (Sun)
62477053200, #  local_start 1980-10-26 01:00:00 (Sun)
62492781600, #    local_end 1981-04-26 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62492810400, #    utc_start 1981-04-26 10:00:00 (Sun)
62508531600, #      utc_end 1981-10-25 09:00:00 (Sun)
62492785200, #  local_start 1981-04-26 03:00:00 (Sun)
62508506400, #    local_end 1981-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62508531600, #    utc_start 1981-10-25 09:00:00 (Sun)
62524260000, #      utc_end 1982-04-25 10:00:00 (Sun)
62508502800, #  local_start 1981-10-25 01:00:00 (Sun)
62524231200, #    local_end 1982-04-25 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62524260000, #    utc_start 1982-04-25 10:00:00 (Sun)
62540586000, #      utc_end 1982-10-31 09:00:00 (Sun)
62524234800, #  local_start 1982-04-25 03:00:00 (Sun)
62540560800, #    local_end 1982-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62540586000, #    utc_start 1982-10-31 09:00:00 (Sun)
62555709600, #      utc_end 1983-04-24 10:00:00 (Sun)
62540557200, #  local_start 1982-10-31 01:00:00 (Sun)
62555680800, #    local_end 1983-04-24 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62555709600, #    utc_start 1983-04-24 10:00:00 (Sun)
62572035600, #      utc_end 1983-10-30 09:00:00 (Sun)
62555684400, #  local_start 1983-04-24 03:00:00 (Sun)
62572010400, #    local_end 1983-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62572035600, #    utc_start 1983-10-30 09:00:00 (Sun)
62587764000, #      utc_end 1984-04-29 10:00:00 (Sun)
62572006800, #  local_start 1983-10-30 01:00:00 (Sun)
62587735200, #    local_end 1984-04-29 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62587764000, #    utc_start 1984-04-29 10:00:00 (Sun)
62603485200, #      utc_end 1984-10-28 09:00:00 (Sun)
62587738800, #  local_start 1984-04-29 03:00:00 (Sun)
62603460000, #    local_end 1984-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62603485200, #    utc_start 1984-10-28 09:00:00 (Sun)
62619213600, #      utc_end 1985-04-28 10:00:00 (Sun)
62603456400, #  local_start 1984-10-28 01:00:00 (Sun)
62619184800, #    local_end 1985-04-28 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62619213600, #    utc_start 1985-04-28 10:00:00 (Sun)
62634934800, #      utc_end 1985-10-27 09:00:00 (Sun)
62619188400, #  local_start 1985-04-28 03:00:00 (Sun)
62634909600, #    local_end 1985-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62634934800, #    utc_start 1985-10-27 09:00:00 (Sun)
62650663200, #      utc_end 1986-04-27 10:00:00 (Sun)
62634906000, #  local_start 1985-10-27 01:00:00 (Sun)
62650634400, #    local_end 1986-04-27 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62650663200, #    utc_start 1986-04-27 10:00:00 (Sun)
62666384400, #      utc_end 1986-10-26 09:00:00 (Sun)
62650638000, #  local_start 1986-04-27 03:00:00 (Sun)
62666359200, #    local_end 1986-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62666384400, #    utc_start 1986-10-26 09:00:00 (Sun)
62680298400, #      utc_end 1987-04-05 10:00:00 (Sun)
62666355600, #  local_start 1986-10-26 01:00:00 (Sun)
62680269600, #    local_end 1987-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62680298400, #    utc_start 1987-04-05 10:00:00 (Sun)
62697834000, #      utc_end 1987-10-25 09:00:00 (Sun)
62680273200, #  local_start 1987-04-05 03:00:00 (Sun)
62697808800, #    local_end 1987-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62697834000, #    utc_start 1987-10-25 09:00:00 (Sun)
62711748000, #      utc_end 1988-04-03 10:00:00 (Sun)
62697805200, #  local_start 1987-10-25 01:00:00 (Sun)
62711719200, #    local_end 1988-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62711748000, #    utc_start 1988-04-03 10:00:00 (Sun)
62729888400, #      utc_end 1988-10-30 09:00:00 (Sun)
62711722800, #  local_start 1988-04-03 03:00:00 (Sun)
62729863200, #    local_end 1988-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62729888400, #    utc_start 1988-10-30 09:00:00 (Sun)
62743197600, #      utc_end 1989-04-02 10:00:00 (Sun)
62729859600, #  local_start 1988-10-30 01:00:00 (Sun)
62743168800, #    local_end 1989-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62743197600, #    utc_start 1989-04-02 10:00:00 (Sun)
62761338000, #      utc_end 1989-10-29 09:00:00 (Sun)
62743172400, #  local_start 1989-04-02 03:00:00 (Sun)
62761312800, #    local_end 1989-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62761338000, #    utc_start 1989-10-29 09:00:00 (Sun)
62774647200, #      utc_end 1990-04-01 10:00:00 (Sun)
62761309200, #  local_start 1989-10-29 01:00:00 (Sun)
62774618400, #    local_end 1990-04-01 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62774647200, #    utc_start 1990-04-01 10:00:00 (Sun)
62792787600, #      utc_end 1990-10-28 09:00:00 (Sun)
62774622000, #  local_start 1990-04-01 03:00:00 (Sun)
62792762400, #    local_end 1990-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62792787600, #    utc_start 1990-10-28 09:00:00 (Sun)
62806701600, #      utc_end 1991-04-07 10:00:00 (Sun)
62792758800, #  local_start 1990-10-28 01:00:00 (Sun)
62806672800, #    local_end 1991-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62806701600, #    utc_start 1991-04-07 10:00:00 (Sun)
62824237200, #      utc_end 1991-10-27 09:00:00 (Sun)
62806676400, #  local_start 1991-04-07 03:00:00 (Sun)
62824212000, #    local_end 1991-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62824237200, #    utc_start 1991-10-27 09:00:00 (Sun)
62838151200, #      utc_end 1992-04-05 10:00:00 (Sun)
62824208400, #  local_start 1991-10-27 01:00:00 (Sun)
62838122400, #    local_end 1992-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62838151200, #    utc_start 1992-04-05 10:00:00 (Sun)
62855686800, #      utc_end 1992-10-25 09:00:00 (Sun)
62838126000, #  local_start 1992-04-05 03:00:00 (Sun)
62855661600, #    local_end 1992-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62855686800, #    utc_start 1992-10-25 09:00:00 (Sun)
62869600800, #      utc_end 1993-04-04 10:00:00 (Sun)
62855658000, #  local_start 1992-10-25 01:00:00 (Sun)
62869572000, #    local_end 1993-04-04 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62869600800, #    utc_start 1993-04-04 10:00:00 (Sun)
62887741200, #      utc_end 1993-10-31 09:00:00 (Sun)
62869575600, #  local_start 1993-04-04 03:00:00 (Sun)
62887716000, #    local_end 1993-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62887741200, #    utc_start 1993-10-31 09:00:00 (Sun)
62901050400, #      utc_end 1994-04-03 10:00:00 (Sun)
62887712400, #  local_start 1993-10-31 01:00:00 (Sun)
62901021600, #    local_end 1994-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62901050400, #    utc_start 1994-04-03 10:00:00 (Sun)
62919190800, #      utc_end 1994-10-30 09:00:00 (Sun)
62901025200, #  local_start 1994-04-03 03:00:00 (Sun)
62919165600, #    local_end 1994-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62919190800, #    utc_start 1994-10-30 09:00:00 (Sun)
62932500000, #      utc_end 1995-04-02 10:00:00 (Sun)
62919162000, #  local_start 1994-10-30 01:00:00 (Sun)
62932471200, #    local_end 1995-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62932500000, #    utc_start 1995-04-02 10:00:00 (Sun)
62950640400, #      utc_end 1995-10-29 09:00:00 (Sun)
62932474800, #  local_start 1995-04-02 03:00:00 (Sun)
62950615200, #    local_end 1995-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62950640400, #    utc_start 1995-10-29 09:00:00 (Sun)
62956166400, #      utc_end 1996-01-01 08:00:00 (Mon)
62950611600, #  local_start 1995-10-29 01:00:00 (Sun)
62956137600, #    local_end 1996-01-01 00:00:00 (Mon)
-28800,
0,
'PST',
    ],
    [
62956166400, #    utc_start 1996-01-01 08:00:00 (Mon)
62964554400, #      utc_end 1996-04-07 10:00:00 (Sun)
62956137600, #  local_start 1996-01-01 00:00:00 (Mon)
62964525600, #    local_end 1996-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62964554400, #    utc_start 1996-04-07 10:00:00 (Sun)
62982090000, #      utc_end 1996-10-27 09:00:00 (Sun)
62964529200, #  local_start 1996-04-07 03:00:00 (Sun)
62982064800, #    local_end 1996-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
62982090000, #    utc_start 1996-10-27 09:00:00 (Sun)
62996004000, #      utc_end 1997-04-06 10:00:00 (Sun)
62982061200, #  local_start 1996-10-27 01:00:00 (Sun)
62995975200, #    local_end 1997-04-06 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
62996004000, #    utc_start 1997-04-06 10:00:00 (Sun)
63013539600, #      utc_end 1997-10-26 09:00:00 (Sun)
62995978800, #  local_start 1997-04-06 03:00:00 (Sun)
63013514400, #    local_end 1997-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63013539600, #    utc_start 1997-10-26 09:00:00 (Sun)
63027453600, #      utc_end 1998-04-05 10:00:00 (Sun)
63013510800, #  local_start 1997-10-26 01:00:00 (Sun)
63027424800, #    local_end 1998-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63027453600, #    utc_start 1998-04-05 10:00:00 (Sun)
63044989200, #      utc_end 1998-10-25 09:00:00 (Sun)
63027428400, #  local_start 1998-04-05 03:00:00 (Sun)
63044964000, #    local_end 1998-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63044989200, #    utc_start 1998-10-25 09:00:00 (Sun)
63058903200, #      utc_end 1999-04-04 10:00:00 (Sun)
63044960400, #  local_start 1998-10-25 01:00:00 (Sun)
63058874400, #    local_end 1999-04-04 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63058903200, #    utc_start 1999-04-04 10:00:00 (Sun)
63077043600, #      utc_end 1999-10-31 09:00:00 (Sun)
63058878000, #  local_start 1999-04-04 03:00:00 (Sun)
63077018400, #    local_end 1999-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63077043600, #    utc_start 1999-10-31 09:00:00 (Sun)
63090352800, #      utc_end 2000-04-02 10:00:00 (Sun)
63077014800, #  local_start 1999-10-31 01:00:00 (Sun)
63090324000, #    local_end 2000-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63090352800, #    utc_start 2000-04-02 10:00:00 (Sun)
63108493200, #      utc_end 2000-10-29 09:00:00 (Sun)
63090327600, #  local_start 2000-04-02 03:00:00 (Sun)
63108468000, #    local_end 2000-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63108493200, #    utc_start 2000-10-29 09:00:00 (Sun)
63114019200, #      utc_end 2001-01-01 08:00:00 (Mon)
63108464400, #  local_start 2000-10-29 01:00:00 (Sun)
63113990400, #    local_end 2001-01-01 00:00:00 (Mon)
-28800,
0,
'PST',
    ],
    [
63114019200, #    utc_start 2001-01-01 08:00:00 (Mon)
63121802400, #      utc_end 2001-04-01 10:00:00 (Sun)
63113990400, #  local_start 2001-01-01 00:00:00 (Mon)
63121773600, #    local_end 2001-04-01 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63121802400, #    utc_start 2001-04-01 10:00:00 (Sun)
63139942800, #      utc_end 2001-10-28 09:00:00 (Sun)
63121777200, #  local_start 2001-04-01 03:00:00 (Sun)
63139917600, #    local_end 2001-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63139942800, #    utc_start 2001-10-28 09:00:00 (Sun)
63149875200, #      utc_end 2002-02-20 08:00:00 (Wed)
63139914000, #  local_start 2001-10-28 01:00:00 (Sun)
63149846400, #    local_end 2002-02-20 00:00:00 (Wed)
-28800,
0,
'PST',
    ],
    [
63149875200, #    utc_start 2002-02-20 08:00:00 (Wed)
63153856800, #      utc_end 2002-04-07 10:00:00 (Sun)
63149846400, #  local_start 2002-02-20 00:00:00 (Wed)
63153828000, #    local_end 2002-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63153856800, #    utc_start 2002-04-07 10:00:00 (Sun)
63171392400, #      utc_end 2002-10-27 09:00:00 (Sun)
63153831600, #  local_start 2002-04-07 03:00:00 (Sun)
63171367200, #    local_end 2002-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63171392400, #    utc_start 2002-10-27 09:00:00 (Sun)
63185306400, #      utc_end 2003-04-06 10:00:00 (Sun)
63171363600, #  local_start 2002-10-27 01:00:00 (Sun)
63185277600, #    local_end 2003-04-06 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63185306400, #    utc_start 2003-04-06 10:00:00 (Sun)
63202842000, #      utc_end 2003-10-26 09:00:00 (Sun)
63185281200, #  local_start 2003-04-06 03:00:00 (Sun)
63202816800, #    local_end 2003-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63202842000, #    utc_start 2003-10-26 09:00:00 (Sun)
63216756000, #      utc_end 2004-04-04 10:00:00 (Sun)
63202813200, #  local_start 2003-10-26 01:00:00 (Sun)
63216727200, #    local_end 2004-04-04 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63216756000, #    utc_start 2004-04-04 10:00:00 (Sun)
63234896400, #      utc_end 2004-10-31 09:00:00 (Sun)
63216730800, #  local_start 2004-04-04 03:00:00 (Sun)
63234871200, #    local_end 2004-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63234896400, #    utc_start 2004-10-31 09:00:00 (Sun)
63248205600, #      utc_end 2005-04-03 10:00:00 (Sun)
63234867600, #  local_start 2004-10-31 01:00:00 (Sun)
63248176800, #    local_end 2005-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63248205600, #    utc_start 2005-04-03 10:00:00 (Sun)
63266346000, #      utc_end 2005-10-30 09:00:00 (Sun)
63248180400, #  local_start 2005-04-03 03:00:00 (Sun)
63266320800, #    local_end 2005-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63266346000, #    utc_start 2005-10-30 09:00:00 (Sun)
63279655200, #      utc_end 2006-04-02 10:00:00 (Sun)
63266317200, #  local_start 2005-10-30 01:00:00 (Sun)
63279626400, #    local_end 2006-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63279655200, #    utc_start 2006-04-02 10:00:00 (Sun)
63297795600, #      utc_end 2006-10-29 09:00:00 (Sun)
63279630000, #  local_start 2006-04-02 03:00:00 (Sun)
63297770400, #    local_end 2006-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63297795600, #    utc_start 2006-10-29 09:00:00 (Sun)
63311104800, #      utc_end 2007-04-01 10:00:00 (Sun)
63297766800, #  local_start 2006-10-29 01:00:00 (Sun)
63311076000, #    local_end 2007-04-01 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63311104800, #    utc_start 2007-04-01 10:00:00 (Sun)
63329245200, #      utc_end 2007-10-28 09:00:00 (Sun)
63311079600, #  local_start 2007-04-01 03:00:00 (Sun)
63329220000, #    local_end 2007-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63329245200, #    utc_start 2007-10-28 09:00:00 (Sun)
63343159200, #      utc_end 2008-04-06 10:00:00 (Sun)
63329216400, #  local_start 2007-10-28 01:00:00 (Sun)
63343130400, #    local_end 2008-04-06 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63343159200, #    utc_start 2008-04-06 10:00:00 (Sun)
63360694800, #      utc_end 2008-10-26 09:00:00 (Sun)
63343134000, #  local_start 2008-04-06 03:00:00 (Sun)
63360669600, #    local_end 2008-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63360694800, #    utc_start 2008-10-26 09:00:00 (Sun)
63374608800, #      utc_end 2009-04-05 10:00:00 (Sun)
63360666000, #  local_start 2008-10-26 01:00:00 (Sun)
63374580000, #    local_end 2009-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63374608800, #    utc_start 2009-04-05 10:00:00 (Sun)
63392144400, #      utc_end 2009-10-25 09:00:00 (Sun)
63374583600, #  local_start 2009-04-05 03:00:00 (Sun)
63392119200, #    local_end 2009-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63392144400, #    utc_start 2009-10-25 09:00:00 (Sun)
63406058400, #      utc_end 2010-04-04 10:00:00 (Sun)
63392115600, #  local_start 2009-10-25 01:00:00 (Sun)
63406029600, #    local_end 2010-04-04 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63406058400, #    utc_start 2010-04-04 10:00:00 (Sun)
63424198800, #      utc_end 2010-10-31 09:00:00 (Sun)
63406033200, #  local_start 2010-04-04 03:00:00 (Sun)
63424173600, #    local_end 2010-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63424198800, #    utc_start 2010-10-31 09:00:00 (Sun)
63437508000, #      utc_end 2011-04-03 10:00:00 (Sun)
63424170000, #  local_start 2010-10-31 01:00:00 (Sun)
63437479200, #    local_end 2011-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63437508000, #    utc_start 2011-04-03 10:00:00 (Sun)
63455648400, #      utc_end 2011-10-30 09:00:00 (Sun)
63437482800, #  local_start 2011-04-03 03:00:00 (Sun)
63455623200, #    local_end 2011-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63455648400, #    utc_start 2011-10-30 09:00:00 (Sun)
63468957600, #      utc_end 2012-04-01 10:00:00 (Sun)
63455619600, #  local_start 2011-10-30 01:00:00 (Sun)
63468928800, #    local_end 2012-04-01 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63468957600, #    utc_start 2012-04-01 10:00:00 (Sun)
63487098000, #      utc_end 2012-10-28 09:00:00 (Sun)
63468932400, #  local_start 2012-04-01 03:00:00 (Sun)
63487072800, #    local_end 2012-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63487098000, #    utc_start 2012-10-28 09:00:00 (Sun)
63501012000, #      utc_end 2013-04-07 10:00:00 (Sun)
63487069200, #  local_start 2012-10-28 01:00:00 (Sun)
63500983200, #    local_end 2013-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63501012000, #    utc_start 2013-04-07 10:00:00 (Sun)
63518547600, #      utc_end 2013-10-27 09:00:00 (Sun)
63500986800, #  local_start 2013-04-07 03:00:00 (Sun)
63518522400, #    local_end 2013-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63518547600, #    utc_start 2013-10-27 09:00:00 (Sun)
63532461600, #      utc_end 2014-04-06 10:00:00 (Sun)
63518518800, #  local_start 2013-10-27 01:00:00 (Sun)
63532432800, #    local_end 2014-04-06 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63532461600, #    utc_start 2014-04-06 10:00:00 (Sun)
63549997200, #      utc_end 2014-10-26 09:00:00 (Sun)
63532436400, #  local_start 2014-04-06 03:00:00 (Sun)
63549972000, #    local_end 2014-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63549997200, #    utc_start 2014-10-26 09:00:00 (Sun)
63563911200, #      utc_end 2015-04-05 10:00:00 (Sun)
63549968400, #  local_start 2014-10-26 01:00:00 (Sun)
63563882400, #    local_end 2015-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63563911200, #    utc_start 2015-04-05 10:00:00 (Sun)
63581446800, #      utc_end 2015-10-25 09:00:00 (Sun)
63563886000, #  local_start 2015-04-05 03:00:00 (Sun)
63581421600, #    local_end 2015-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63581446800, #    utc_start 2015-10-25 09:00:00 (Sun)
63595360800, #      utc_end 2016-04-03 10:00:00 (Sun)
63581418000, #  local_start 2015-10-25 01:00:00 (Sun)
63595332000, #    local_end 2016-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63595360800, #    utc_start 2016-04-03 10:00:00 (Sun)
63613501200, #      utc_end 2016-10-30 09:00:00 (Sun)
63595335600, #  local_start 2016-04-03 03:00:00 (Sun)
63613476000, #    local_end 2016-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63613501200, #    utc_start 2016-10-30 09:00:00 (Sun)
63626810400, #      utc_end 2017-04-02 10:00:00 (Sun)
63613472400, #  local_start 2016-10-30 01:00:00 (Sun)
63626781600, #    local_end 2017-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63626810400, #    utc_start 2017-04-02 10:00:00 (Sun)
63644950800, #      utc_end 2017-10-29 09:00:00 (Sun)
63626785200, #  local_start 2017-04-02 03:00:00 (Sun)
63644925600, #    local_end 2017-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63644950800, #    utc_start 2017-10-29 09:00:00 (Sun)
63658260000, #      utc_end 2018-04-01 10:00:00 (Sun)
63644922000, #  local_start 2017-10-29 01:00:00 (Sun)
63658231200, #    local_end 2018-04-01 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63658260000, #    utc_start 2018-04-01 10:00:00 (Sun)
63676400400, #      utc_end 2018-10-28 09:00:00 (Sun)
63658234800, #  local_start 2018-04-01 03:00:00 (Sun)
63676375200, #    local_end 2018-10-28 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63676400400, #    utc_start 2018-10-28 09:00:00 (Sun)
63690314400, #      utc_end 2019-04-07 10:00:00 (Sun)
63676371600, #  local_start 2018-10-28 01:00:00 (Sun)
63690285600, #    local_end 2019-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63690314400, #    utc_start 2019-04-07 10:00:00 (Sun)
63707850000, #      utc_end 2019-10-27 09:00:00 (Sun)
63690289200, #  local_start 2019-04-07 03:00:00 (Sun)
63707824800, #    local_end 2019-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63707850000, #    utc_start 2019-10-27 09:00:00 (Sun)
63721764000, #      utc_end 2020-04-05 10:00:00 (Sun)
63707821200, #  local_start 2019-10-27 01:00:00 (Sun)
63721735200, #    local_end 2020-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63721764000, #    utc_start 2020-04-05 10:00:00 (Sun)
63739299600, #      utc_end 2020-10-25 09:00:00 (Sun)
63721738800, #  local_start 2020-04-05 03:00:00 (Sun)
63739274400, #    local_end 2020-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63739299600, #    utc_start 2020-10-25 09:00:00 (Sun)
63753213600, #      utc_end 2021-04-04 10:00:00 (Sun)
63739270800, #  local_start 2020-10-25 01:00:00 (Sun)
63753184800, #    local_end 2021-04-04 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63753213600, #    utc_start 2021-04-04 10:00:00 (Sun)
63771354000, #      utc_end 2021-10-31 09:00:00 (Sun)
63753188400, #  local_start 2021-04-04 03:00:00 (Sun)
63771328800, #    local_end 2021-10-31 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63771354000, #    utc_start 2021-10-31 09:00:00 (Sun)
63784663200, #      utc_end 2022-04-03 10:00:00 (Sun)
63771325200, #  local_start 2021-10-31 01:00:00 (Sun)
63784634400, #    local_end 2022-04-03 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63784663200, #    utc_start 2022-04-03 10:00:00 (Sun)
63802803600, #      utc_end 2022-10-30 09:00:00 (Sun)
63784638000, #  local_start 2022-04-03 03:00:00 (Sun)
63802778400, #    local_end 2022-10-30 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63802803600, #    utc_start 2022-10-30 09:00:00 (Sun)
63816112800, #      utc_end 2023-04-02 10:00:00 (Sun)
63802774800, #  local_start 2022-10-30 01:00:00 (Sun)
63816084000, #    local_end 2023-04-02 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63816112800, #    utc_start 2023-04-02 10:00:00 (Sun)
63834253200, #      utc_end 2023-10-29 09:00:00 (Sun)
63816087600, #  local_start 2023-04-02 03:00:00 (Sun)
63834228000, #    local_end 2023-10-29 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63834253200, #    utc_start 2023-10-29 09:00:00 (Sun)
63848167200, #      utc_end 2024-04-07 10:00:00 (Sun)
63834224400, #  local_start 2023-10-29 01:00:00 (Sun)
63848138400, #    local_end 2024-04-07 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63848167200, #    utc_start 2024-04-07 10:00:00 (Sun)
63865702800, #      utc_end 2024-10-27 09:00:00 (Sun)
63848142000, #  local_start 2024-04-07 03:00:00 (Sun)
63865677600, #    local_end 2024-10-27 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63865702800, #    utc_start 2024-10-27 09:00:00 (Sun)
63879616800, #      utc_end 2025-04-06 10:00:00 (Sun)
63865674000, #  local_start 2024-10-27 01:00:00 (Sun)
63879588000, #    local_end 2025-04-06 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63879616800, #    utc_start 2025-04-06 10:00:00 (Sun)
63897152400, #      utc_end 2025-10-26 09:00:00 (Sun)
63879591600, #  local_start 2025-04-06 03:00:00 (Sun)
63897127200, #    local_end 2025-10-26 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
    [
63897152400, #    utc_start 2025-10-26 09:00:00 (Sun)
63911066400, #      utc_end 2026-04-05 10:00:00 (Sun)
63897123600, #  local_start 2025-10-26 01:00:00 (Sun)
63911037600, #    local_end 2026-04-05 02:00:00 (Sun)
-28800,
0,
'PST',
    ],
    [
63911066400, #    utc_start 2026-04-05 10:00:00 (Sun)
63928602000, #      utc_end 2026-10-25 09:00:00 (Sun)
63911041200, #  local_start 2026-04-05 03:00:00 (Sun)
63928576800, #    local_end 2026-10-25 02:00:00 (Sun)
-25200,
1,
'PDT',
    ],
];

sub olson_version {'2015a'}

sub has_dst_changes {62}

sub _max_year {2025}

sub _new_instance {
    return shift->_init( @_, spans => $spans );
}

sub _last_offset { -28800 }

my $last_observance = bless( {
  'format' => 'P%sT',
  'gmtoff' => '-8:00',
  'local_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 730901,
    'local_rd_secs' => 0,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 730901,
    'utc_rd_secs' => 0,
    'utc_year' => 2003
  }, 'DateTime' ),
  'offset_from_std' => 0,
  'offset_from_utc' => -28800,
  'until' => [],
  'utc_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 730901,
    'local_rd_secs' => 28800,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 730901,
    'utc_rd_secs' => 28800,
    'utc_year' => 2003
  }, 'DateTime' )
}, 'DateTime::TimeZone::OlsonDB::Observance' )
;
sub _last_observance { $last_observance }

my $rules = [
  bless( {
    'at' => '2:00',
    'from' => '2002',
    'in' => 'Oct',
    'letter' => 'S',
    'name' => 'Mexico',
    'offset_from_std' => 0,
    'on' => 'lastSun',
    'save' => '0',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' ),
  bless( {
    'at' => '2:00',
    'from' => '2002',
    'in' => 'Apr',
    'letter' => 'D',
    'name' => 'Mexico',
    'offset_from_std' => 3600,
    'on' => 'Sun>=1',
    'save' => '1:00',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' )
]
;
sub _rules { $rules }


1;
