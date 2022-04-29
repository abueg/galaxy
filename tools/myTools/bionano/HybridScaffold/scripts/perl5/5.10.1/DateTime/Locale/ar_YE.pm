###########################################################################
#
# This file is auto-generated by the Perl DateTime Suite locale
# generator (0.05).  This code generator comes with the
# DateTime::Locale distribution in the tools/ directory, and is called
# generate-from-cldr.
#
# This file as generated from the CLDR XML locale data.  See the
# LICENSE.cldr file included in this distribution for license details.
#
# This file was generated from the source file ar_YE.xml
# The source file version number was 1.48, generated on
# 2009/05/05 23:06:34.
#
# Do not edit this file directly.
#
###########################################################################

package DateTime::Locale::ar_YE;

use strict;
use warnings;
use utf8;

use base 'DateTime::Locale::ar';

sub cldr_version { return "1\.7\.1" }

{
    my $day_format_abbreviated = [ "الاثنين", "الثلاثاء", "الأربعاء", "الخميس", "الجمعة", "السبت", "الأحد" ];
    sub day_format_abbreviated { return $day_format_abbreviated }
}

sub day_format_narrow { $_[0]->day_format_abbreviated() }


sub day_stand_alone_abbreviated { $_[0]->day_format_abbreviated() }

{
    my $first_day_of_week = "6";
    sub first_day_of_week { return $first_day_of_week }
}

{
    my $glibc_date_format = "\%d\ \%b\,\ \%Y";
    sub glibc_date_format { return $glibc_date_format }
}

{
    my $glibc_date_1_format = "\%a\ \%b\ \%e\ \%H\:\%M\:\%S\ \%Z\ \%Y";
    sub glibc_date_1_format { return $glibc_date_1_format }
}

{
    my $glibc_datetime_format = "\%d\ \%b\,\ \%Y\ \%Z\ \%I\:\%M\:\%S\ \%p";
    sub glibc_datetime_format { return $glibc_datetime_format }
}

{
    my $glibc_time_format = "\%Z\ \%I\:\%M\:\%S\ ";
    sub glibc_time_format { return $glibc_time_format }
}

{
    my $glibc_time_12_format = "\%Z\ \%I\:\%M\:\%S\ \%p";
    sub glibc_time_12_format { return $glibc_time_12_format }
}

1;

__END__


=pod

=encoding utf8

=head1 NAME

DateTime::Locale::ar_YE

=head1 SYNOPSIS

  use DateTime;

  my $dt = DateTime->now( locale => 'ar_YE' );
  print $dt->month_name();

=head1 DESCRIPTION

This is the DateTime locale package for Arabic Yemen.

=head1 DATA

This locale inherits from the L<DateTime::Locale::ar> locale.

It contains the following data.

=head2 Days

=head3 Wide (format)

  الإثنين
  الثلاثاء
  الأربعاء
  الخميس
  الجمعة
  السبت
  الأحد

=head3 Abbreviated (format)

  الاثنين
  الثلاثاء
  الأربعاء
  الخميس
  الجمعة
  السبت
  الأحد

=head3 Narrow (format)

  الاثنين
  الثلاثاء
  الأربعاء
  الخميس
  الجمعة
  السبت
  الأحد

=head3 Wide (stand-alone)

  الإثنين
  الثلاثاء
  الأربعاء
  الخميس
  الجمعة
  السبت
  الأحد

=head3 Abbreviated (stand-alone)

  الاثنين
  الثلاثاء
  الأربعاء
  الخميس
  الجمعة
  السبت
  الأحد

=head3 Narrow (stand-alone)

  ن
  ث
  ر
  خ
  ج
  س
  ح

=head2 Months

=head3 Wide (format)

  يناير
  فبراير
  مارس
  أبريل
  مايو
  يونيو
  يوليو
  أغسطس
  سبتمبر
  أكتوبر
  نوفمبر
  ديسمبر

=head3 Abbreviated (format)

  يناير
  فبراير
  مارس
  أبريل
  مايو
  يونيو
  يوليو
  أغسطس
  سبتمبر
  أكتوبر
  نوفمبر
  ديسمبر

=head3 Narrow (format)

  ي
  ف
  م
  أ
  و
  ن
  ل
  غ
  س
  ك
  ب
  د

=head3 Wide (stand-alone)

  يناير
  فبراير
  مارس
  أبريل
  مايو
  يونيو
  يوليو
  أغسطس
  سبتمبر
  أكتوبر
  نوفمبر
  ديسمبر

=head3 Abbreviated (stand-alone)

  يناير
  فبراير
  مارس
  أبريل
  مايو
  يونيو
  يوليو
  أغسطس
  سبتمبر
  أكتوبر
  نوفمبر
  ديسمبر

=head3 Narrow (stand-alone)

  ي
  ف
  م
  أ
  و
  ن
  ل
  غ
  س
  ك
  ب
  د

=head2 Quarters

=head3 Wide (format)

  الربع الأول
  الربع الثاني
  الربع الثالث
  الربع الرابع

=head3 Abbreviated (format)

  الربع الأول
  الربع الثاني
  الربع الثالث
  الربع الرابع

=head3 Narrow (format)

  ١
  ٢
  ٣
  ٤

=head3 Wide (stand-alone)

  الربع الأول
  الربع الثاني
  الربع الثالث
  الربع الرابع

=head3 Abbreviated (stand-alone)

  الربع الأول
  الربع الثاني
  الربع الثالث
  الربع الرابع

=head3 Narrow (stand-alone)

  ١
  ٢
  ٣
  ٤

=head2 Eras

=head3 Wide

  قبل الميلاد
  ميلادي

=head3 Abbreviated

  ق.م
  م

=head3 Narrow

  ق.م
  م

=head2 Date Formats

=head3 Full

   2008-02-05T18:30:30 = الثلاثاء، 5 فبراير، 2008
   1995-12-22T09:05:02 = الجمعة، 22 ديسمبر، 1995
  -0010-09-15T04:44:23 = السبت، 15 سبتمبر، -10

=head3 Long

   2008-02-05T18:30:30 = 5 فبراير، 2008
   1995-12-22T09:05:02 = 22 ديسمبر، 1995
  -0010-09-15T04:44:23 = 15 سبتمبر، -10

=head3 Medium

   2008-02-05T18:30:30 = 05‏/02‏/2008
   1995-12-22T09:05:02 = 22‏/12‏/1995
  -0010-09-15T04:44:23 = 15‏/09‏/-010

=head3 Short

   2008-02-05T18:30:30 = 5‏/2‏/2008
   1995-12-22T09:05:02 = 22‏/12‏/1995
  -0010-09-15T04:44:23 = 15‏/9‏/-010

=head3 Default

   2008-02-05T18:30:30 = 05‏/02‏/2008
   1995-12-22T09:05:02 = 22‏/12‏/1995
  -0010-09-15T04:44:23 = 15‏/09‏/-010

=head2 Time Formats

=head3 Full

   2008-02-05T18:30:30 = UTC 6:30:30 م
   1995-12-22T09:05:02 = UTC 9:05:02 ص
  -0010-09-15T04:44:23 = UTC 4:44:23 ص

=head3 Long

   2008-02-05T18:30:30 = UTC 6:30:30 م
   1995-12-22T09:05:02 = UTC 9:05:02 ص
  -0010-09-15T04:44:23 = UTC 4:44:23 ص

=head3 Medium

   2008-02-05T18:30:30 = 6:30:30 م
   1995-12-22T09:05:02 = 9:05:02 ص
  -0010-09-15T04:44:23 = 4:44:23 ص

=head3 Short

   2008-02-05T18:30:30 = 6:30 م
   1995-12-22T09:05:02 = 9:05 ص
  -0010-09-15T04:44:23 = 4:44 ص

=head3 Default

   2008-02-05T18:30:30 = 6:30:30 م
   1995-12-22T09:05:02 = 9:05:02 ص
  -0010-09-15T04:44:23 = 4:44:23 ص

=head2 Datetime Formats

=head3 Full

   2008-02-05T18:30:30 = الثلاثاء، 5 فبراير، 2008 UTC 6:30:30 م
   1995-12-22T09:05:02 = الجمعة، 22 ديسمبر، 1995 UTC 9:05:02 ص
  -0010-09-15T04:44:23 = السبت، 15 سبتمبر، -10 UTC 4:44:23 ص

=head3 Long

   2008-02-05T18:30:30 = 5 فبراير، 2008 UTC 6:30:30 م
   1995-12-22T09:05:02 = 22 ديسمبر، 1995 UTC 9:05:02 ص
  -0010-09-15T04:44:23 = 15 سبتمبر، -10 UTC 4:44:23 ص

=head3 Medium

   2008-02-05T18:30:30 = 05‏/02‏/2008 6:30:30 م
   1995-12-22T09:05:02 = 22‏/12‏/1995 9:05:02 ص
  -0010-09-15T04:44:23 = 15‏/09‏/-010 4:44:23 ص

=head3 Short

   2008-02-05T18:30:30 = 5‏/2‏/2008 6:30 م
   1995-12-22T09:05:02 = 22‏/12‏/1995 9:05 ص
  -0010-09-15T04:44:23 = 15‏/9‏/-010 4:44 ص

=head3 Default

   2008-02-05T18:30:30 = 05‏/02‏/2008 6:30:30 م
   1995-12-22T09:05:02 = 22‏/12‏/1995 9:05:02 ص
  -0010-09-15T04:44:23 = 15‏/09‏/-010 4:44:23 ص

=head2 Available Formats

=head3 d (d)

   2008-02-05T18:30:30 = 5
   1995-12-22T09:05:02 = 22
  -0010-09-15T04:44:23 = 15

=head3 EEEd (d EEE)

   2008-02-05T18:30:30 = 5 الثلاثاء
   1995-12-22T09:05:02 = 22 الجمعة
  -0010-09-15T04:44:23 = 15 السبت

=head3 Hm (H:mm)

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 9:05
  -0010-09-15T04:44:23 = 4:44

=head3 hm (h:mm a)

   2008-02-05T18:30:30 = 6:30 م
   1995-12-22T09:05:02 = 9:05 ص
  -0010-09-15T04:44:23 = 4:44 ص

=head3 Hms (H:mm:ss)

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 9:05:02
  -0010-09-15T04:44:23 = 4:44:23

=head3 hms (h:mm:ss a)

   2008-02-05T18:30:30 = 6:30:30 م
   1995-12-22T09:05:02 = 9:05:02 ص
  -0010-09-15T04:44:23 = 4:44:23 ص

=head3 M (L)

   2008-02-05T18:30:30 = 2
   1995-12-22T09:05:02 = 12
  -0010-09-15T04:44:23 = 9

=head3 Md (d/‏M)

   2008-02-05T18:30:30 = 5/‏2
   1995-12-22T09:05:02 = 22/‏12
  -0010-09-15T04:44:23 = 15/‏9

=head3 MEd (E، d-M)

   2008-02-05T18:30:30 = الثلاثاء، 5-2
   1995-12-22T09:05:02 = الجمعة، 22-12
  -0010-09-15T04:44:23 = السبت، 15-9

=head3 MMdd (dd‏/MM)

   2008-02-05T18:30:30 = 05‏/02
   1995-12-22T09:05:02 = 22‏/12
  -0010-09-15T04:44:23 = 15‏/09

=head3 MMM (LLL)

   2008-02-05T18:30:30 = فبراير
   1995-12-22T09:05:02 = ديسمبر
  -0010-09-15T04:44:23 = سبتمبر

=head3 MMMd (d MMM)

   2008-02-05T18:30:30 = 5 فبراير
   1995-12-22T09:05:02 = 22 ديسمبر
  -0010-09-15T04:44:23 = 15 سبتمبر

=head3 MMMEd (E d MMM)

   2008-02-05T18:30:30 = الثلاثاء 5 فبراير
   1995-12-22T09:05:02 = الجمعة 22 ديسمبر
  -0010-09-15T04:44:23 = السبت 15 سبتمبر

=head3 MMMMd (d MMMM)

   2008-02-05T18:30:30 = 5 فبراير
   1995-12-22T09:05:02 = 22 ديسمبر
  -0010-09-15T04:44:23 = 15 سبتمبر

=head3 MMMMEd (E d MMMM)

   2008-02-05T18:30:30 = الثلاثاء 5 فبراير
   1995-12-22T09:05:02 = الجمعة 22 ديسمبر
  -0010-09-15T04:44:23 = السبت 15 سبتمبر

=head3 ms (mm:ss)

   2008-02-05T18:30:30 = 30:30
   1995-12-22T09:05:02 = 05:02
  -0010-09-15T04:44:23 = 44:23

=head3 y (y)

   2008-02-05T18:30:30 = 2008
   1995-12-22T09:05:02 = 1995
  -0010-09-15T04:44:23 = -10

=head3 yM (M‏/yyyy)

   2008-02-05T18:30:30 = 2‏/2008
   1995-12-22T09:05:02 = 12‏/1995
  -0010-09-15T04:44:23 = 9‏/-010

=head3 yMEd (EEE، d/‏M/‏yyyy)

   2008-02-05T18:30:30 = الثلاثاء، 5/‏2/‏2008
   1995-12-22T09:05:02 = الجمعة، 22/‏12/‏1995
  -0010-09-15T04:44:23 = السبت، 15/‏9/‏-010

=head3 yMMM (MMM y)

   2008-02-05T18:30:30 = فبراير 2008
   1995-12-22T09:05:02 = ديسمبر 1995
  -0010-09-15T04:44:23 = سبتمبر -10

=head3 yMMMEd (EEE، d MMMM y)

   2008-02-05T18:30:30 = الثلاثاء، 5 فبراير 2008
   1995-12-22T09:05:02 = الجمعة، 22 ديسمبر 1995
  -0010-09-15T04:44:23 = السبت، 15 سبتمبر -10

=head3 yMMMM (MMMM y)

   2008-02-05T18:30:30 = فبراير 2008
   1995-12-22T09:05:02 = ديسمبر 1995
  -0010-09-15T04:44:23 = سبتمبر -10

=head3 yQ (yyyy Q)

   2008-02-05T18:30:30 = 2008 1
   1995-12-22T09:05:02 = 1995 4
  -0010-09-15T04:44:23 = -010 3

=head3 yQQQ (y QQQ)

   2008-02-05T18:30:30 = 2008 الربع الأول
   1995-12-22T09:05:02 = 1995 الربع الرابع
  -0010-09-15T04:44:23 = -10 الربع الثالث

=head3 yyQ (Q yy)

   2008-02-05T18:30:30 = 1 08
   1995-12-22T09:05:02 = 4 95
  -0010-09-15T04:44:23 = 3 -10

=head3 yyyyMM (MM‏/yyyy)

   2008-02-05T18:30:30 = 02‏/2008
   1995-12-22T09:05:02 = 12‏/1995
  -0010-09-15T04:44:23 = 09‏/-010

=head3 yyyyMMMM (MMMM، y)

   2008-02-05T18:30:30 = فبراير، 2008
   1995-12-22T09:05:02 = ديسمبر، 1995
  -0010-09-15T04:44:23 = سبتمبر، -10

=head2 Miscellaneous

=head3 Prefers 24 hour time?

No

=head3 Local first day of the week

السبت


=head1 SUPPORT

See L<DateTime::Locale>.

=head1 AUTHOR

Dave Rolsky <autarch@urth.org>

=head1 COPYRIGHT

Copyright (c) 2008 David Rolsky. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same
terms as Perl itself.

This module was generated from data provided by the CLDR project, see
the LICENSE.cldr in this distribution for details on the CLDR data's
license.

=cut
