// -----------------------------------------------------------------------------
// datetime: because sometimes my own time is much easier to manage!
// -----------------------------------------------------------------------------
use IO.FormattedIO;
// -----------------------------------------------------------------------------
// number of days in each month, non-leap/leap years
// -----------------------------------------------------------------------------
const NNdays:
[0..1,1..12] int = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31;
                    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] ;
// -----------------------------------------------------------------------------
// --> leap: is this a leap year?
// -----------------------------------------------------------------------------
inline proc leap(
   const in year: int): bool {
   return ((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0) ;
}
// -----------------------------------------------------------------------------
// --> ymdh2str: convert ymdh to string
// -----------------------------------------------------------------------------
inline proc ymdh2str(
   const in year: int,
   const in month: int,
   const in day: int,
   const in hour: int
   ): string {
   var s = "%4i-%02i-%02i:%02i".format(year,month,day,hour);
   return s;
}
// -----------------------------------------------------------------------------
// --> ymd2str: convert ymdh to string
// -----------------------------------------------------------------------------
inline proc ymd2str(
   const in year: int,
   const in month: int,
   const in day: int
   ): string {
   var s = "%4i-%02i-%02i".format(year,month,day);
   return s;
}



