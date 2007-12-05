/* Chooses the best routine between Karatsuba and Toom-Cook variants.

  Copyright 2007 Richard P. Brent.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

#ifndef MUL_KARA_THRESHOLD
#define MUL_KARA_THRESHOLD 8
#endif

#ifndef MUL_TOOM_THRESHOLD
#define MUL_TOOM_THRESHOLD 17
#endif

#ifndef MUL_TOOMW_THRESHOLD
#define MUL_TOOMW_THRESHOLD 8
#endif

#ifndef MUL_TOOMU_THRESHOLD
#define MUL_TOOMU_THRESHOLD 33
#endif

#ifndef MUL_TOOM4_THRESHOLD
#define MUL_TOOM4_THRESHOLD 30
#endif

#ifndef TUNING
 
#ifdef BEST_TOOM_TABLE
  short best_tab[] = BEST_TOOM_TABLE;
#else
  short best_tab[] = {1,1,1,1,1,1,1,1};
#endif

#ifdef BEST_UTOOM_TABLE
  short best_utab[] = BEST_UTOOM_TABLE;
#else
  short best_utab[] = {0}
#endif

#endif    

/* Returns 1 for KarMul, 2 for Toom3Mul, 3 for Toom3WMul, 4 for Toom4Mul
   depending on which is predicted to be fastest for the given degree n.

   RPB, 20070511 */

static inline long BestToom (long n)

{
// BEST_TOOM_TABLE should be generated by the tuning program tunetoom.
//
// The n-th entry in the list gives the code for the fastest algorithm for
// input size n.  For example: 
// #define BEST_TOOM_TABLE {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,3,2,3,3,2,3}
// would be reasonable if Toom3Mul was fastest for n = 18, 21, 24.

  if (n < MUL_TOOMW_THRESHOLD)
    return 1;					// KarMul

  if (n <= sizeof(best_tab)/sizeof(short))	// In range of table
   return (long)best_tab[n-1];			// Return table entry
   
// Here we are outside the bounds of the the table best_tab so return the
// largest feasible value

  if (n >= MUL_TOOM4_THRESHOLD)
    return 4;					// Toom4Mul
  else
    return 3;					// Toom3WMul
}

static inline long BestuToom (long n)

{
// BEST_UTOOM_TABLE should be generated by the tuning program tuneutoom.
//
// The n-th entry in the list gives the code for the fastest algorithm for
// input size n.  0 means the obvious splitting algorithm and 1 means
// Toom3uMul.

  if (n < MUL_TOOMU_THRESHOLD)
    return 0;					// Default

  if (n <= sizeof(best_utab)/sizeof(short))	// In range of table
   return (long)best_utab[n-1];			// Return table entry
   
// Here we are outside the bounds of the the table best_tab so return 1

    return 1;					// Toom3uMul
}
