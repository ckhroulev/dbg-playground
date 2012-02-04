#include "drainagebasin.hh"
#include <string.h>

#include <gsl/gsl_matrix.h>

void init_mask(int Mx, int My,
               double *thickness,
               double *mask,
               double *tmp) {

  memset(mask, 0, Mx*My*sizeof(double));

#define THK(i, j)  thickness[(j)*Mx + (i)]
#define MASK(i, j) mask[(j)*Mx + (i)]
#define TMP(i, j)  tmp[(j)*Mx + (i)]

  double thk_eps = 1;

  int marker = 1;

  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {

      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        TMP(i, j) = ICE_FREE;
        continue;
      }

      double thk = THK(i,     j),
        thk_w    = THK(i - 1, j),
        thk_nw   = THK(i - 1, j + 1),
        thk_n    = THK(i,     j + 1),
        thk_ne   = THK(i + 1, j + 1),
        thk_e    = THK(i + 1, j),
        thk_se   = THK(i + 1, j - 1),
        thk_s    = THK(i,     j - 1),
        thk_sw   = THK(i - 1, j - 1);

      if (thk > thk_eps) {
        // icy cell

        if (thk_w <= thk_eps || thk_nw <= thk_eps || thk_n <= thk_eps || thk_ne <= thk_eps ||
            thk_e <= thk_eps || thk_se <= thk_eps || thk_s <= thk_eps || thk_sw <= thk_eps) {
          // ice margin
          TMP(i, j) = ++marker;
        } else {
          // interior ice
          TMP(i, j) = NO_VALUE;
        }

      } else {
        // ice-free
        TMP(i, j) = ICE_FREE;
      }

    } // inner for loop
  } // outer for loop

  // second pass
  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {

      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        MASK(i, j) = TMP(i, j);
        continue;
      }

      double m  = TMP(i,     j),
        mask_w  = TMP(i - 1, j),
        mask_nw = TMP(i - 1, j + 1),
        mask_n  = TMP(i,     j + 1),
        mask_ne = TMP(i + 1, j + 1),
        mask_e  = TMP(i + 1, j),
        mask_se = TMP(i + 1, j - 1),
        mask_s  = TMP(i,     j - 1),
        mask_sw = TMP(i - 1, j - 1);

      double thk = THK(i,     j),
        thk_w    = THK(i - 1, j),
        thk_nw   = THK(i - 1, j + 1),
        thk_n    = THK(i,     j + 1),
        thk_ne   = THK(i + 1, j + 1),
        thk_e    = THK(i + 1, j),
        thk_se   = THK(i + 1, j - 1),
        thk_s    = THK(i,     j - 1),
        thk_sw   = THK(i - 1, j - 1);

      // ice-free cell next to an icy cell
      if (thk < thk_eps &&
          (thk_w >= thk_eps || thk_nw >= thk_eps || thk_n >= thk_eps || thk_ne >= thk_eps ||
           thk_e >= thk_eps || thk_se >= thk_eps || thk_s >= thk_eps || thk_sw >= thk_eps)) {

        if (mask_w > 0)
          MASK(i, j) =  mask_w;
        else if (mask_nw > 0)
          MASK(i, j) =  mask_nw;
        else if (mask_n > 0)
          MASK(i, j) =  mask_n;
        else if (mask_ne > 0)
          MASK(i, j) =  mask_ne;
        else if (mask_e > 0)
          MASK(i, j) =  mask_e;
        else if (mask_se > 0)
          MASK(i, j) =  mask_se;
        else if (mask_s > 0)
          MASK(i, j) =  mask_s;
        else if (mask_sw > 0)
          MASK(i, j) =  mask_sw;
      } else {
        MASK(i, j) =  m;
      }

    } // inner for loop
  } // outer for loop

#undef THK(i,j)
#undef MASK(i,j)
#undef TMP(i,j)

  memcpy(tmp, mask, Mx*My*sizeof(double));
}
