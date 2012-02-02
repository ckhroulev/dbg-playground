#include "drainagebasin.hh"
#include <string.h>

void init_mask(int Mx, int My,
               double *thickness,
               double *mask,
               double *new_mask) {

  memset(mask, 0, Mx*My*sizeof(double));

  gsl_matrix_view mask_view = gsl_matrix_view_array(mask, Mx, My);
  gsl_matrix * m_mask = &mask_view.matrix;

  gsl_matrix_view tmp_view = gsl_matrix_view_array(new_mask, Mx, My);
  gsl_matrix * m_tmp = &tmp_view.matrix;

  gsl_matrix_view thickness_view = gsl_matrix_view_array(thickness, Mx, My);
  gsl_matrix * m_thickness = &thickness_view.matrix;

  double thk_eps = 1;

  int marker = 1;

  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {

      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        gsl_matrix_set(m_tmp, i, j, -1);
        continue;
      }

      double thk = gsl_matrix_get(m_thickness, i, j),
        thk_w = gsl_matrix_get(m_thickness, i - 1, j),
        thk_nw = gsl_matrix_get(m_thickness, i - 1, j + 1),
        thk_n = gsl_matrix_get(m_thickness, i, j + 1),
        thk_ne = gsl_matrix_get(m_thickness, i + 1, j + 1),
        thk_e = gsl_matrix_get(m_thickness, i + 1, j),
        thk_se = gsl_matrix_get(m_thickness, i + 1, j - 1),
        thk_s = gsl_matrix_get(m_thickness, i, j - 1),
        thk_sw = gsl_matrix_get(m_thickness, i - 1, j - 1);

      if (thk > thk_eps) {
        // icy cell

        if (thk_w <= thk_eps || thk_nw <= thk_eps || thk_n <= thk_eps || thk_ne <= thk_eps ||
            thk_e <= thk_eps || thk_se <= thk_eps || thk_s <= thk_eps || thk_sw <= thk_eps) {
          // ice margin
          gsl_matrix_set(m_tmp, i, j, ++marker);
        } else {
          // interior ice
          gsl_matrix_set(m_tmp, i, j, -1);
        }

      } else {
        // ice-free
        gsl_matrix_set(m_tmp, i, j, -2);
      }


    } // inner for loop
  } // outer for loop

  // second pass
  for (int i = 1; i < Mx-1; ++i) {
    for (int j = 1; j < My-1; ++j) {

      double m = gsl_matrix_get(m_tmp, i, j),
        mask_w = gsl_matrix_get(m_tmp, i - 1, j),
        mask_nw = gsl_matrix_get(m_tmp, i - 1, j + 1),
        mask_n = gsl_matrix_get(m_tmp, i, j + 1),
        mask_ne = gsl_matrix_get(m_tmp, i + 1, j + 1),
        mask_e = gsl_matrix_get(m_tmp, i + 1, j),
        mask_se = gsl_matrix_get(m_tmp, i + 1, j - 1),
        mask_s = gsl_matrix_get(m_tmp, i, j - 1),
        mask_sw = gsl_matrix_get(m_tmp, i - 1, j - 1);

      double thk = gsl_matrix_get(m_thickness, i, j),
        thk_w = gsl_matrix_get(m_thickness, i - 1, j),
        thk_nw = gsl_matrix_get(m_thickness, i - 1, j + 1),
        thk_n = gsl_matrix_get(m_thickness, i, j + 1),
        thk_ne = gsl_matrix_get(m_thickness, i + 1, j + 1),
        thk_e = gsl_matrix_get(m_thickness, i + 1, j),
        thk_se = gsl_matrix_get(m_thickness, i + 1, j - 1),
        thk_s = gsl_matrix_get(m_thickness, i, j - 1),
        thk_sw = gsl_matrix_get(m_thickness, i - 1, j - 1);

      // ice-free cell next to an icy cell
      if (thk < thk_eps &&
          (thk_w >= thk_eps || thk_nw >= thk_eps || thk_n >= thk_eps || thk_ne >= thk_eps ||
           thk_e >= thk_eps || thk_se >= thk_eps || thk_s >= thk_eps || thk_sw >= thk_eps)) {

        if (mask_w > 0)
          gsl_matrix_set(m_mask, i, j, mask_w);
        else if (mask_nw > 0)
          gsl_matrix_set(m_mask, i, j, mask_nw);
        else if (mask_n > 0)
          gsl_matrix_set(m_mask, i, j, mask_n);
        else if (mask_ne > 0)
          gsl_matrix_set(m_mask, i, j, mask_ne);
        else if (mask_e > 0)
          gsl_matrix_set(m_mask, i, j, mask_e);
        else if (mask_se > 0)
          gsl_matrix_set(m_mask, i, j, mask_se);
        else if (mask_s > 0)
          gsl_matrix_set(m_mask, i, j, mask_s);
        else if (mask_sw > 0)
          gsl_matrix_set(m_mask, i, j, mask_sw);
      } else {
        gsl_matrix_set(m_mask, i, j, m);
      }

    } // inner for loop
  } // outer for loop

  memcpy(new_mask, mask, Mx*My*sizeof(double));
}
