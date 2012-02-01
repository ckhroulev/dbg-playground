#include "drainagebasin.hh"

void init_mask(int Mx, int My,
               double *elevation,
               double *thickness,
               double *mask) {

  gsl_matrix_view mask_view = gsl_matrix_view_array(mask, Mx, My);
  gsl_matrix * m_mask = &mask_view.matrix;

  gsl_matrix_view elevation_view = gsl_matrix_view_array(elevation, Mx, My);
  gsl_matrix * m_elevation = &elevation_view.matrix;

  gsl_matrix_view thickness_view = gsl_matrix_view_array(thickness, Mx, My);
  gsl_matrix * m_thickness = &thickness_view.matrix;

  int marker = 1;

  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {
      if (i == 0 || i == Mx - 1 ||
          j == 0 || j == My - 1) {
        gsl_matrix_set(m_mask, i, j, -1);
        continue;
      }

      double h = gsl_matrix_get(m_elevation, i, j),
        h_w = gsl_matrix_get(m_elevation, i - 1, j),
        h_nw = gsl_matrix_get(m_elevation, i - 1, j + 1),
        h_n = gsl_matrix_get(m_elevation, i, j + 1),
        h_ne = gsl_matrix_get(m_elevation, i + 1, j + 1),
        h_e = gsl_matrix_get(m_elevation, i + 1, j),
        h_se = gsl_matrix_get(m_elevation, i + 1, j - 1),
        h_s = gsl_matrix_get(m_elevation, i, j - 1),
        h_sw = gsl_matrix_get(m_elevation, i - 1, j - 1);

      double thk = gsl_matrix_get(m_thickness, i, j),
        thk_w = gsl_matrix_get(m_thickness, i - 1, j),
        thk_nw = gsl_matrix_get(m_thickness, i - 1, j + 1),
        thk_n = gsl_matrix_get(m_thickness, i, j + 1),
        thk_ne = gsl_matrix_get(m_thickness, i + 1, j + 1),
        thk_e = gsl_matrix_get(m_thickness, i + 1, j),
        thk_se = gsl_matrix_get(m_thickness, i + 1, j - 1),
        thk_s = gsl_matrix_get(m_thickness, i, j - 1),
        thk_sw = gsl_matrix_get(m_thickness, i - 1, j - 1);

      // icy cell next to an ice-free cell
      if (thk > 1 &&
          (thk_w <= 1 || thk_nw <= 1 || thk_n <= 1 || thk_ne <= 1 ||
           thk_e <= 1 || thk_se <= 1 || thk_s <= 1 || thk_sw <= 1)) {

        // there is an ice-free cell next to the current cell with the same or
        // lower elevation
        if (h_w <= h || h_nw <= h || h_n <= h || h_ne <= h ||
            h_e <= h || h_se <= h || h_s <= h || h_sw <= h) {
          gsl_matrix_set(m_mask, i, j, ++marker);
        } else {
          gsl_matrix_set(m_mask, i, j, -2);
        }

      } // icy next to ice-free

    } // inner for loop
  } // outer for loop

}
