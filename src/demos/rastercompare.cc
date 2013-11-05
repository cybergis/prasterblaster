//
// Copyright 0000 <Nobody>
// @file
// @author David Matthew Mattli <dmattli@usgs.gov>
//
// @section LICENSE
//
// This software is in the public domain, furnished "as is", without
// technical support, and with no warranty, express or implied, as to
// its usefulness for any purpose.
//
// @section DESCRIPTION
//
//
//

#include <gdal_priv.h>
#include <gdal.h>

const double delta = 0.001;

int main(int argc, char *argv[]) {
  // Check command-line arguments
  if (argc < 3) {
    fprintf(stderr, "USAGE: rastercompare <control raster filename> "
            "<test raster filename>\n");
    return 1;
  }

  GDALAllRegister();

  GDALDataset *control =
      static_cast<GDALDataset*>(GDALOpen(argv[1], GA_ReadOnly));
  GDALDataset *test =
      static_cast<GDALDataset*>(GDALOpen(argv[2], GA_ReadOnly));

  if (control == NULL) {
    fprintf(stderr, "Error opening control raster!\n");
    return 1;
  }

  if (test == NULL) {
    fprintf(stderr, "Error opening test raster!\n");
    return 1;
  }

  // Basic checks
  if (control->GetRasterYSize() != test->GetRasterYSize()) {
    printf("Control raster has %d rows and test has %d\n",
           control->GetRasterYSize(),
           test->GetRasterYSize());
    return 1;
  }

  if (control->GetRasterXSize() != test->GetRasterXSize()) {
    printf("Control raster has %d columns and test has %d\n",
           control->GetRasterXSize(),
           test->GetRasterXSize());
    return 1;
  }

  // Check pixel values
  const int band_count = control->GetRasterCount();
  GDALRasterBand *band = control->GetRasterBand(1);
  int block_x_size, block_y_size;
  band->GetBlockSize(&block_x_size, &block_y_size);

  double *control_pixels, *test_pixels;

  control_pixels = new double[band_count
                               * control->GetRasterXSize() * sizeof(double)];
  test_pixels = new double[band_count
                            * control->GetRasterXSize() * sizeof(double)];

  // Loop over rows of blocks
  //   Loop over columns of blocks
  //     Loop over rows in block
  //      Loop over columns in block
  int bad_pixels = 0;
  const int y_size = control->GetRasterYSize();
  const int x_size = control->GetRasterXSize();
  for (int y = 0; y < y_size; ++y) {
    // Read a row
    control->RasterIO(GF_Read,
                      0,
                      y,
                      control->GetRasterXSize(),
                      1,
                      static_cast<void*>(control_pixels),
                      control->GetRasterXSize(),
                      1,
                      GDT_Float64,
                      band_count,
                      NULL,
                      0,
                      0,
                      0);
    test->RasterIO(GF_Read,
                   0,
                   y,
                   test->GetRasterXSize(),
                   1,
                   static_cast<void*>(test_pixels),
                   test->GetRasterXSize(),
                   1,
                   GDT_Float64,
                   band_count,
                   NULL,
                   0,
                   0,
                   0);
    if (memcmp(control_pixels, test_pixels, x_size * sizeof(double)) != 0) {
      bad_pixels++;
      continue;
      for (int x = 0; x < x_size; ++x) {
        //      printf("Control Value: %f, Test Value: %f, counter %d\n",
        //             control_pixels[x], test_pixels[x], counter++);
        if (fabs(control_pixels[x] - test_pixels[x]) > delta) {
          printf("Values at (%d, %d) are too different! %f vs %f\n",
                 x, y, control_pixels[x], test_pixels[x]);
          bad_pixels++;
        }
      }
    }
  }

  // Cleanup
  delete control_pixels;
  delete test_pixels;
  GDALClose(control);
  GDALClose(test);
  if (bad_pixels == 0) {
    return 0;
  } else {
    if (bad_pixels == 1) {
      printf("One pixel was different\n");
    } else {
      printf("%d pixels were different\n", bad_pixels);
    }
    return 1;
  }
}
