

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

namespace librasterblaster {
enum PRB_ERROR {
  PRB_NOERROR,
  PRB_IOERROR,
  PRB_BADARG,
  PRB_PROJERROR,
};

struct Area {
  Area() {}
  Area(double ulx, double uly, double lrx, double lry) : ul(ulx, uly, UNDEF), lr(lrx, lry, UNDEF) {}
		
  Coordinate ul;
  Coordinate lr;
  ProjUnit units;
	
};
}

#endif // SRC_UTILS_H_
