//
// sharedptr.hh for pekwm
// Copyright © 2012 Andreas Schlick <ioerror{@}lavabit{.}com>
//
// This program is licensed under the GNU GPL.
// See the LICENSE file for more information.
//

#ifndef _SHAREDPTR_HH_
#define _SHAREDPTR_HH_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STD_SHARED_PTR
#include <memory>
using std::shared_ptr;
using std::dynamic_pointer_cast;
#elif HAVE_TR1_SHARED_PTR
#include <tr1/memory>
using std::tr1::shared_ptr;
using std::tr1::dynamic_pointer_cast;
#elif HAVE_BOOST_SHARED_PTR
#include <boost/memory>
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
#endif

class WinLayouter;
typedef shared_ptr<WinLayouter> WLayouter;

#endif // _SHAREDPTR_HH_
