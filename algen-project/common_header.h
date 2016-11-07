#pragma once

#ifndef _COMMON_HEADER_H_
#define _COMMON_HEADER_H_
///////////////////////////////////////////////////////////////////////////////////////////////
// Include:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// DEFINED:
///////////////////////////////////////////////////////////////////////////////////////////////
#include "common_defined.h"
///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// DATA TYPES:
///////////////////////////////////////////////////////////////////////////////////////////////
#include "common_types.h"
///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// MACRO:
///////////////////////////////////////////////////////////////////////////////////////////////
#include "common_macro.h"
///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// GLOBALS:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////

#endif







//// checks for pragma once support
//// Sci¹gn¹³em z to z jednej stronki z schematami...
//// Chyba odpuszcze to sobie, bo raczej nie ma to sensu...
//#if (defined(__SUNPRO_C)&&(__SUNPRO_C>=0x5140))||(defined(__SUNPRO_CC)&&(__SUNPRO_CC>=0x5140))
////Oracle Developer Studio C/C++ (12.5 or later)
//#    define _pragma_once_support 1
//
////MSVC
//#elif defined(_MSC_VER)&&(_MSC_VER>=1020)
//#    define _pragma_once_support 1
//
//// clang:
//#elif defined(__clang__)
//#    define _pragma_once_support 1
//
//// comeau
//#elif defined(__COMO__)
//#    define _pragma_once_support 1
//
//// C++Builder (XE3 or greater)
//#elif defined(__CODEGEARC__)&&(__CODEGEARC__ >=650)
//#    define _pragma_once_support 1
//
//// Digital Mars
//#elif defined(__DMC__)
//#    define _pragma_once_support 1
//
//// GCC
//#elif defined(__GNUC__)&&((__GNUC__ >3)||(defined(__GNUC_MINOR__)&&(__GNUC__ ==3)&&(__GNUC_MINOR__ >=4)))
//#    define _pragma_once_support 1
//
//// HP aC++ (A.06.12)
//#elif defined(__HP_aCC)&&(__HP_aCC >=61200)
//#    define _pragma_once_support 1
//
//// IBM
//#elif defined(__xlC__)&&((__xlC__ >1301)||((__xlC__ ==1301)&&(__xlC_ver__ >0100)))
//#    define _pragma_once_support 1
//
//// intel
//#elif defined(__INTEL_COMPILER)||defined(__ICC)||defined(__ECC)||defined(__ICL)
//#    define _pragma_once_support 1
//
//// Pelles C
//#elif defined(__POCC__)			
//#    define _pragma_once_support 1
//
//// ARM compiler
//#elif defined(__CC_ARM)			
//#    define _pragma_once_support 1
//
//// IAR C/C++
//#elif defined(__IAR_SYSTEMS_ICC__)
//#    define _pragma_once_support 1
//
//// Portland Group C/C++
//#elif defined(__PGI)
//#    define _pragma_once_support 0
//
//#endif
//
//// if pragma once support then use it in addition to include guard
//#if defined(_pragma_once_support)
//#    pragma once
//#endif