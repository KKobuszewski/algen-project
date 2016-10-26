#pragma once

#ifndef _COMMON_MACROS_H_
#define _COMMON_MACROS_H_
///////////////////////////////////////////////////////////////////////////////////////////////
// INCLUDE:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// My:
#include "common_defined.h"
#include "common_types.h"

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// PRINT-ING:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// Clasic C:

#define	CPRINT(_msg) {																		  \
	fprintf(stdout, _msg);																	  \
}

///////////////////////////////////
// C++ 

//#define PRINT(_msg) {																		  \
//	std::cout _msg;																			  \
//}

//#define PRINT_ERROR {																		  \
//	PRINT(<< "File: " << __FILE__                 << std::endl );							  \
//	PRINT(<< "Line: " << __LINE__                 << std::endl );							  \
//	PRINT(<< "Func: " << __MY_FUNC_NAME__ << "()" << std::endl );							  \
//}

//#define PRINT_ERROR_C(_my_text) {															  \
//	PRINT(<< "File: " << __FILE__                 << std::endl );							  \
//	PRINT(<< "Line: " << __LINE__                 << std::endl );							  \
//	PRINT(<< "Func: " << __MY_FUNC_NAME__ << "()" << std::endl );							  \
//	PRINT(<< "Text: " << _my_text                 << std::endl );							  \
//}

///////////

#endif