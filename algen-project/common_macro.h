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
	printf(_msg);																			  \
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






///////////////////////////////////////////////////////////////////////////////////////////////
// MATH:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// Cplx math:

#define CPLX_ADD(_z1, _z2, _out)															  \
	_out.x = _z1.x + _z2.x;																	  \
	_out.y = _z1.y + _z2.y;

#define CPLX_SUBTRACT(_z1, _z2, _out) {														  \
	_out.x = _z1.x - _z2.x;																	  \
	_out.y = _z1.y - _z2.y;																	  \
}

#define CPLX_MULTIPLY(_z1, _z2, _out) {														  \
	Data _tmp_out;																			  \
	_out.x = _z1.x * _z2.x - _z1.y * _z2.y;													  \
	_out.y = _z1.x * _z2.y + _z1.y * _z2.x;													  \
	_out.x = _tmp_out.x;																	  \
	_out.y = _tmp_out.y;																	  \
}

#define CPLX_DIVIDE(_z1, _z2, _out) {														  \
	Data _tmp_out;																			  \
	_tmp_out.x = (_z1.x * _z2.x + _z1.y * _z2.y) / (_z2.x * _z2.x + _z2.y * _z2.y);			  \
	_tmp_out.y = (_z1.y * _z2.x - _z1.x * _z2.y) / (_z2.x * _z2.x + _z2.y * _z2.y);			  \
	_out.x = _tmp_out.x;																	  \
	_out.y = _tmp_out.y;																	  \
}

#define CPLX_SCALE(_z, _scale, _out) {														  \
	Data _tmp_out;																			  \
	_tmp_out.x = _z.x * _scale;																  \
	_tmp_out.y = _z.y * _scale;																  \
	_out.x = _tmp_out.x;																	  \
	_out.y = _tmp_out.y;																	  \
}

#define CPLX_MODULE_2(_z)																	  \
	_z.x * _z.x + _z.y * _z.y

#define CPLX_MODULE(_z)																		  \
	sqrt( CPLX_MODULE_2(_z) )

#define CPLX_CONJUGATE(_z)																	  \
	_z.y = (-1) * _z.y;

///////////

#endif