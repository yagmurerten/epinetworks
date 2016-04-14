#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#ifdef _DEBUG
#  include <cassert>
#endif // _DEBUG

#ifdef _DEBUG
#  define DEBUG_CODE(debug_code) debug_code
#  define DEBUG_OR_RELEASE_CODE(debug_code, release_code) debug_code
#else
#  define DEBUG_CODE(debug_code)
#  define DEBUG_OR_RELEASE_CODE(debug_code, release_code) release_code
#endif // _DEBUG

#define DEBUG_ASSERT(predicate) DEBUG_CODE(assert(predicate))

#define DEBUG_ERROR DEBUG_ASSERT(false)

#endif // DEBUG_H_INCLUDED
