#include <algorithm>
#include <iostream>
#include <string>
#include <iterator>

/**
 * C++ version std::string style "itoa":
 */

std::string itoa_base (const int value, const int base) {
  const char digitMap[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  std::string buf;

  /* Translating number to string with base 26: */
  int new_value = value;
  do {
	int index = new_value % base;
	buf += digitMap[index];
	new_value /= base;
  } while (new_value);
  std::string ret;
  std::copy(buf.rbegin(), buf.rend(), std::back_inserter(ret));
  return ret;
}

