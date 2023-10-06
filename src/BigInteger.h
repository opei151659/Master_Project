#ifndef BIGINTEGER_H
#define BIGINTEGER_H
#include <string>

using namespace std;

#define ullg unsigned long long

class BigInt {
private:
	string value;
public:
	BigInt() : value("0") { }
	BigInt(const string str) {
		value.append(str.size(), '0');
		for (int i = 0; i < str.size(); i++) {
			value[i] = str[str.size() - i - 1];
		}
	}
	BigInt(const ullg v) {
		string str = to_string(v);
		value.append(str.size(), '0');
		for (int i = 0; i < str.size(); i++) {
			value[i] = str[str.size() - i - 1];
		}
	}
	BigInt(const int v) {
		string str = to_string(v);
		value.append(str.size(), '0');
		for (int i = 0; i < str.size(); i++) {
			value[i] = str[str.size() - i - 1];
		}
	}
	BigInt(const char* a) {
		string str = a;
		value.append(str.size(), '0');
		for (int i = 0; i < str.size(); i++) {
			value[i] = str[str.size() - i - 1];
		}
	}
	BigInt(const BigInt& a) {
		value = a.value;
	}


	BigInt(BigInt&& other) {
		value = other.value;
	}

	void print_value() {
		for (int i = value.size(); i >= 0; i--) {
			std::cout << value[i];
		}
		//std::cout << std::endl;
		//std::cout << value << std::endl;
	}
	string get_value() {
		string str;
		str.append(value.size(), '0');
		for (int i = 0; i < str.size(); i++) {
			str[i] = value[str.size() - i - 1];
		}
		return str;
	}

	BigInt& operator= (const BigInt& a) {
		this->value = a.value;
		return *this;
	}
	BigInt& operator= (const int a) {
		BigInt tmp(a);
		*this = tmp;
		return *this;
	}
	BigInt& operator= (const ullg a) {
		BigInt tmp(a);
		*this = tmp;
		return *this;
	}
	BigInt& operator= (const string a) {
		BigInt tmp(a);
		*this = tmp;
		return *this;
	}
	BigInt& operator= (const char* a) {
		BigInt tmp(a);
		*this = tmp;
		return *this;
	}

	BigInt operator+ (const BigInt& a) {
		BigInt tmp = *this;
		tmp += a;
		return tmp;
	}
	BigInt operator+ (const ullg a) {
		BigInt tmp = *this;
		tmp += a;
		return tmp;
	}

	BigInt operator* (const BigInt& a) {
		BigInt tmp = *this;
		tmp *= a;
		return tmp;
	}
	BigInt operator* (const ullg a) {
		BigInt tmp(a);
		tmp *= *this;
		return tmp;
	}

	bool operator> (const BigInt& other) const {
		if (value.size() > other.value.size()) {
			return true;
		}
		else if (value.size() < other.value.size()) {
			return false;
		}
		else {
			for (int i = value.size(); i >= 0; i--) {
				if (value[i] > other.value[i]) return true;
				else if (value[i] < other.value[i]) return false;
			}
			return false;
		}
	}

	bool operator< (const BigInt& other) const {
		if (value.size() < other.value.size()) {
			return true;
		}
		else if (value.size() > other.value.size()) {
			return false;
		}
		else {
			for (int i = value.size(); i >= 0; i--) {
				if (value[i] < other.value[i]) return true;
				else if (value[i] > other.value[i]) return false;
			}
			return false;
		}
	}

	friend BigInt& operator*= (BigInt& a, const BigInt& b) {
		if (a.value == "0" || b.value == "0") {
			a.value = "0";
			return a;
		}
		BigInt tmp, r;
		int carry, v, add;
		for (int i = 0; i < b.value.size(); i++) {
			tmp.value.assign(a.value.size() + i, '0');
			for (int j = 0; j < a.value.size(); j++) {
				v = (a.value[j] - '0') * (b.value[i] - '0');
				carry = 0;
				for (int p = j + i; v > 0; p++) {
					if (p >= tmp.value.size()) {
						tmp.value.append("0");
					}
					add = v % 10 + (tmp.value[p] - '0') + carry;
					if (add > 10) {
						carry = add / 10;
						add %= 10;
					}
					else {
						carry = 0;
					}
					v /= 10;
					tmp.value[p] = add + '0';
				}
			}
			r += tmp;
		}
		a = r;
		/*int v, carry = 0;
		BigInt tmp, add;
		for (int i = 0; i < b.value.size(); i++) {
			for (int j = 0; j < a.value.size(); j++) {
				tmp += ((a.value[j] - '0') * (b.value[i] - '0') + carry) * pow(10, i + j);
			}
		}
		a = tmp;*/
		return a;
	}

	friend BigInt& operator+= (BigInt& a, const BigInt& b) {
		int len1 = a.value.size(), len2 = b.value.size();
		int carry = 0, v;
		if (len1 < len2)
			a.value.append(len2 - len1, '0');
		int len = a.value.size();
		for (int i = 0; i < len; i++) {
			if (i < len2) {
				v = (a.value[i] - '0') + (b.value[i] - '0') + carry;
				if (v > 9) {
					carry = 1;
					v -= 10;
				}
				else carry = 0;
			}
			else {
				v = (a.value[i] - '0') + carry;
				if (v > 9) {
					carry = 1;
					v -= 10;
				}
				else carry = 0;
			}
			a.value[i] = v + '0';
		}
		if (carry) {
			a.value.append("1");
		}
		return a;
	}
	friend BigInt& operator+= (BigInt& a, const ullg b) {
		a += BigInt(b);
		return a;
	}
	/*if use pow (int) -> (double), error. so this situation use ullg*/
	//friend BigInt& operator+= (BigInt& a, const int b) {
	//	a += BigInt(b);
	//	return a;
	//}
	friend bool operator== (const BigInt& a, const BigInt& b) {
		if (a.value.size() != b.value.size())
			return false;
		for (int i = 0; i < a.value.size(); i++) {
			if (a.value[i] != b.value[i])
				return false;
		}
		return true;
	}
	friend bool operator!= (const BigInt& a, const BigInt& b) {
		if (a.value.size() != b.value.size())
			return true;
		for (int i = 0; i < a.value.size(); i++) {
			if (a.value[i] != b.value[i])
				return true;
		}
		return false;
	}


};

#endif BIGINTEGER_H