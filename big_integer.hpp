#pragma once
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class BigInt {
 private:
  // 12345678901 = {345678901, 12}
  std::vector<int64_t> big_int_;
  // it will store int32_t partitions of used number from

  enum Signs { MINUS = 0, PLUS = 1 };

  int sign_;
  const size_t kSizeOfCell = 9;
  const int64_t kMod = 1e9;

 public:
  BigInt() { sign_ = PLUS; }

  BigInt(int64_t num) {
    sign_ = num < 0LL ? MINUS : PLUS;
    big_int_.push_back(
        llabs(num % kMod));  // to store only absolute value of num
    if (num >= kMod * kMod || num <= -(kMod * kMod)) {
      big_int_.push_back(llabs(num / kMod) % kMod);
      big_int_.push_back(llabs(num / kMod / kMod));
    } else if (num >= kMod || num <= -kMod) {
      big_int_.push_back(llabs(num / kMod));
    }
  }

  BigInt(std::string string_to_init) {
    sign_ = string_to_init[0] == '-' ? MINUS : PLUS;
    int number_of_partitions =
        static_cast<int>(string_to_init.size() / kSizeOfCell) - 1;
    int remaining_chars = static_cast<int>(string_to_init.size() % kSizeOfCell);

    // extracting substrings of size kSizeOfCell from string
    for (; number_of_partitions >= 0; number_of_partitions--) {
      const int kPartition = std::stoi(string_to_init.substr(
          number_of_partitions * kSizeOfCell + remaining_chars, kSizeOfCell));
      big_int_.push_back(kPartition);
    }

    // adding remaining number
    if (remaining_chars != 0) {
      const int kIncludeFirstElem = sign_ == MINUS ? 1 : 0;
      big_int_.push_back(std::stoi(string_to_init.substr(
          kIncludeFirstElem, remaining_chars - kIncludeFirstElem)));
    }

    // fixing "-0" :)
    if (big_int_.size() == 1 && big_int_[0] == 0) {
      sign_ = PLUS;
    }
  }

  BigInt(const BigInt& obj_to_copy) {
    sign_ = obj_to_copy.sign_;
    big_int_.clear();
    std::copy(obj_to_copy.big_int_.begin(), obj_to_copy.big_int_.end(),
              std::back_inserter(big_int_));
  }

  void RemoveLeadingZeroes() {
    while (big_int_.size() > 1 && big_int_.back() == 0) {
      big_int_.erase(big_int_.end() - 1);
    }
    if (big_int_.size() == 1 && big_int_[0] == 0) {
      sign_ = PLUS;
    }
  }

  friend BigInt& operator-(BigInt& obj_to_reverse) {
    if (obj_to_reverse != 0) {
      obj_to_reverse.sign_ = obj_to_reverse.sign_ == PLUS ? MINUS : PLUS;
    }
    return obj_to_reverse;
  }

  BigInt& operator=(const BigInt& obj_to_assign) {
    sign_ = obj_to_assign.sign_;
    big_int_.clear();

    std::copy(obj_to_assign.big_int_.begin(), obj_to_assign.big_int_.end(),
              std::back_inserter(big_int_));
    return *this;
  }

  friend BigInt AddDifferentSigns(BigInt& left, BigInt& right) {
    BigInt result;

    if (right.sign_ == MINUS) {
      // l + (-r) = l - r
      right.sign_ = PLUS;
      result = left - right;
      right.sign_ = MINUS;
    } else {
      // (-l) + r = r - l
      left.sign_ = PLUS;
      result = right - left;
      left.sign_ = MINUS;
    }

    return result;
  }

  BigInt operator+(BigInt& obj_to_add) {
    BigInt result;
    result.sign_ = sign_;
    int64_t remainder = 0;
    if (sign_ != obj_to_add.sign_) {
      return AddDifferentSigns(*this, obj_to_add);
    }

    // I want to add a smaller number to a bigger
    if (*this < obj_to_add) {
      return obj_to_add + *this;
    }

    // adding common parts of objects
    for (size_t i = 0; i < obj_to_add.big_int_.size(); i++) {
      const int64_t kAddValue =
          big_int_[i] + obj_to_add.big_int_[i] + remainder;
      result.big_int_.push_back(kAddValue % kMod);
      remainder = kAddValue / kMod;
    }

    // adding additional parts of objects
    for (size_t i = obj_to_add.big_int_.size(); i < big_int_.size(); i++) {
      const int64_t kAddValue = big_int_[i] + remainder;
      result.big_int_.push_back(kAddValue % kMod);
      remainder = kAddValue / kMod;
    }

    if (remainder != 0) {
      result.big_int_.push_back(remainder);
    }

    result.RemoveLeadingZeroes();
    return result;
  }

  friend BigInt SubOneMinus(BigInt& left, BigInt& right) {
    BigInt result;

    if (left.sign_ != right.sign_) {
      if (right.sign_ == MINUS) {
        // +l - (-r) = l + r
        right.sign_ = PLUS;
        result = left + right;
        right.sign_ = MINUS;
      } else {
        // (-l) - r = -(l + r)
        left.sign_ = PLUS;
        result = right + left;
        result.sign_ = MINUS;
        left.sign_ = MINUS;
      }
    } else if (left.sign_ == MINUS) {
      // sign_ == obj_to_sub.sign_ == MINUS
      right.sign_ = PLUS;
      result = right + left;
      right.sign_ = MINUS;
    }

    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt operator-(BigInt& obj_to_sub) {
    BigInt result;
    if (sign_ != obj_to_sub.sign_ || sign_ == MINUS) {
      return SubOneMinus(*this, obj_to_sub);
    }

    if (*this < obj_to_sub) {  // Want to substract a smaller number from bigger
      BigInt temp = obj_to_sub - *this;
      result = -temp;
      result.RemoveLeadingZeroes();
      return result;
    }

    // if sizes are equal to one we don't need to use additional decimal unit
    if (big_int_.size() == obj_to_sub.big_int_.size() && big_int_.size() == 1) {
      result.big_int_.push_back(llabs(big_int_[0] - obj_to_sub.big_int_[0]));
      result.sign_ = *this < obj_to_sub ? MINUS : PLUS;
      return result;
    }

    for (size_t i = 0, tmp = obj_to_sub.big_int_.size();
         i < big_int_.size() - tmp; i++) {  // to make sizes equal
      obj_to_sub.big_int_.push_back(0LL);
    }

    int64_t rem = 0;  // stands for remainder
    for (size_t i = 0; i < big_int_.size(); i++) {
      const int64_t kSub = big_int_[i] - obj_to_sub.big_int_[i] + rem;
      rem = kSub < 0 ? -1 : 0;
      result.big_int_.push_back((rem < 0 ? (kSub + kMod) : llabs(kSub)) % kMod);
    }

    result.sign_ = *this < obj_to_sub ? MINUS : PLUS;
    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt MultSingle(int64_t num, size_t shift = 0) const {
    BigInt result;

    // assume shift to be >= 0
    int64_t remainder = 0;

    // multiplying a BigInt object by an int64_t number
    for (const auto& item : big_int_) {
      const int64_t kMulValue = item * num + remainder;
      remainder = kMulValue / kMod;
      result.big_int_.push_back(kMulValue % kMod);
    }

    if (remainder != 0) {
      result.big_int_.push_back(remainder);
    }

    // shifting received object by 'shift' positions
    if (shift != 0) {
      for (size_t i = 0; i < shift; i++) {
        result.Unshift(0);
      }
    }

    return result;
  }

  BigInt operator*(BigInt& obj_to_mul) const {
    BigInt result;

    // as column multiplying
    for (size_t i = 0; i < obj_to_mul.big_int_.size(); i++) {
      BigInt layer = MultSingle(obj_to_mul.big_int_[i], i);
      result += layer;
    }

    result.sign_ = sign_ == obj_to_mul.sign_ ? PLUS : MINUS;

    result.RemoveLeadingZeroes();
    return result;
  }

  // works as unshift in JS
  BigInt& Unshift(int64_t num_to_insert) {
    big_int_.push_back(num_to_insert);
    std::rotate(big_int_.begin(), big_int_.end() - 1, big_int_.end());
    return *this;
  }

  friend int64_t BinSearch(const BigInt& obj_to_div, BigInt& number_to_sub) {
    int64_t required_num = 0;
    int64_t left = 0;
    int64_t right = obj_to_div.kMod;

    while (left <= right) {
      int64_t mid = (left + right) / 2;
      BigInt multed = obj_to_div.MultSingle(mid);
      number_to_sub.sign_ = PLUS;
      if (multed <= number_to_sub) {
        required_num = mid;
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }

    return required_num;
  }

  BigInt operator/(BigInt& obj_to_div) {
    BigInt result;
    result.big_int_.resize(big_int_.size());

    if (obj_to_div == 0) {
      std::cerr << "Who are you such sophisticated in science?\n";
      exit(1);
    }

    auto temp_sign = obj_to_div.sign_;
    obj_to_div.sign_ = PLUS;
    BigInt number_to_substract;
    for (int i = (int)big_int_.size() - 1; i >= 0; i--) {
      number_to_substract.Unshift(big_int_[i]);
      number_to_substract.RemoveLeadingZeroes();

      // bin search of needed multiplyer
      int64_t required_num = BinSearch(obj_to_div, number_to_substract);

      result.big_int_[i] = required_num;
      BigInt multed = obj_to_div.MultSingle(required_num);
      multed.RemoveLeadingZeroes();
      number_to_substract -= multed;
    }

    obj_to_div.sign_ = temp_sign;
    result.sign_ = sign_ == obj_to_div.sign_ ? PLUS : MINUS;
    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt operator%(BigInt& obj_to_mod) {
    BigInt whole_part = (*this / obj_to_mod) * obj_to_mod;
    BigInt result = *this - whole_part;
    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt& operator+=(BigInt& obj_to_add) {
    *this = *this + obj_to_add;
    return *this;
  }

  BigInt operator++(int) {
    BigInt copy = *this;
    ++*this;
    return copy;
  }

  BigInt& operator++() {
    BigInt one = 1;
    *this += one;

    this->RemoveLeadingZeroes();
    return *this;
  }

  BigInt operator--(int) {
    BigInt copy = *this;
    --*this;
    return copy;
  }

  BigInt& operator--() {
    BigInt one = 1;
    *this -= one;

    this->RemoveLeadingZeroes();
    return *this;
  }

  BigInt& operator-=(BigInt& obj_to_sub) {
    *this = *this - obj_to_sub;
    return *this;
  }

  BigInt& operator*=(BigInt& obj_to_mul) {
    *this = *this * obj_to_mul;
    return *this;
  }

  BigInt& operator/=(BigInt& obj_to_div) {
    *this = *this / obj_to_div;
    return *this;
  }

  BigInt& operator%=(BigInt& obj_to_mod) {
    *this = *this % obj_to_mod;
    return *this;
  }

  // -1 if this < compared object, 0 if equal, 1 else
  int CompareTo(const BigInt& obj_to_compare) const {
    if (sign_ != obj_to_compare.sign_) {
      return sign_ == MINUS ? -1 : 1;
    }

    // if sign is minus we need to compare numbers vice versa as if it was plus
    const int kMinusMultiplier = sign_ == MINUS ? -1 : 1;

    if (big_int_.size() != obj_to_compare.big_int_.size()) {
      return (big_int_.size() < obj_to_compare.big_int_.size() ? -1 : 1) *
             kMinusMultiplier;
    }

    // comparing the numbers with equal sign and body size
    for (int i = (int)big_int_.size() - 1; i >= 0; i--) {
      if (big_int_[i] != obj_to_compare.big_int_[i]) {
        return (big_int_[i] < obj_to_compare.big_int_[i] ? -1 : 1) *
               kMinusMultiplier;
      }
    }

    return 0;
  }

  bool operator>(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == 1;
  }

  bool operator<(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == -1;
  }

  bool operator>=(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != -1;
  }

  bool operator<=(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != 1;
  }

  bool operator==(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == 0;
  }

  bool operator!=(const BigInt& obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != 0;
  }

  friend std::istream& operator>>(std::istream& is, BigInt& num_to_scan) {
    std::string string_to_scan;
    is >> string_to_scan;
    num_to_scan = BigInt(string_to_scan);
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const BigInt& num_to_print) {
    if (num_to_print.sign_ == MINUS) {
      os << '-';
    }

    if (num_to_print.big_int_.empty()) {
      return os << 0;
    }

    os << num_to_print.big_int_.back();

    // 'i' is not size_t 'cause it will be an infinite loop
    for (int i = (int)(num_to_print.big_int_.size()) - 2; i >= 0; i--) {
      const int64_t kNumToPrint = num_to_print.big_int_[i];
      os << std::setw(num_to_print.kSizeOfCell) << std::setfill('0')
         << kNumToPrint;
    }

    return os;
  }
};