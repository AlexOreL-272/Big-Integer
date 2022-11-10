#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-stack-address"
#pragma ide diagnostic ignored "misc-no-recursion"
#pragma once
#include <cmath>
#include <string>
#include <vector>

#include <iostream>
using namespace std;

class BigInt {
 private:
  // 12345678901 = {345678901, 12}
  vector<int64_t> big_int_; // it will store int32_t partitions of used number from

  enum signs {
    MINUS = 0,
    PLUS = 1
  };

  int sign_;
  const int kSizeOfCell_ = 9;
  const int64_t kMod_ = static_cast<int>(pow(10, 9));
 public:
  BigInt() {
    sign_ = PLUS;
  }

  BigInt(int64_t num) {
    sign_ = num < 0ll ? MINUS : PLUS;
    big_int_.push_back(llabs(num % kMod_)); // to store only absolute value of num
    if (num >= kMod_ * kMod_ || num <= -(kMod_ * kMod_)) {
      big_int_.push_back(llabs(num / kMod_) % kMod_);
      big_int_.push_back(llabs(num / kMod_ / kMod_));
    } else if (num >= kMod_ || num <= -kMod_) {
      big_int_.push_back(llabs(num / kMod_));
    }
  }

  BigInt(std::string string_to_init) {
    sign_ = string_to_init[0] == '-' ? MINUS : PLUS;
    int number_of_partitions = static_cast<int>(string_to_init.size() / kSizeOfCell_) - 1;
    int remaining_chars = static_cast<int>(string_to_init.size() % kSizeOfCell_);

    // extracting substrings of size kSizeOfCell_ from string
    for (; number_of_partitions >= 0; number_of_partitions--) {
      const int kPartition =
          std::stoi(string_to_init.substr(number_of_partitions * kSizeOfCell_ + remaining_chars, kSizeOfCell_));
      big_int_.push_back(kPartition);
    }

    // adding remaining number
    if (remaining_chars) {
      const int kIncludeFirstElem = sign_ == MINUS ? 1 : 0;
      big_int_.push_back(std::stoi(string_to_init.substr(kIncludeFirstElem, remaining_chars - kIncludeFirstElem)));
    }

    // fixing "-0" :)
    if (big_int_.size() == 1 && big_int_[0] == 0) {
      sign_ = PLUS;
    }
  }

  BigInt(const BigInt &obj_to_copy) {
    sign_ = obj_to_copy.sign_;
    std::copy(obj_to_copy.big_int_.begin(), obj_to_copy.big_int_.end(), std::back_inserter(big_int_));
  }

  void RemoveLeadingZeroes() {
    while (big_int_.size() > 1 && big_int_.back() == 0) {
      big_int_.erase(big_int_.end() - 1);
    }
    if (big_int_.size() == 1 && big_int_[0] == 0) {
      sign_ = PLUS;
    }
  }

  friend BigInt &operator-(BigInt &obj_to_reverse) {
    if (obj_to_reverse != 0) {
      obj_to_reverse.sign_ = obj_to_reverse.sign_ == PLUS ? MINUS : PLUS;
    }
    return obj_to_reverse;
  }

  BigInt &operator=(const BigInt &obj_to_assign) {
    sign_ = obj_to_assign.sign_;
    big_int_.clear();
    std::copy(obj_to_assign.big_int_.begin(), obj_to_assign.big_int_.end(), std::back_inserter(big_int_));
    return *this;
  }

  BigInt operator+(BigInt &obj_to_add) {
    BigInt result;

    if (sign_ != obj_to_add.sign_) {
      if (obj_to_add.sign_ == MINUS) {
        // l + (-r) = l - r
        obj_to_add.sign_ = PLUS;
        result = *this - obj_to_add;
        obj_to_add.sign_ = MINUS;
      } else {
        // (-l) + r = r - l
        sign_ = PLUS;
        result = obj_to_add - *this;
        sign_ = MINUS;
      }
      return result;
    }

    result.sign_ = sign_;
    int64_t remainder = 0;

    // adding common parts of objects
    for (size_t i = 0; i < min(big_int_.size(), obj_to_add.big_int_.size()); i++) {
      const int64_t kAddValue = big_int_[i] + obj_to_add.big_int_[i] + remainder;
      result.big_int_.push_back(kAddValue % kMod_);
      remainder = kAddValue / kMod_;
    }

    // adding additional parts of objects
    if (big_int_.size() != obj_to_add.big_int_.size()) {
      if (big_int_.size() < obj_to_add.big_int_.size()) {
        for (size_t i = big_int_.size(); i < obj_to_add.big_int_.size(); i++) {
          const int64_t kAddValue = obj_to_add.big_int_[i] + remainder;
          result.big_int_.push_back(kAddValue % kMod_);
          remainder = kAddValue / kMod_;
        }
      } else {
        for (size_t i = obj_to_add.big_int_.size(); i < big_int_.size(); i++) {
          const int64_t kAddValue = big_int_[i] + remainder;
          result.big_int_.push_back(kAddValue % kMod_);
          remainder = kAddValue / kMod_;
        }
      }
    }

    if (remainder) {
      result.big_int_.push_back(remainder);
    }

    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt operator-(BigInt &obj_to_sub) {
    BigInt result;

    if (sign_ != obj_to_sub.sign_) {
      if (obj_to_sub.sign_ == MINUS) {
        // +l - (-r) = l + r
        obj_to_sub.sign_ = PLUS;
        result = *this + obj_to_sub;
        obj_to_sub.sign_ = MINUS;
      } else {
        // (-l) - r = -(l + r)
        sign_ = PLUS;
        result = obj_to_sub + *this;
        result.sign_ = MINUS;
        sign_ = MINUS;
      }

      result.RemoveLeadingZeroes();
      return result;
    }

    if (sign_ == MINUS) {
      // sign_ == obj_to_sub.sign_ == MINUS
      obj_to_sub.sign_ = PLUS;
      result = obj_to_sub + *this;
      obj_to_sub.sign_ = MINUS;

      result.RemoveLeadingZeroes();
      return result;
    }

    // now sign_ == obj_to_sub.sign_ == PLUS
    if (big_int_.size() < obj_to_sub.big_int_.size()) { // I am doing that to substract a smaller number from bigger
      result = obj_to_sub - *this;
      result.sign_ = MINUS;

      result.RemoveLeadingZeroes();
      return result;
    }

    // if sizes are equal to one we don't need to use additional decimal unit
    if (big_int_.size() == obj_to_sub.big_int_.size() && big_int_.size() == 1) {
      result.big_int_.push_back(llabs(big_int_[0] - obj_to_sub.big_int_[0]));
      result.sign_ = *this < obj_to_sub ? MINUS : PLUS;
      return result;
    }

    int64_t remainder = 0;

    // substracting common parts of numbers
    for (size_t i = 0; i < obj_to_sub.big_int_.size(); i++) {
      const int64_t kSubValue = big_int_[i] - obj_to_sub.big_int_[i] + remainder;
      remainder = kSubValue < 0 ? -1 : 0;

      // need to determine whether to use additional decimal unit as in 25 - 7 we use 1 from 2 to substract 7 from 15
      result.big_int_.push_back((remainder == -1 ? (kSubValue + kMod_) : llabs(kSubValue)) % kMod_);
    }

    // substracting additional parts of object
    for (size_t i = obj_to_sub.big_int_.size(); i < big_int_.size(); i++) {
      const int64_t kSubValue = big_int_[i] + remainder;
      remainder = kSubValue < 0 ? -1 : 0;

      // need to determine whether to use additional decimal unit
      result.big_int_.push_back((remainder == -1 ? (kSubValue + kMod_) : llabs(kSubValue)) % kMod_);
    }

    result.sign_ = *this < obj_to_sub ? MINUS : PLUS;

    result.RemoveLeadingZeroes();
    return result;
  }

  BigInt MultSingle(int64_t num, size_t shift) {
    BigInt result;

    // assume shift to be >= 0
    int64_t remainder = 0;

    // multiplying a BigInt object by an int64_t number
    for (long long &item : big_int_) {
      int64_t kMulValue = item * num + remainder;
      remainder = kMulValue / kMod_;
      result.big_int_.push_back(kMulValue % kMod_);
    }

    if (remainder != 0) {
      result.big_int_.push_back(remainder);
    }

    // shifting received object by 'shift' positions
    if (shift != 0) {
      result.big_int_.resize(big_int_.size() + shift);
      for (size_t i = result.big_int_.size() - 1; i >= shift; i--) {
        result.big_int_[i] = result.big_int_[i - shift];
      }
      for (int i = 0; i < shift; i++) {
        result.big_int_[i] = 0;
      }
    }

    return result;
  }

  BigInt operator*(BigInt &obj_to_mul) {
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

  BigInt operator/(BigInt &obj_to_div) {
    BigInt result;

    return result;
  }

  BigInt operator%(BigInt &obj_to_mod) {
    BigInt result;

    return result;
  }

  BigInt &operator+=(BigInt &obj_to_add) {
    *this = *this + obj_to_add;
    return *this;
  }

  BigInt operator++(int) {
    BigInt copy = *this;
    ++*this;
    return copy;
  }

  BigInt &operator++() {
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

  BigInt &operator--() {
    BigInt one = 1;
    *this -= one;

    this->RemoveLeadingZeroes();
    return *this;
  }

  BigInt &operator-=(BigInt &obj_to_sub) {
    *this = *this - obj_to_sub;
    return *this;
  }

  BigInt &operator*=(BigInt &obj_to_mul) {
    *this = *this * obj_to_mul;
    return *this;
  }

  BigInt &operator/=(BigInt &obj_to_div) {
    *this = *this / obj_to_div;
    return *this;
  }

  BigInt &operator%=(BigInt &obj_to_mod) {
    *this = *this % obj_to_mod;
    return *this;
  }

  // -1 if this < compared object, 0 if equal, 1 else
  int CompareTo(const BigInt &obj_to_compare) const {
    if (sign_ != obj_to_compare.sign_) {
      return sign_ == MINUS ? -1 : 1;
    }

    // if sign is minus we need to compare numbers vice versa as if it was plus
    const int kMinusMultiplier = sign_ == MINUS ? -1 : 1;

    if (big_int_.size() != obj_to_compare.big_int_.size()) {
      return (big_int_.size() < obj_to_compare.big_int_.size() ? -1 : 1) * kMinusMultiplier;
    }

    // comparing the numbers with equal sign and body size
    for (int i = 0; i < big_int_.size(); i++) {
      if (big_int_[i] != obj_to_compare.big_int_[i]) {
        return (big_int_[i] < obj_to_compare.big_int_[i] ? -1 : 1) * kMinusMultiplier;
      }
    }

    return 0;
  }

  bool operator>(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == 1;
  }

  bool operator<(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == -1;
  }

  bool operator>=(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != -1;
  }

  bool operator<=(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != 1;
  }

  bool operator==(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag == 0;
  }

  bool operator!=(const BigInt &obj_to_compare) const {
    const int kCmpFlag = this->CompareTo(obj_to_compare);
    return kCmpFlag != 0;
  }

  friend std::ostream &operator<<(std::ostream &os, const BigInt &num_to_print) {
    if (num_to_print.sign_ == MINUS) {
      os << '-';
    }

    os << num_to_print.big_int_.back();

    // 'i' is not size_t 'cause it will be an infinite loop
    for (int i = static_cast<int>(num_to_print.big_int_.size()) - 2; i >= 0; i--) {
      const int64_t kNumToPrint = num_to_print.big_int_[i];
      os << std::string(num_to_print.kSizeOfCell_ - (int64_t) (log10(kNumToPrint)) - 1, '0') << kNumToPrint;
    }

    return os;
  }
};