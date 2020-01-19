#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <vector>

#include "obl_primitives.h"

namespace obl {

struct Foo {
  float key;
  float payload;

  static bool ogreater(const Foo &a, const Foo &b) { return a.key > b.key; }
};

struct Entry {
  float value;

  static std::string ToString(const Entry &a) {
    return std::to_string(a.value);
  }

  static bool ogreater(const Entry &a, const Entry &b) {
    return ObliviousGreater(a.value, b.value);
  }
};

struct Item {
  Entry entry;
  float rank;
  bool has_entry;

  inline bool operator<(const Item &rhs) const {
    return rank < rhs.rank ||
           (rank == rhs.rank && entry.value < rhs.entry.value);
  }

  static bool ogreater(const Item &a, const Item &b) {
    bool b0 = ogreater_impl(a, b);
    EXPECT_EQ(b0, b < a) << "a=" << ToString(a) << ", b=" << ToString(b);
    return b0;
  }

  static bool ogreater_impl(const Item &a, const Item &b) {
    bool b0 = ObliviousGreater(a.rank, b.rank);
    bool same_rank = ObliviousEqual(a.rank, b.rank);
    bool b1 = ObliviousChoose(same_rank,
                              (bool)Entry::ogreater(a.entry, b.entry), false);
    return ObliviousChoose(b0, true, b1);
  }

  static std::string ToString(const Item &a) {
    return Entry::ToString(a.entry) + ", rank=" + std::to_string(a.rank);
  }
};

bool operator>=(const Foo &a, const Foo &b) { return a.key >= b.key; }

std::ostream &operator<<(std::ostream &out, const Foo &a) {
  return out << "{ " << a.key << ", " << a.payload << " }";
}

std::ostream &operator<<(std::ostream &out, const std::vector<Foo> &vec) {
  for (const auto &foo : vec) {
    out << foo << ", ";
  }
  return out;
}

template <typename T>
void CheckSortedPod(const T *array, size_t n, bool ascending) {
  for (size_t idx = 1; idx < n; ++idx) {
    EXPECT_TRUE(array[idx] >= array[idx - 1]);
  }
}

float GetRandomKey() {
  std::random_device rd;
  return float(rd() % 100) / 10;
}

std::vector<Foo> GenerateRandomFoos(size_t n) {
  std::vector<Foo> ret;
  ret.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    Foo foo;
    foo.key = GetRandomKey();
    foo.payload = i;
    ret.push_back(foo);
  }
  return ret;
}

TEST(ObliviousSortPODTest, Works) {
  std::vector<Foo> array{{1, 2}, {2, 3}, {-1, 0}, {2, 3},
                         {1, 9}, {3, 4}, {2, 3}};
  ObliviousSortPOD(array.data(), 0, array.size(), true);
  CheckSortedPod(array.data(), array.size(), true);
  // std::cout << array << std::endl;
}

TEST(ObliviousSortPODTest, WorksRandom) {
  std::vector<Foo> array = GenerateRandomFoos(4096);
  ObliviousSortPOD(array.data(), 0, array.size(), true);
  CheckSortedPod(array.data(), array.size(), true);
  // std::cout << array << std::endl;
}

TEST(GT, Works) {
  {
    double a = 17.6;
    double b = 13.7;
    EXPECT_TRUE(ObliviousGreater(a, b));
  }

  {
    float a = 17.6;
    float b = 13.7;
    EXPECT_TRUE(ObliviousGreater(a, b));
  }

  // {
  //   Item a{Entry{std::numeric_limits<float>::max()}, 5.905882, false};
  //   Item b{Entry{std::numeric_limits<float>::max()}, 9.811765, false};
  //   EXPECT_FALSE(Item::ogreater(a, b));
  //   std::cout << "a=" << Item::ToString(a) << ", b=" << Item::ToString(b) <<
  //   std::endl;
  // }

  // {
  //   double a1 = -0.18;
  //   double b1 = -0.21;
  //   EXPECT_TRUE(ObliviousGreater(a1, b1));
  // }
}

} // namespace obl
