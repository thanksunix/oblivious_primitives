#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace obl {

constexpr size_t CACHE_LINE_SIZE = 64;

template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
uint8_t ObliviousGreater(T x, T y);

template <typename T, typename std::enable_if<std::is_floating_point<T>::value,
                                              int>::type = 0>
uint8_t ObliviousGreater(T x, T y);

template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
uint8_t ObliviousGreaterOrEqual(T x, T y);

uint8_t ObliviousGreaterOrEqual(double x, double y);

template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
uint8_t ObliviousLess(T x, T y);

uint8_t ObliviousLess(double x, double y);

template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
uint8_t ObliviousLessOrEqual(T x, T y);

uint8_t ObliviousLessOrEqual(double x, double y);

template <typename T,
          typename std::enable_if<std::is_scalar<T>::value, int>::type = 0>
uint8_t ObliviousEqual(T x, T y);

template <typename T>
T ObliviousChoose(bool pred, const T &t_val, const T &f_val);

template <typename T,
          typename std::enable_if<std::is_trivially_destructible<T>::value,
                                  int>::type = 0>
inline void ObliviousAssign(bool pred, const T &t_val, const T &f_val, T *out);

template <typename T>
inline void ObliviousMerge(T *arr, uint32_t low, uint32_t len, bool ascending);

template <typename T>
inline void ObliviousSort(T *arr, uint32_t low, uint32_t len, bool ascending);

template <typename T> inline T ObliviousArrayAccess(T *arr, size_t i, size_t n);

template <typename T>
void ObliviousArrayAssign(T *arr, size_t i, size_t n, const T &val);

//------------------------------------------------------
// Implementation.

// Helper functions

//------------------------------------------------------
// Oblivious primitives

// Return x > y
template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type>
inline uint8_t ObliviousGreater(T x, T y) {
  uint8_t result;
  __asm__ volatile("cmp %2, %1;"
                   "setg %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

template <typename T,
          typename std::enable_if<std::is_floating_point<T>::value, int>::type>
inline uint8_t ObliviousGreater(T _x, T _y) {
  double x = _x;
  double y = _y;
  uint8_t result;
  __asm__ volatile("comisd %%xmm1, %%xmm0;"
                   "seta %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

// Return x >= y
template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type>
inline uint8_t ObliviousGreaterOrEqual(T x, T y) {
  uint8_t result;
  __asm__ volatile("cmp %2, %1;"
                   "setge %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

inline uint8_t ObliviousGreaterOrEqual(double x, double y) {
  uint8_t result;
  __asm__ volatile("comisd %%xmm1, %%xmm0;"
                   "setae %0;"
                   : "=r"(result)
                   :
                   : "cc");
  return result;
}

// Return x == y
template <typename T,
          typename std::enable_if<std::is_scalar<T>::value, int>::type>
inline uint8_t ObliviousEqual(T x, T y) {
  uint8_t result;
  __asm__ volatile("cmp %2, %1;"
                   "sete %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

// Return x < y
template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type>
inline uint8_t ObliviousLess(T x, T y) {
  uint8_t result;
  __asm__ volatile("cmp %2, %1;"
                   "setl %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

inline uint8_t ObliviousLess(double x, double y) {
  uint8_t result;
  __asm__ volatile("comisd %%xmm1, %%xmm0;"
                   "setb %0;"
                   : "=r"(result)
                   :
                   : "cc");
  return result;
}

// Return x <= y
template <typename T,
          typename std::enable_if<std::is_integral<T>::value, int>::type>
inline uint8_t ObliviousLessOrEqual(T x, T y) {
  uint8_t result;
  __asm__ volatile("cmp %2, %1;"
                   "setle %0;"
                   : "=r"(result)
                   : "r"(x), "r"(y)
                   : "cc");
  return result;
}

inline uint8_t ObliviousLessOrEqual(double x, double y) {
  uint8_t result;
  __asm__ volatile("comisd %%xmm1, %%xmm0;"
                   "setbe %0;"
                   : "=r"(result)
                   :
                   : "cc");
  return result;
}

// Fill `out` with (pred ? t_val : f_val). Supports 16, 32, and 64 bit types
template <typename T,
          typename std::enable_if<std::is_scalar<T>::value, int>::type = 0>
inline void ObliviousAssignHelper(bool pred, const T &t_val, const T &f_val,
                                  T *out) {
  T result;
  __asm__ volatile("mov %2, %0;"
                   "test %1, %1;"
                   "cmovz %3, %0;"
                   : "=&r"(result)
                   : "r"(pred), "r"(t_val), "r"(f_val), "m"(out)
                   : "cc");
  *out = result;
}

// Iteratively apply `ObliviousAssign` to fill generic types of size > 1 byte
template <typename T, typename std::enable_if<
                          std::is_trivially_destructible<T>::value, int>::type>
inline void ObliviousAssign(bool pred, const T &t_val, const T &f_val, T *out) {
  size_t bytes = sizeof(T);
  char *res = (char *)out;
  char *t = (char *)&t_val;
  char *f = (char *)&f_val;

  // Obliviously assign 8 bytes at a time
  size_t num_8_iter = bytes / 8;
#pragma omp simd
  for (int i = 0; i < num_8_iter; i++) {
    ObliviousAssignHelper(pred, *((uint64_t *)t), *((uint64_t *)f),
                          (uint64_t *)res);
    res += 8;
    t += 8;
    f += 8;
  }

  // Obliviously assign 4 bytes
  if ((bytes % 8) / 4) {
    ObliviousAssignHelper(pred, *((uint32_t *)t), *((uint32_t *)f),
                          (uint32_t *)res);
    res += 4;
    t += 4;
    f += 4;
  }

  // Obliviously assign 2 bytes
  if ((bytes % 4) / 2) {
    ObliviousAssignHelper(pred, *((uint16_t *)t), *((uint16_t *)f),
                          (uint16_t *)res);
    res += 2;
    t += 2;
    f += 2;
  }

  if ((bytes % 2)) {
    ObliviousAssignHelper(pred, *((uint16_t *)t), *((uint16_t *)f),
                          (uint16_t *)res);
  }
}

// Return (pred ? t_val : f_val). Supports types of size > 1 byte
template <typename T>
inline T ObliviousChoose(bool pred, const T &t_val, const T &f_val) {
  T result;
  ObliviousAssign(pred, t_val, f_val, &result);
  return result;
}

// Return arr[i]
template <typename T> inline T ObliviousArrayAccess(T *arr, int i, size_t n) {
  T result = arr[0];
  int step = sizeof(T) < CACHE_LINE_SIZE ? CACHE_LINE_SIZE / sizeof(T) : 1;
  for (int j = 0; j < n; j += step) {
    bool cond = ObliviousEqual(j / step, i / step);
    int pos = ObliviousChoose(cond, i, j);
    result = ObliviousChoose(cond, arr[pos], result);
  }
  return result;
}

// Set arr[i] = val
template <typename T>
inline void ObliviousArrayAssign(T *arr, int i, size_t n, T val) {
  int step = sizeof(T) < CACHE_LINE_SIZE ? CACHE_LINE_SIZE / sizeof(T) : 1;
  for (int j = 0; j < n; j += step) {
    bool cond = ObliviousEqual(j / step, i / step);
    int pos = ObliviousChoose(cond, i, j);
    arr[pos] = ObliviousChoose(cond, val, arr[pos]);
  }
}

namespace detail {

inline uint32_t greatest_power_of_two_less_than(uint32_t n) {
  uint32_t k = 1;
  while (k < n)
    k = k << 1;
  return k >> 1;
}

inline uint32_t log2_ceil(uint32_t n) {
  uint32_t k = 0;
  uint32_t _n = n;
  while (n > 1) {
    k++;
    n /= 2;
  }
  if ((1 << k) < _n)
    k++;
  return k;
}

/***************************************************************************************
 * Oblivious bitonic sort (imperative version)
 **************************************************************************************/

// Imperative implementation of bitonic merge network
template <typename T>
inline void imperative_o_merge(T *arr, uint32_t low, uint32_t len,
                               bool ascending) {
  uint32_t i, j, k;
  uint32_t l = log2_ceil(len);
  uint32_t n = 1 << l;
  for (i = 0; i < l; i++) {
    for (j = 0; j<n; j += n>> i) {
      for (k = 0; k < (n >> i) / 2; k++) {
        uint32_t i1 = low + k + j;
        uint32_t i2 = i1 + (n >> i) / 2;
        if (i2 >= low + len)
          break;
        bool pred = ObliviousGreater(arr[i1], arr[i2]);
        pred = ObliviousEqual(pred, ascending);
        // These array accesses are oblivious because the indices are
        // deterministic
        T tmp = arr[i1];
        arr[i1] = ObliviousChoose(pred, arr[i2], arr[i1]);
        arr[i2] = ObliviousChoose(pred, tmp, arr[i2]);
      }
    }
  }
}

// Imperative implementation of bitonic merge network for POD type
// Assumes T implements T::ogreater
template <typename T>
inline void imperative_o_merge_pod(T *arr, uint32_t low, uint32_t len,
                                   bool ascending) {
  uint32_t i, j, k;
  uint32_t l = log2_ceil(len);
  uint32_t n = 1 << l;
  for (i = 0; i < l; i++) {
    for (j = 0; j<n; j += n>> i) {
      for (k = 0; k < (n >> i) / 2; k++) {
        uint32_t i1 = low + k + j;
        uint32_t i2 = i1 + (n >> i) / 2;
        if (i2 >= low + len)
          break;
        bool pred = T::ogreater(arr[i1], arr[i2]);
        pred = ObliviousEqual(pred, ascending);
        // These array accesses are oblivious because the indices are
        // deterministic
        T tmp = arr[i1];
        arr[i1] = ObliviousChoose(pred, arr[i2], arr[i1]);
        arr[i2] = ObliviousChoose(pred, tmp, arr[i2]);
      }
    }
  }
}

// Imperative implementation of bitonic sorting network -- works only for powers
// of 2
template <typename T>
inline void imperative_o_sort(T *arr, size_t n, bool ascending) {
  uint32_t i, j, k;
  for (k = 2; k <= n; k = 2 * k) {
    for (j = k >> 1; j > 0; j = j >> 1) {
      for (i = 0; i < n; i++) {
        uint32_t ij = i ^ j;
        if (ij > i) {
          if ((i & k) == 0) {
            bool pred = ObliviousGreater(arr[i], arr[ij]);
            pred = ObliviousEqual(pred, ascending);
            // These array accesses are oblivious because the indices are
            // deterministic
            T tmp = arr[i];
            arr[i] = ObliviousChoose(pred, arr[ij], arr[i]);
            arr[ij] = ObliviousChoose(pred, tmp, arr[ij]);
          } else {
            bool pred = ObliviousGreater(arr[ij], arr[i]);
            pred = ObliviousEqual(pred, ascending);
            // These array accesses are oblivious because the indices are
            // deterministic
            T tmp = arr[i];
            arr[i] = ObliviousChoose(pred, arr[ij], arr[i]);
            arr[ij] = ObliviousChoose(pred, tmp, arr[ij]);
          }
        }
      }
    }
  }
}

// Imperative implementation of bitonic sorting network for POD type -- works
// only for powers of 2 Assumes T implements T::ogreater
template <typename T>
inline void imperative_o_sort_pod(T *arr, size_t n, bool ascending) {
  uint32_t i, j, k;
  for (k = 2; k <= n; k = 2 * k) {
    for (j = k >> 1; j > 0; j = j >> 1) {
      for (i = 0; i < n; i++) {
        uint32_t ij = i ^ j;
        if (ij > i) {
          if ((i & k) == 0) {
            bool pred = T::ogreater(arr[i], arr[ij]);
            pred = ObliviousEqual(pred, ascending);
            // These array accesses are oblivious because the indices are
            // deterministic
            T tmp = arr[i];
            arr[i] = ObliviousChoose(pred, arr[ij], arr[i]);
            arr[ij] = ObliviousChoose(pred, tmp, arr[ij]);
          } else {
            bool pred = T::ogreater(arr[ij], arr[i]);
            pred = ObliviousEqual(pred, ascending);
            // These array accesses are oblivious because the indices are
            // deterministic
            T tmp = arr[i];
            arr[i] = ObliviousChoose(pred, arr[ij], arr[i]);
            arr[ij] = ObliviousChoose(pred, tmp, arr[ij]);
          }
        }
      }
    }
  }
}

// Sort <len> elements in arr -- starting from index arr[low]
template <typename T>
inline void o_sort(T *arr, uint32_t low, uint32_t len, bool ascending) {
  if (len > 1) {
    uint32_t m = greatest_power_of_two_less_than(len);
    if (m * 2 == len) {
      imperative_o_sort(arr + low, len, ascending);
    } else {
      imperative_o_sort(arr + low, m, !ascending);
      o_sort(arr, low + m, len - m, ascending);
      imperative_o_merge(arr, low, len, ascending);
    }
  }
}

// Sort <len> elements in arr of POD type -- starting from index arr[low]
// Assumes T implements T::ogreater
template <typename T>
inline void o_sort_pod(T *arr, uint32_t low, uint32_t len, bool ascending) {
  if (len > 1) {
    uint32_t m = greatest_power_of_two_less_than(len);
    if (m * 2 == len) {
      imperative_o_sort_pod(arr + low, len, ascending);
    } else {
      imperative_o_sort_pod(arr + low, m, !ascending);
      o_sort_pod(arr, low + m, len - m, ascending);
      imperative_o_merge_pod(arr, low, len, ascending);
    }
  }
}

} // namespace detail

template <typename T>
inline void ObliviousMerge(T *arr, uint32_t low, uint32_t len, bool ascending) {
  detail::imperative_o_merge(arr, low, len, ascending);
}

template <typename T>
inline void ObliviousMergePOD(T *arr, uint32_t low, uint32_t len,
                              bool ascending) {
  detail::imperative_o_merge_pod(arr, low, len, ascending);
}

template <typename T>
inline void ObliviousSort(T *arr, uint32_t low, uint32_t len, bool ascending) {
  detail::o_sort(arr, low, len, ascending);
}

template <typename T>
inline void ObliviousSortPOD(T *arr, uint32_t low, uint32_t len,
                             bool ascending) {
  detail::o_sort_pod(arr, low, len, ascending);
}

} // namespace obl
