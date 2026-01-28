// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#ifndef PYMATCHING2_VECTORWRAPPER_H
#define PYMATCHING2_VECTORWRAPPER_H

// Provides vector-like functions for an array
#include <type_traits>

template <typename T> struct VectorWrapper {
public:
    T *arr_;
    std::size_t size_;
    std::size_t capacity_;
    VectorWrapper()
        : arr_(nullptr), size_(0), capacity_(0) {};
    VectorWrapper(T *arr, std::size_t capacity)
        : arr_(arr), size_(0), capacity_(capacity) {}
    ~VectorWrapper() {
        destroy_elements(std::integral_constant<bool, std::is_trivially_destructible<T>::value>{});
    }
    VectorWrapper(const VectorWrapper&) = delete;
    VectorWrapper& operator=(const VectorWrapper&) = delete;

    VectorWrapper(VectorWrapper&& other) noexcept
        : arr_(other.arr_), size_(other.size_), capacity_(other.capacity_) {
        other.arr_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }
    VectorWrapper& operator=(VectorWrapper&& other) noexcept {
        if (this != &other) {
            destroy_elements(std::integral_constant<bool, std::is_trivially_destructible<T>::value>{});
            arr_ = other.arr_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.arr_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    inline std::size_t size() const { return size_; };
    inline bool empty() const { return size_ == 0; };

    inline T& operator[](std::size_t i) { return arr_[i]; }
    inline const T& operator[](std::size_t i) const { return arr_[i]; }

    inline T *begin() { return arr_; };
    inline T* end()   { return arr_ + size_; }
    inline const T* begin() const { return arr_; }
    inline const T* end()   const { return arr_ + size_; }

    inline void push_back(const T& value) {
        if (size_ >= capacity_) {
            throw std::invalid_argument("VectorWrapper::push_back() overflow (size_ not less than capacity_)");
        }
        arr_[size_++] = value;
    }

    template<class... Args>
    inline T& emplace_back(Args&&... args) {
        if (size_ >= capacity_) {
            throw std::invalid_argument("VectorWrapper::emplace_back() overflow (size_ not less than capacity_)");
        }
        T* p = &arr_[size_++];
        new (p) T(std::forward<Args>(args)...);
        return *p;
    }

    // insert() assumes that T supports direct copying.
    // Note: num_elems != 1 is not supported
    inline T* insert(T* pos, size_t num_elems, T data) {
        if (num_elems != 1) {
            throw std::invalid_argument("VectorWrapper::insert() num_elems != 1 is not supported");
        }
        if (size_ >= capacity_) {
            throw std::invalid_argument("VectorWrapper::insert() overflow (size_ not less than capacity_)");
        }
        if (pos < arr_ || pos > arr_+size_) {
            throw std::invalid_argument("VectorWrapper::insert() OOB (pos is out of bounds)");
        }
        for (T* p = arr_+size_; p > pos; --p) {
            *p = *(p-1); 
        }
        *pos = data;
        size_++;
        return pos;
    }

private:
    void destroy_elements(std::true_type) {}
    void destroy_elements(std::false_type) {
        for (std::size_t i = 0; i < size_; ++i) {
            arr_[i].~T();
        }
        size_ = 0;
    }
};

#endif