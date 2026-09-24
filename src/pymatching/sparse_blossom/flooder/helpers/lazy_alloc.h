#ifndef PYMATCHING2_LAZY_ALLOC_H
#define PYMATCHING2_LAZY_ALLOC_H

#include <sys/mman.h>

#include <cstddef>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

namespace pm {

/// Allocator for the global-indexed node arrays. In its default (non-lazy) mode it behaves like
/// std::allocator. In lazy mode:
///   - storage comes from an anonymous mmap(MAP_NORESERVE): a page is only backed by RAM once written, and
///     reads of never-written pages see zeros;
///   - value-initialising construct(p) is a no-op, so vector::resize(n) writes nothing.
/// The array is still one contiguous block, so element i lives at data() + i exactly as with std::vector.
/// A slot is only valid after it has been explicitly placement-new'd; untouched slots are all-zero bytes
/// (which destroy safely: every member is a vector or pointer). Only the slots a caller constructs become
/// resident.
template <typename T>
struct LazyNodeAllocator {
    using value_type = T;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::false_type;

    bool lazy{false};

    LazyNodeAllocator() = default;
    explicit LazyNodeAllocator(bool lazy_) : lazy(lazy_) {
    }
    template <typename U>
    LazyNodeAllocator(const LazyNodeAllocator<U>& other) noexcept : lazy(other.lazy) {
    }

    T* allocate(std::size_t n) {
        if (!lazy) {
            return std::allocator<T>().allocate(n);
        }
        std::size_t bytes = n * sizeof(T);
        void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        if (p == MAP_FAILED) {
            throw std::bad_alloc();
        }
        return static_cast<T*>(p);
    }

    void deallocate(T* p, std::size_t n) noexcept {
        if (!lazy) {
            std::allocator<T>().deallocate(p, n);
        } else {
            munmap(p, n * sizeof(T));
        }
    }

    /// Default-construction: skipped in lazy mode (the slot stays all-zero until placement-new'd).
    template <typename U>
    void construct(U* p) {
        if (!lazy) {
            ::new (static_cast<void*>(p)) U();
        }
    }
    template <typename U, typename... Args>
    void construct(U* p, Args&&... args) {
        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
    }

    template <typename U>
    bool operator==(const LazyNodeAllocator<U>& other) const noexcept {
        return lazy == other.lazy;
    }
    template <typename U>
    bool operator!=(const LazyNodeAllocator<U>& other) const noexcept {
        return lazy != other.lazy;
    }
};

/// A std::vector of nodes that can be created lazily (see LazyNodeAllocator).
template <typename T>
using NodeVec = std::vector<T, LazyNodeAllocator<T>>;

}  // namespace pm

#endif  // PYMATCHING2_LAZY_ALLOC_H
