#pragma once

#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

template<typename T>
class Span {
public:
    using element_type = T;
    using value_type = typename std::remove_cv<T>::type;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using iterator = pointer;
    using const_iterator = const T*;

    // Constructors
    constexpr Span() noexcept : ptr_(nullptr), size_(0) {}
    constexpr Span(pointer ptr, size_type count) noexcept : ptr_(ptr), size_(count) {}
    constexpr Span(pointer first, pointer last) noexcept : ptr_(first), size_(last - first) {}

	template<typename It,
			typename = typename std::enable_if<
				std::is_pointer<typename std::iterator_traits<It>::pointer>::value
			>::type>
	constexpr Span(It first, It last) noexcept
		: ptr_(&*first), size_(static_cast<size_type>(last - first)) {}

    template<typename Container,
             typename = typename std::enable_if<
                 std::is_same<typename std::remove_pointer<decltype(std::declval<Container>().data())>::type,
                              value_type>::value
             >::type>
    constexpr Span(Container& cont) noexcept : ptr_(cont.data()), size_(cont.size()) {}

    // Iterators
    constexpr iterator begin() const noexcept { return ptr_; }
    constexpr iterator end() const noexcept { return ptr_ + size_; }
    constexpr const_iterator cbegin() const noexcept { return ptr_; }
    constexpr const_iterator cend() const noexcept { return ptr_ + size_; }

    // Element access
    constexpr reference operator[](size_type idx) const noexcept { return ptr_[idx]; }
    constexpr reference front() const noexcept { return ptr_[0]; }
    constexpr reference back() const noexcept { return ptr_[size_ - 1]; }
    constexpr pointer data() const noexcept { return ptr_; }

    // Capacity
    constexpr size_type size() const noexcept { return size_; }
    constexpr bool empty() const noexcept { return size_ == 0; }

    // Subviews
    constexpr Span<T> first(size_type count) const noexcept { return Span<T>(ptr_, count); }
    constexpr Span<T> last(size_type count) const noexcept { return Span<T>(ptr_ + (size_ - count), count); }
    constexpr Span<T> subSpan(size_type offset, size_type count = static_cast<size_type>(-1)) const noexcept {
        if (count == static_cast<size_type>(-1)) count = size_ - offset;
        return Span<T>(ptr_ + offset, count);
    }

private:
    pointer ptr_;
    size_type size_;
};

template <typename T>
auto safe_subspan(std::vector<T>& vec, size_t start, size_t end) {
    if (start >= end || start >= vec.size()) {
        return Span<T>(); // Return empty span
    }
    return Span<T>(&vec[start], end - start);
}
