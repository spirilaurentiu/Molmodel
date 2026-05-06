#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>

#include "molmodel/internal/DuMMForceFieldSubsystem.h"
#include "molmodel/internal/common.h"

using namespace SimTK;

// Define unique index types that are only used internally.

// This is the index type for the subset of mobilized bodies that have been
// mentioned at all to this instance of DuMMForceField. These will not
// necessarily all be used in force calculations since there may be bodies
// that don't have any included atoms attached to them.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(DuMMBodyIndex);

// This is the index type for the subset of mobilized bodies that is actually
// involved in force calculations because they have "included atoms" attached.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(DuMMIncludedBodyIndex);

// This is the index type for the subset of atoms that are "atom 1" for any
// bond for which we are going to compute a bond force at run time.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(DuMMBondStarterIndex);

//-----------------------------------------------------------------------------
//                                INDEX PAIR
//-----------------------------------------------------------------------------
template <class T>
class IndexPair {
    public:
    IndexPair() = default;
    IndexPair(T index1, T index2, bool canon = false) {
        ixs[0] = index1;
        ixs[1] = index2;
        if (canon) {
            canonicalize();
        }
    }
    auto operator[](int index) const -> const T& {
        assert(0 <= index && index < 2);
        return ixs[index];
    }
    auto operator[](int index) -> T& {
        assert(0 <= index && index < 2);
        return ixs[index];
    }
    [[nodiscard]] auto isValid() const -> bool {
        return ixs[0].isValid() && ixs[1].isValid();
    }
    void invalidate() {
        ixs[0].invalidate();
        ixs[1].invalidate();
    }
    // canonical is low,high
    void canonicalize() {
        if (ixs[0] > ixs[1]) {
            std::swap(ixs[0], ixs[1]);
        }
    }

    private:
    T ixs[2];
};

template <class T>
static inline auto operator<<(std::ostream& ostream, const IndexPair<T>& indexPair) -> std::ostream& {
    ostream << "(" << (int)indexPair[0] << "," << (int)indexPair[1] << ")";
    return ostream;
}

template <class T>
static inline auto operator<(const IndexPair<T>& indexPair1, const IndexPair<T>& indexPair2) -> bool {
    assert(indexPair1.isValid() && indexPair2.isValid());
    if (indexPair1[0] < indexPair2[0]) {
        return true;
    }
    if (indexPair1[0] > indexPair2[0]) {
        return false;
    }
    return indexPair1[1] < indexPair2[1];
}

using AtomIndexPair = IndexPair<DuMM::AtomIndex>;
using IncludedAtomIndexPair = IndexPair<DuMM::IncludedAtomIndex>;
using AtomClassIndexPair = IndexPair<DuMM::AtomClassIndex>;
using MobodIndexPair = IndexPair<MobilizedBodyIndex>;

//-----------------------------------------------------------------------------
//                               INDEX TRIPLE
//-----------------------------------------------------------------------------
template <class T>
class IndexTriple {
    public:
    IndexTriple() = default;
    IndexTriple(T index1, T index2, T index3, bool canon = false) {
        ixs[0] = index1;
        ixs[1] = index2;
        ixs[2] = index3;
        if (canon) {
            canonicalize();
        }
    }
    auto operator[](int index) const -> const T& {
        assert(0 <= index && index < 3);
        return ixs[index];
    }
    auto operator[](int index) -> T& {
        assert(0 <= index && index < 3);
        return ixs[index];
    }
    [[nodiscard]] auto isValid() const -> bool {
        return ixs[0].isValid() && ixs[1].isValid() && ixs[2].isValid();
    }
    void invalidate() {
        ixs[0].invalidate();
        ixs[1].invalidate();
        ixs[2].invalidate();
    }
    // canonical has 1st number <= last number; middle stays put
    void canonicalize() {
        if (ixs[0] > ixs[2]) {
            std::swap(ixs[0], ixs[2]);
        }
    }

    private:
    T ixs[3];
};

template <class T>
static inline auto operator<<(std::ostream& ostream, const IndexTriple<T>& indexTriple) -> std::ostream& {
    ostream << "(" << (int)indexTriple[0] << "," << (int)indexTriple[1] << "," << (int)indexTriple[2] << ")";
    return ostream;
}

template <class T>
static inline auto operator<(const IndexTriple<T>& indexTriple1, const IndexTriple<T>& indexTriple2) -> bool {
    assert(indexTriple1.isValid() && indexTriple2.isValid());
    if (indexTriple1[0] < indexTriple2[0]) {
        return true;
    }
    if (indexTriple1[0] > indexTriple2[0]) {
        return false;
    }
    if (indexTriple1[1] < indexTriple2[1]) {
        return true;
    }
    if (indexTriple1[1] > indexTriple2[1]) {
        return false;
    }
    return indexTriple1[2] < indexTriple2[2];
}

using AtomIndexTriple = IndexTriple<DuMM::AtomIndex>;
using IncludedAtomIndexTriple = IndexTriple<DuMM::IncludedAtomIndex>;
using AtomClassIndexTriple = IndexTriple<DuMM::AtomClassIndex>;

//-----------------------------------------------------------------------------
//                              INDEX QUAD
//-----------------------------------------------------------------------------
template <class T>
class IndexQuad {
    public:
    IndexQuad() = default;
    IndexQuad(T index1, T index2, T index3, T index4, bool canon = false) {
        ixs[0] = index1;
        ixs[1] = index2;
        ixs[2] = index3;
        ixs[3] = index4;
        if (canon) {
            canonicalize();
        }
    }
    auto operator[](int index) const -> const T& {
        assert(0 <= index && index < 4);
        return ixs[index];
    }
    auto operator[](int index) -> T& {
        assert(0 <= index && index < 4);
        return ixs[index];
    }
    [[nodiscard]] auto isValid() const -> bool {
        return ixs[0].isValid() && ixs[1].isValid() && ixs[2].isValid() && ixs[3].isValid();
    }
    void invalidate() {
        ixs[0].invalidate();
        ixs[1].invalidate();
        ixs[2].invalidate();
        ixs[3].invalidate();
    }

    // canonical has 1st number <= last number; middle two must swap
    // if the outside ones do
    void canonicalize() {
        // Index quad has additional case where 1 == 4 and 2 differs from 3
        if ((ixs[0] > ixs[3]) || ((ixs[0] == ixs[3]) && (ixs[1] > ixs[2]))) {
            std::swap(ixs[0], ixs[3]);
            std::swap(ixs[1], ixs[2]);
        }
    }

    private:
    T ixs[4];
};

template <class T>
static inline auto operator<<(std::ostream& ostream, const IndexQuad<T>& indexQuad) -> std::ostream& {
    ostream << "(" << (int)indexQuad[0] << "," << (int)indexQuad[1] << "," << (int)indexQuad[2] << ","
            << (int)indexQuad[3] << ")";
    return ostream;
}

template <class T>
static inline auto operator<(const IndexQuad<T>& indexQuad1, const IndexQuad<T>& indexQuad2) -> bool {
    assert(indexQuad1.isValid() && indexQuad2.isValid());
    if (indexQuad1[0] < indexQuad2[0]) {
        return true;
    }
    if (indexQuad1[0] > indexQuad2[0]) {
        return false;
    }
    if (indexQuad1[1] < indexQuad2[1]) {
        return true;
    }
    if (indexQuad1[1] > indexQuad2[1]) {
        return false;
    }
    if (indexQuad1[2] < indexQuad2[2]) {
        return true;
    }
    if (indexQuad1[2] > indexQuad2[2]) {
        return false;
    }
    return indexQuad1[3] < indexQuad2[3];
}

using AtomIndexQuad = IndexQuad<DuMM::AtomIndex>;
using IncludedAtomIndexQuad = IndexQuad<DuMM::IncludedAtomIndex>;
using AtomClassIndexQuad = IndexQuad<DuMM::AtomClassIndex>;
