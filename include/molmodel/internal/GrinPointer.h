/* GrinPointer.h */

/* Portions copyright (c) 2007 Stanford University and Christopher M. Bruns.
 * Contributors:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIdED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef SimTK_QUALIFIEDSMARTPOINTER_H_
#define SimTK_QUALIFIEDSMARTPOINTER_H_

#include "Clonable.h"

// Templated class that uses semantics of a pointer.
// Used as a helper class in the "Cheshire Cat", or PIMPL,
// paradigm used to separate interface from implementation in 
// C++.  

// Very similar to grin_ptr<> class by Alan Griffiths
// http://www.octopull.demon.co.uk/arglib/TheGrin.html

// Uses concept of qualified smart pointer from 
// Kevlin Henney ("Coping with Copying in C++" - Overload 33 ISSN 1354 3172)

// "Qualified" means that the access operators, operator* and
// operator->, are overloaded with respect to const-ness, so
// their results follow the const qualification.

// From Alan Griffiths' article:
// http://accu.org/index.php/journals/482
// 
//        "Let us start with the last point mentioned in the discussion of auto_ptr<>
//        - how to cope with deleting an incomplete type. The destructor can't be a 
//        simple "delete p;" because at the point of instantiation the pointer is to 
//        an "incomplete type" and the compiler won't call the destructor.
//
//        To avoid this I make use of the fact that the constructor for grin_ptr<> 
//        is instantiated in the implementation file, where the class definition for 
//        implementation resides. At this point I force the compiler to generate a 
//        deletion function using a trick I first saw in Richard Hickey's article 
//        "Callbacks in C++ Using Template Functors" (C ++ Gems ISBN 1 884842 37 2): 
//        the constructor stores a pointer to this function in the grin_ptr<>. This 
//        provides a safe method for the destructor to delete the object it owns. 
//        The point of passing around function pointers instead of the apparently 
//        more natural use of virtual member functions is that everything can be 
//        done "by value" and no dynamic allocation is required.
//
//        A similar function is used for copying the object..."

namespace SimTK {

// Overloaded deepCopy() methods will use "clone()" method if 
// referent is derived from Clonable.  Otherwise it calls new
// with its copy constructor, which should work for non-derived
// classes

// Fallback deep copy for non-Clonable objects
template<class PointeeType>
inline PointeeType* deepCopy(const PointeeType* ptr, const void*) 
{
	return ptr ? new PointeeType(*ptr) : 0;
}

// Deep copy for objects that derive from Clonable
template<class PointeeType>
inline PointeeType* deepCopy(const PointeeType *ptr, const Clonable *)
{
	return ptr ? ptr->clone() : 0;
}

// At instantiation time, this method will decide which of the above to use.
template<class PointeeType>
inline PointeeType* deepCopy(const PointeeType* ptr) 
{
	return deepCopy(ptr, ptr);
}


template<typename PointeeType>
class GrinPointer {
public:
    explicit GrinPointer(PointeeType* pointee)
         : doCopy(&myCopyFunction), ptr(pointee), doDelete(&myDeleteFunction) 
    {
        // doCopy = &myCopyFunction;
        // ptr = pointee;
        // doDelete = &myDeleteFunction;
    }

    GrinPointer(const GrinPointer & src);

    ~GrinPointer() throw() { doDelete(ptr); }

    const PointeeType* operator->() const {return ptr;}
    PointeeType* operator->() {return ptr;}

    const PointeeType& operator*() const {return *ptr;}
    PointeeType& operator*() {return *ptr;}

    GrinPointer& operator=(const GrinPointer& src);

private:
    typedef void (*DeleteFunctionPtr)(PointeeType* ptr);
    typedef PointeeType* (*CopyFunctionPtr)(const PointeeType* ptr);

    CopyFunctionPtr doCopy;
    PointeeType* ptr;
    DeleteFunctionPtr doDelete;

    static void myDeleteFunction(PointeeType* ptr) {
        delete ptr;
    }
    static PointeeType*  myCopyFunction(const PointeeType* ptr) {
        return deepCopy(ptr);
    }

};

// Deep copy uses clone() method if wrapped type inherits the
// "Clonable" interface.

template<typename PointeeType>
inline GrinPointer<PointeeType>::GrinPointer(const GrinPointer& src)
:
	doCopy(src.doCopy),
	ptr(doCopy(src.ptr)),
	doDelete(src.doDelete)
{}
    
template<typename PointeeType>
inline GrinPointer<PointeeType>& GrinPointer<PointeeType>::operator=(const GrinPointer& src)
{
    // Allocate...
	PointeeType* tmp = doCopy(src.ptr);

    // ...before release...
	doDelete(ptr);

    // ...and update
	ptr = tmp;
	return *this;
}

}

#endif // SimTK_QUALIFIEDSMARTPOINTER_H_
