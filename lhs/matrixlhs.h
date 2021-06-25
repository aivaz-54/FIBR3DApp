/**
 * @file matrixlhs.h
 * @author Robert Carnell
 * @copyright Copyright (c) 2013, Robert Carnell
 * 
 * @license <a href="http://www.gnu.org/licenses/lgpl.html">GNU Lesser General Public License (LGPL v3)</a>
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MATRIXLHS_H
#define MATRIXLHS_H

#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cstddef>

/**
 * @namespace bclib The bertcarnell template library namespace
 */
namespace bclib {

    // forward declare the iterator
    template <class T, bool ISROWWISE> class matrixlhsIter;
    template <class T, bool ISROWWISE> class matrixlhsConstIter;

    /**
     * Matrixlhs Class
     * @tparam T a generic type of the kind that can be used in std::vector
     */
    template<class T>
    class matrixlhs {
        friend class matrixlhsIter<T, true>; /**< make the class a friend of the row-wise iterator */
        friend class matrixlhsIter<T, false>; /**< make the class a friend of the column-wise iterator */
        friend class matrixlhsConstIter<T, true>; /**< make the class a friend of the row-wise iterator */
        friend class matrixlhsConstIter<T, false>; /**< make the class a friend of the column-wise iterator */
    public:
        typedef typename std::vector<T>::size_type size_type; /**< define the size_type as std::vector */
        typedef typename std::vector<T>::iterator iterator; /**< define iterator from the std::vector internals */
        typedef typename std::vector<T>::const_iterator const_iterator; /**< define the const iterator from the std::vector */
        typedef matrixlhsIter<T, true> rowwise_iterator; /**< an iterator that iterates across rows then down columns */
        typedef matrixlhsIter<T, false> columnwise_iterator; /**< an iterator that iterates down columns then across rows */
        typedef matrixlhsConstIter<T, true> const_rowwise_iterator; /**< a const row-wise iterator */
        typedef matrixlhsConstIter<T, false> const_columnwise_iterator; /**< a const column-wise iterator */
        typedef ptrdiff_t difference_type; /**< define difference_type for consistency with stdlib */
        typedef T value_type; /**< define value_type for consistency with stdlib */
        typedef T * pointer; /**< define a pointer type for consistency with stdlib */
        typedef T & reference; /**< define a reference type for consistency with stdlib */

    private:
        size_type rows; /**< number of rows */
        size_type cols; /**< number of columns */
        std::vector<T> elements; /**< array of elements */
        bool bTranspose; /**< is the matrixlhs transposed from creation */

        /**
         * calculate tne location of the value in the vector holding the matrixlhs values
         * @param row the row location
         * @param col the column location
         * @return the location of the value in the vector holding the matrixlhs values
         */
        size_type calcLocation(const size_type row, const size_type col) {
            return (!bTranspose) ? (row * cols + col) : (col * rows + row);
        }

        /**
         * calculate tne location of the value in the vector holding the matrixlhs values
         * @param row the row location
         * @param col the column location
         * @return the location fo the value in the vector holding the matrixlhs values
         */
        size_type calcLocation(const size_type row, const size_type col) const {
            return (!bTranspose) ? (row * cols + col) : (col * rows + row);
        }
    public:
        /// The number of rows in the matrixlhs

        size_type rowsize() const {
            return rows;
        };

        /// The number of columns in the matrixlhs 

        size_type colsize() const {
            return cols;
        };

        /**
         * matrixlhs element access
         * @note does not check for index in range
         * @param row row index (zero based)
         * @param col column index (zero based)
         * @return a reference to the requested element
         */
        T& operator()(size_type row, size_type col) {
            return elements[calcLocation(row, col)];
        }

        /**
         * matrixlhs element access
         * @note does not check for arguments out of range
         * @param row row index (zero based)
         * @param col column index (zero based)
         * @return a const reference to the requested element
         */
        const T& operator()(size_type row, size_type col) const {
            return elements[calcLocation(row, col)];
        }

        /**
         * matrixlhs element access
         * @throws std::out_of_range from the internal std::vector
         * @param row row index (zero based)
         * @param col column index (zero based)
         * @return a const reference to the requested element
         */
        const T& at(size_type row, size_type col) const {
            return elements.at(calcLocation(row, col));
        }

        /**
         * matrixlhs element access
         * @throws std::out_of_range from the internal std::vector
         * @param row row index (zero based)
         * @param col column index (zero based)
         * @return a reference to the requested element
         */
        T& at(size_type row, size_type col) {
            return elements.at(calcLocation(row, col));
        }

        /**
         * matrixlhs element access
         * @throws std::out_of_range from the internal std::vector
         * @param i vector index (zero based)
         * @return a reference to the requested element
         */
        T& at(size_type loc) {
            return elements.at(loc);
        }

        /**
         * matrixlhs element access
         * @throws std::out_of_range from the internal std::vector
         * @param i vector index (zero based)
         * @return const a reference to the requested element
         */
        const T& at(size_type loc) const {
            return elements.at(loc);
        }

        /// a pointer to the internal data array

        T* data() {
            return elements.data();
        };

        /// get the internal data vector

        std::vector<T> getDataVector() const {
            return elements;
        };

        /// Default Constructor with zero rows and zero columns 
        matrixlhs();

        /**
         * Constructor
         * @param rows the number of rows in the matrixlhs
         * @param cols the number of columns in the matrixlhs
         */
        matrixlhs(size_type rows, size_type cols);

        /**
         * Constructor
         * @param rows the number of rows in the matrixlhs
         * @param cols the number of columns in the matrixlhs
         * @param elementArray an array to use as the initial values
         */
        matrixlhs(size_type rows, size_type cols, const T* elementArray);

        /**
         * Constructor
         * @param rows the number of rows in the matrixlhs
         * @param cols the number of columns in the matrixlhs
         * @param elementVector a std::vector to use as the initial values
         */
        matrixlhs(size_type rows, size_type cols, const std::vector<T> & elementVector);

        /**
         * Copy Constructor
         * @param the matrixlhs to be copied
         */
        matrixlhs(const matrixlhs<T> &);

        /// Destructor
        ~matrixlhs();

        /**
         * Matrixlhs assignment
         * @param right hand side matrixlhs
         * @return the left hand side matrixlhs
         */
        matrixlhs<T>& operator=(const matrixlhs<T>&);

        /**
         * Equality comparison operator
         * @param rhs the right hand side matrixlhs
         * @return true if the matrices are equivalent
         */
        bool operator==(const matrixlhs<T> & rhs) const;

        /**
         * Inequality comparison operator
         * @param rhs the right hand side matrixlhs
         * @return true if the matrices are not equivalent
         */
        bool operator!=(const matrixlhs<T> & rhs) const;

        /**
         * Get a row of the matrixlhs as a std::vector
         * @note does not check to ensure the row is in range
         * @param row the row number
         * @return a vector representation of that row
         */
        std::vector<T> getrow(size_type row) const;

        /**
         * Get a row of the matrixlhs as a std::vector
         * @throws std::out_of_range when the row is not in range
         * @param row the row number
         * @return a vector representation of that row
         */
        std::vector<T> getrow_at(size_type row) const;

        /**
         * Get a row of the matrixlhs as a row matrixlhs
         * @note does not check to ensure argument is in range
         * @param row the row number
         * @return a matrixlhs representation of that row
         */
        matrixlhs<T> getRowMatrixlhs(size_type row) const;

        /**
         * Get a row of the matrixlhs as a row matrixlhs
         * @throws an out of range exception for an argument out of range
         * @param row the row number
         * @return a matrixlhs representation of that row
         */
        matrixlhs<T> getRowMatrixlhs_at(size_type row) const;

        /**
         * get a column of the matrixlhs as a vector
         * @note does not check the array bounds
         * @param col column number
         * @return a vector of the requested column
         */
        std::vector<T> getcol(size_type col) const;

        /**
         * Get a column of the matrixlhs as a vector
         * @throws out_of_range error if the column requested is out of bounds
         * @param col the column number
         * @return a vector of the requested column
         */
        std::vector<T> getcol_at(size_type col) const;

        /**
         * Get a column of the matrixlhs as a column matrixlhs
         * @note does not check if the requested column is in bounds
         * @param col the column number
         * @return a column matrixlhs of the requested column
         */
        matrixlhs<T> getColumnMatrixlhs(size_type col) const;

        /**
         * Get a column of the matrixlhs as a column matrixlhs
         * @throws if the requested column is out of range
         * @param col the column number
         * @return a column matrixlhs of the requested column
         */
        matrixlhs<T> getColumnMatrixlhs_at(size_type col) const;

        /**
         * fill the matrixlhs with a value
         * @param x the value to fill the matrixlhs with
         */
        void fill(const T & x) {
            elements.assign(rows*cols, x);
        };

        /**
         * fill the matrixlhs with a value
         * @param x the value to fill the matrixlhs with
         */
        //void fill(const T x)
        //{
        //    elements.assign(rows*cols, x);
        //};

        /// Clear the matrixlhs to zero rows and columns
        void clear();

        /// return true if the matrixlhs is empty

        bool isEmpty() const {
            return elements.empty();
        };

        /// return a string representation of the matrixlhs
        std::string toString() const;

        /// Transpose the matrixlhs
        void transpose();

        /// return true if this matrixlhs is operating as a transposed matrixlhs from the original definition

        bool isTransposed() const {
            return bTranspose;
        };

        /********* Matrixlhs Iterators *********/

        /// an iterator for the beginning of the internal vector

        iterator begin() {
            return elements.begin();
        };

        const_iterator begin() const {
            return elements.begin();
        };

        /// An iterator for one iteration past the end of the internal vector

        iterator end() {
            return elements.end();
        };

        const_iterator end() const {
            return elements.end();
        };

        /// An iterator that operates along the matrixlhs rows 

        rowwise_iterator rowwisebegin() {
            return rowwise_iterator(*this, 0, 0);
        };

        const_rowwise_iterator rowwisebegin() const {
            return const_rowwise_iterator(*this, 0, 0);
        };

        /**
         * return a row wise iterator for the beginning of the ith row (0 based)
         * @param irow
         */
        rowwise_iterator rowwisebegin(size_type irow) {
            return rowwise_iterator(*this, irow, 0);
        };

        const_rowwise_iterator rowwisebegin(size_type irow) const {
            return const_rowwise_iterator(*this, irow, 0);
        };

        /// An iterator that operates along the matrixlhs row

        rowwise_iterator rowwiseend() {
            return rowwise_iterator(*this, rows, 0);
        };

        const_rowwise_iterator rowwiseend() const {
            return const_rowwise_iterator(*this, rows, 0);
        };

        /**
         * return a row wise iterator for the end of the ith row (0 based)
         * @param irow
         */
        rowwise_iterator rowwiseend(size_type irow) {
            return rowwise_iterator(*this, irow + 1, 0);
        };

        const_rowwise_iterator rowwiseend(size_type irow) const {
            return const_rowwise_iterator(*this, irow + 1, 0);
        };

        /// An iterator that operates along the matrixlhs columns

        columnwise_iterator columnwisebegin() {
            return columnwise_iterator(*this, 0, 0);
        };

        const_columnwise_iterator columnwisebegin() const {
            return const_columnwise_iterator(*this, 0, 0);
        };

        /**
         * return a column wise iterator for the beginning of the jth column (0 based)
         * @param irow
         */
        columnwise_iterator columnwisebegin(size_type jcol) {
            return columnwise_iterator(*this, 0, jcol);
        };

        const_columnwise_iterator columnwisebegin(size_type jcol) const {
            return const_columnwise_iterator(*this, 0, jcol);
        };

        /// An iterator that operates along the matrixlhs columns

        columnwise_iterator columnwiseend() {
            return columnwise_iterator(*this, 0, cols);
        };

        const_columnwise_iterator columnwiseend() const {
            return const_columnwise_iterator(*this, 0, cols);
        };

        /**
         * return a column wise iterator for the end of the jth column (0 based)
         * @param irow
         */
        columnwise_iterator columnwiseend(size_type jcol) {
            return columnwise_iterator(*this, 0, jcol + 1);
        };

        const_columnwise_iterator columnwiseend(size_type jcol) const {
            return const_columnwise_iterator(*this, 0, jcol + 1);
        };
    };

    /******************************************************************************/

    /**
     * An iterator class for the <code>matrixlhs</code> class
     * @tparam T the type of object stored in the matrixlhs
     * @tparam ISROWWISE a boolean to indicate if the matrixlhs is iterated row-wise
     */
    template <class T, bool ISROWWISE>
    class matrixlhsIter : public std::iterator<std::forward_iterator_tag, T> {
        friend class matrixlhsConstIter<T, ISROWWISE>;
    private:
        matrixlhs<T> & myMatrixlhs; /**< The object that the iterator is referencing */
        typename matrixlhs<T>::size_type rows; /**< the row being pointed to */
        typename matrixlhs<T>::size_type cols; /**< the column being pointed to */
    public:

        /**
         * Constructor
         * @param mat the matrixlhs being indexed
         * @param r the row location of the iterator
         * @param c the column location of the iterator
         */
        matrixlhsIter(matrixlhs<T> & mat, typename matrixlhs<T>::size_type r,
                typename matrixlhs<T>::size_type c)
        : myMatrixlhs(mat), rows(r), cols(c) {
        }
        /// Equality operator
        bool operator==(const matrixlhsIter<T, ISROWWISE> & other) const;
        /// Inequality operator

        bool operator!=(const matrixlhsIter<T, ISROWWISE> & other) const {
            return !(*this == other);
        }
        /// pre-increment operator
        matrixlhsIter<T, ISROWWISE> & operator++();
        /// post-increment operator
        matrixlhsIter<T, ISROWWISE> operator++(int);
        /// assignment operator
        matrixlhsIter<T, ISROWWISE> & operator=(const matrixlhsIter<T, ISROWWISE> & rhs);
        /// de-reference operator

        T & operator*() {
            return myMatrixlhs(rows, cols);
        }
    };

    /**
     * An const_iterator class for the <code>matrixlhs</code> class
     * @tparam T the type of object stored in the matrixlhs
     * @tparam ISROWWISE a boolean to indicate if the matrixlhs is iterated row-wise
     */
    template <class T, bool ISROWWISE>
    class matrixlhsConstIter : public std::iterator<std::forward_iterator_tag, T> {
        friend class matrixlhsIter<T, ISROWWISE>;
    private:
        const matrixlhs<T> & myMatrixlhs; /**< The object that the iterator is referencing */
        typename matrixlhs<T>::size_type rows; /**< the row being pointed to */
        typename matrixlhs<T>::size_type cols; /**< the column being pointed to */
    public:

        /**
         * Constructor
         * @param mat the matrixlhs being indexed
         * @param r the row location of the iterator
         * @param c the column location of the iterator
         */
        matrixlhsConstIter(const matrixlhs<T> & mat, typename matrixlhs<T>::size_type r,
                typename matrixlhs<T>::size_type c)
        : myMatrixlhs(mat), rows(r), cols(c) {
        }

        /**
         * Copy constructor from non-const to const
         * @param mi the matrixlhs being copied
         */
        matrixlhsConstIter(const matrixlhsIter<T, ISROWWISE> & mi)
        : myMatrixlhs(mi.myMatrixlhs), rows(mi.rows), cols(mi.cols) {
        }
        /// Equality operator
        bool operator==(const matrixlhsConstIter<T, ISROWWISE> & other) const;
        /// Inequality operator

        bool operator!=(const matrixlhsConstIter<T, ISROWWISE> & other) const {
            return !(*this == other);
        }
        /// pre-increment operator
        matrixlhsConstIter<T, ISROWWISE> & operator++();
        /// post-increment operator
        matrixlhsConstIter<T, ISROWWISE> operator++(int);
        /// Assignment operator
        /** @TODO:  does an assignment operator make sense for a const iterator? */
        matrixlhsConstIter<T, ISROWWISE> & operator=(const matrixlhsConstIter<T, ISROWWISE> & rhs);
        /// de-reference operator

        const T & operator*() {
            return myMatrixlhs(rows, cols);
        }
    };

    // heavily influenced by:  http://www.sj-vs.net/c-implementing-const_iterator-and-non-const-iterator-without-code-duplication/

    /******************************************************************************/

    template<class T>
    matrixlhs<T>::matrixlhs(size_type rows, size_type cols)
    : rows(rows), cols(cols), bTranspose(false) {
        if (rows == 0 || cols == 0) {
            throw std::range_error("attempt to create a degenerate matrixlhs");
        }
        elements = std::vector<T>(rows * cols);
    }

    template<class T>
    matrixlhs<T>::matrixlhs(size_type rows, size_type cols, const T* elementArray)
    : rows(rows), cols(cols), bTranspose(false) {
        if (rows == 0 || cols == 0) {
            throw std::range_error("attempt to create a degenerate matrixlhs");
        }
        // initialize from array
        elements = std::vector<T>(rows * cols);

        for (size_t i = 0; i < rows * cols; i++) {
            elements[i] = elementArray[i];
        }
    }

    template<class T>
    matrixlhs<T>::matrixlhs(size_type rows, size_type cols, const std::vector<T> & elementVector)
    : rows(rows), cols(cols), bTranspose(false) {
        if (rows == 0 || cols == 0) {
            throw std::range_error("attempt to create a degenerate matrixlhs");
        }
        if (elementVector.size() != rows * cols) {
            throw std::range_error("Input element Vector is not the right size");
        }
        elements.assign(elementVector.begin(), elementVector.end());
    }

    template<class T>
    matrixlhs<T>::matrixlhs(const matrixlhs<T> & cp)
    : rows(cp.rows), cols(cp.cols), elements(cp.elements), bTranspose(cp.bTranspose) {
    }

    template<class T>
    matrixlhs<T>::~matrixlhs() {
    }

    template<class T>
    matrixlhs<T>& matrixlhs<T>::operator=(const matrixlhs<T>& cp) {
        if (cp.rows != rows || cp.cols != cols) {
            rows = cp.rows;
            cols = cp.cols;
        }
        elements = cp.elements;
        bTranspose = cp.bTranspose;
        return *this;
    }

    template<class T>
    bool matrixlhs<T>::operator==(const matrixlhs<T>& cp) const {
        if (cp.rows != rows || cp.cols != cols) {
            return false;
        }
        return std::equal(elements.begin(), elements.end(), cp.elements.begin());
    }

    template<class T>
    bool matrixlhs<T>::operator!=(const matrixlhs<T> & cp) const {
        if (*this == cp) {
            return false;
        }
        return true;
    }

    template<class T>
    std::vector<T> matrixlhs<T>::getrow(size_type row) const {
        std::vector<T> a = std::vector<T>(cols);
        for (size_type j = 0; j < cols; j++) {
            a[j] = elements[calcLocation(row, j)];
        }
        return a;
    }

    template<class T>
    std::vector<T> matrixlhs<T>::getrow_at(size_type row) const {
        if (row >= rows) {
            std::ostringstream msg;
            msg << "row " << row << " was requested, but the matrixlhs has " << rows << " rows";
            throw std::out_of_range(msg.str().c_str());
        }
        return getrow(row);
    }

    template<class T>
    matrixlhs<T> matrixlhs<T>::getRowMatrixlhs(size_type row) const {
        // the simple method has an extra loop of assignment
        //std::vector<T> a = this->getrow(i);
        //return matrixlhs<T>(1,cols,a);
        matrixlhs<T> a(1, cols);
        for (size_type j = 0; j < cols; j++) {
            a(0, j) = elements[calcLocation(row, j)];
        }
        return a;
    }

    template<class T>
    matrixlhs<T> matrixlhs<T>::getRowMatrixlhs_at(size_type row) const {
        if (row >= rows) {
            std::ostringstream msg;
            msg << "Row " << row << " was requested, but the matrixlhs has " << rows << " rows";
            throw std::out_of_range(msg.str().c_str());
        }
        return getRowMatrixlhs(row);
    }

    template<class T>
    std::vector<T> matrixlhs<T>::getcol(size_type col) const {
        std::vector<T> a = std::vector<T>(rows);
        for (size_type i = 0; i < rows; i++) {
            a[i] = elements[calcLocation(i, col)];
        }
        return a;
    }

    template<class T>
    std::vector<T> matrixlhs<T>::getcol_at(size_type col) const {
        if (col >= cols) {
            std::ostringstream msg;
            msg << "Column " << col << " was requested, but the matrixlhs has " << cols << " columns";
            throw std::out_of_range(msg.str().c_str());
        }
        return getcol(col);
    }

    template<class T>
    matrixlhs<T> matrixlhs<T>::getColumnMatrixlhs(size_type col) const {
        matrixlhs<T> a(rows, 1);
        for (size_type i = 0; i < rows; i++) {
            a(i, 0) = elements[calcLocation(i, col)];
        }
        return a;
    }

    template<class T>
    matrixlhs<T> matrixlhs<T>::getColumnMatrixlhs_at(size_type col) const {
        if (col >= cols) {
            std::ostringstream msg;
            msg << "Column " << col << " was requested, but the matrixlhs has " << cols << " columns";
            throw std::out_of_range(msg.str().c_str());
        }
        return getColumnMatrixlhs(col);
    }

    template<class T>
    void matrixlhs<T>::clear() {
        elements.clear();
        rows = 0;
        cols = 0;
        bTranspose = false;
    }

    template<class T>
    matrixlhs<T>::matrixlhs() {
        rows = 0;
        cols = 0;
        elements = std::vector<T>();
        bTranspose = false;
    }

    template<class T>
    std::string matrixlhs<T>::toString() const {
        std::ostringstream msg;
        for (size_type irow = 0; irow < rows; irow++) {
            for (size_type jcol = 0; jcol < cols; jcol++) {
                msg << (*this).at(irow, jcol);
                if (cols > 1 && jcol < cols - 1) {
                    msg << ",";
                }
            }
            msg << "\n";
        }
        return msg.str();
    }

    template<class T>
    void matrixlhs<T>::transpose() {
        // decide to not move data during transpose
        bTranspose = !bTranspose;
        size_type oldRows = rows;
        rows = cols;
        cols = oldRows;
    }

    /******************************************************************************/

    template<class T, bool ISROWWISE>
    bool matrixlhsIter<T, ISROWWISE>::operator==(const matrixlhsIter<T, ISROWWISE> & other) const {
        if (this->myMatrixlhs == other.myMatrixlhs &&
                this->rows == other.rows &&
                this->cols == other.cols) {
            return true;
        }
        return false;
    }

    template<class T, bool ISROWWISE>
    matrixlhsIter<T, ISROWWISE> & matrixlhsIter<T, ISROWWISE>::operator++() {
        if (ISROWWISE) {
            if (cols < myMatrixlhs.cols - 1) {
                cols++;
                return *this;
            } else {
                cols = 0;
                rows++;
                return *this;
            }
        } else // ISROWWISE = false
        {
            if (rows < myMatrixlhs.rows - 1) {
                rows++;
                return *this;
            } else {
                rows = 0;
                cols++;
                return *this;
            }
        }
    }

    template<class T, bool ISROWWISE>
    matrixlhsIter<T, ISROWWISE> & matrixlhsIter<T, ISROWWISE>::operator=(const matrixlhsIter<T, ISROWWISE> & rhs) {
        // Check for self-assignment
        if (this == &rhs) {
            return *this;
        } else {
            this->myMatrixlhs = rhs.myMatrixlhs;
            this->rows = rhs.rows;
            this->cols = rhs.cols;
            return *this;
        }
    }

    template<class T, bool ISROWWISE>
    matrixlhsIter<T, ISROWWISE> matrixlhsIter<T, ISROWWISE>::operator++(int) {
        const matrixlhsIter<T, ISROWWISE> clone(*this);
        ++(*this);
        return clone;
    }

    /******************************************************************************/

    template<class T, bool ISROWWISE>
    bool matrixlhsConstIter<T, ISROWWISE>::operator==(const matrixlhsConstIter<T, ISROWWISE> & other) const {
        if (this->myMatrixlhs == other.myMatrixlhs &&
                this->rows == other.rows &&
                this->cols == other.cols) {
            return true;
        }
        return false;
    }

    template<class T, bool ISROWWISE>
    matrixlhsConstIter<T, ISROWWISE> & matrixlhsConstIter<T, ISROWWISE>::operator++() {
        if (ISROWWISE) {
            if (cols < myMatrixlhs.cols - 1) {
                cols++;
                return *this;
            } else {
                cols = 0;
                rows++;
                return *this;
            }
        } else // ISROWWISE = false
        {
            if (rows < myMatrixlhs.rows - 1) {
                rows++;
                return *this;
            } else {
                rows = 0;
                cols++;
                return *this;
            }
        }
    }

    template<class T, bool ISROWWISE>
    matrixlhsConstIter<T, ISROWWISE> & matrixlhsConstIter<T, ISROWWISE>::operator=(const matrixlhsConstIter<T, ISROWWISE> & rhs) {
        // Check for self-assignment
        if (this == &rhs) {
            return *this;
        } else {
            this->myMatrixlhs = rhs.myMatrixlhs;
            this->rows = rhs.rows;
            this->cols = rhs.cols;
            return *this;
        }
    }

    template<class T, bool ISROWWISE>
    matrixlhsConstIter<T, ISROWWISE> matrixlhsConstIter<T, ISROWWISE>::operator++(int) {
        const matrixlhsConstIter<T, ISROWWISE> clone(*this);
        ++(*this);
        return clone;
    }

} // end namespace

#endif /* MATRIXLHS_H */
