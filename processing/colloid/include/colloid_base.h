#ifndef COLLOID_BASE_H
#define COLLOID_BASE_H

#include <iostream>

using namespace std;

/*! base class for colloid
 *
 */
class colloid_base
{
	private:
		int ncol; //!< column size
		int nrow; //!< row size
		int total; //!< ncol*nrow, total size of the data
		float * array; //!< store the array data
	public:
		
		//! constructor
		colloid_base();
		//! create an uninitialized array
		colloid_base(int, int);
		//! create an array initialized to a given value
		colloid_base(int, int, float);
		//! create an array initialized by a given array
		colloid_base(int, int, float*);
		//! copy constructor
		colloid_base(const colloid_base&);

		//! destructor
		~colloid_base();

		/*! operator overload */
		void operator=(const colloid_base&);
		colloid_base operator+(const colloid_base&);
		colloid_base operator-(const colloid_base&);
		
		/* functions */
		
		/*! set value to an initialized colloid_base type variable,
		 *  the variable should be constructed by
		 *  	colloid_base cb;
		 * @param nc
		 * 		column number
		 * @param nr
		 * 		row number
		 * @param value
		 * 		set the array with value
		 */
		void set_value(int, int, float);

		/*! set values to an initialized colloid_base type variable,
		 *  the variable should be constructed by
		 *  	colloid_base cb;
		 * @param nc
		 * 		column number
		 * @param nr
		 * 		row number
		 * @param arr
		 * 		array to store
		 */
		void set_value(int, int, float*);

		/*! reserve memory to an initialized colloid_base type variable,
		 *  the variable should be constructed by
		 *  	colloid_base cb;
		 *  @param nc
		 *  	column number
		 *  @param nr
		 *  	row number
		 */
		void reserve_memory(int, int);

		inline void free_memory()
		{
			delete [] array;
			array=new float [1];
		}
		
		
		/*! get value by index
		 *	@param icol 
		 *		index of column
		 *	@param irow 
		 *		index of row
		 *	@return 
		 *		value stored in the array labled by (icol, irow)
		 */
		float v(int, int);

		/*! get a column values
		 * @param store
		 *		values of a given column will be stored
		 * @param icol
		 * 		index of column
		 */
		void get_column(float*, int);

		/*! get a row values
		 * @param store
		 *		values of a given row will be stored
		 * @param icol
		 * 		index of row
		 */
		void get_row(float*, int);

		/*! get the size
		 * @return
		 * 		size[0]=ncol; size[1]=nrow; size[2]=total;
		 */
		int * get_size() const;
		//void get_size(int**);

		inline int get_total() const
		{
			return total;
		}
		
		/* display size information */
		inline void show_size() const
		{
			cout << "Array: " << ncol << 'x' << nrow
				 << " (=" << total << ")\n";
		}
		
		/*! get the array pointer
		 * @return
		 * 		the array pointer (float)
		 */
		//float* get_array_pointer(); //put inline when defined
		inline float* get_array_pointer() const { return array; }
		// inline function of class member should be put in header file
};




#endif /* COLLOID_H */
