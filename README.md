## covi

The Compressed Overlap Index Data Structure (or covi) implementation. 
<!-- contains the C++11 codes associated with the following reference: 
  R. Canovas, B. Cazaux and E. Rivals. "The Compressed Overlap Index". --> 


This implementation creates an index for all the overlaps of a given 
set of words. The set of words must be given as a single input file 
containing all the words separated by a unique symbol choose by the user.

## Compile

To be able to compile the covi codes: 
- Install sdsl-lite. Follow the installation guide here: (https://github.com/simongog/sdsl-lite)
- Modify the location of the sdsl library in the CMakeLists.txt if necessary.
- Go to the build folder and run: 
	- cmake ..
	- make


## Methods

-[CreateCOvI] (https://github.com/tests):

	Use: ./CreateCOvI <file_name> <opt>
      		<file_name>: Name of the file to be use to create the required data structure 
          	<opt>: 
			-o output_name:  String containing the name of the output file. Default file_name.covi
			-s separator:  (ASCII) Value used as separator between words. Default separator = 0
			-w Index_type. Default = 0
			  |  Index_type
			---+--------------------
			0 | COvI (Compressed Overlap Index)
        		1 | fullAC (Aho-Corasick)

          	output:  covi/fullAC file

		Example: ./CreateCOvI ./data/file -o file.covi -s 10 -w 0
		output:  file.covi  or file.fullac
        

-[TestOperations] (https://github.com/tests):

	Use: ./TestOperations <index_file_name> -w type

		output: Displays times per operation in the output

		Example: ./TestOperations file.covi  -w 0
		output:
			Max-ov(x,y) Time: 5.232 microsec
			Correlation(x,y) Time: 5.212 microsec
			all_right_overlaps(x) Time: 327.2 millis
			all_left_overlaps(x) Time: 615.8 millis
			global_max_ovelap() Time: 2841 millis
		        
Note: These codes assume that the computer have enough RAM memory to create the Trie of the set of words
