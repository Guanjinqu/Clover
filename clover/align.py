"""A module for placing the global matching algorithm.

This module is a global matching algorithm for Clover to enable global matching mode. 
Note that Clover itself does not provide a complete global matching algorithm, 
as there are many mature global matching algorithms available. Therefore, this module 
currently only provides a simple global matching example algorithm to demonstrate the 
most basic global matching functionality. If you want to use global matching mode, we 
recommend that you replace the global matching algorithm with the following function format.

typical usage example:


"""
def global_align(read_1,read_2):
    """Global Matching Algorithm

    Perform a global comparison of the two input sequences and output a list of errors.
    Note: The replaced global comparison algorithm must have the same input and output 
    formats to ensure that Clover runs correctly.

    Args:
        read_1: str,The first sequence that is compared and is in the core set of sequences.
        read_2: str,The second sequence of the comparison and is the post-match sequence.

    Returns:
        Returns a list with elements in tuple format, each tuple contains two elements, the 
        position that does not match, and the base at that position in read_2.
    """
    
    error_list=[]
    for i in range(len(read_1)):
        if read_1[i] == read_2[i]:
            pass
        else :
            error_list.append((i,read_2[i]))
    return error_list


