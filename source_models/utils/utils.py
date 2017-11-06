import numpy as np

def largest_remainder(weights, expected_sum=1, precision=0):
    """Use largest remainder method to round weights such that
    the sum to 1.
    params weights: The raw weights
    params expected_sum: The total the weights should sum to
    params precision: Number of decimal places to round to
    """
    # Ensure is array
    if type(weights) == list:
        weights = np.array(weights)
    total_number = expected_sum*np.power(10,precision)
    weights = weights*np.power(10,precision)
    updated_weights = np.floor(weights)
    remainders = weights - updated_weights
    unallocated_places = total_number - np.sum(updated_weights)
    for i in range(int(unallocated_places)):
        max_remainder_index = np.argmax(remainders)
        updated_weights[max_remainder_index] = updated_weights[max_remainder_index] + 1
        remainders[max_remainder_index] = 0
    updated_weights = updated_weights/np.power(10, precision)    

    return updated_weights
