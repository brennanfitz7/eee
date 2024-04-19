def get_most_common(some_dict):
    """
    Return key/value pair where value is highest. 
    """
    
    keys = list(some_dict.keys())
    values = [some_dict[k] for k in keys]
    to_sort = list(zip(values,keys))
    to_sort.sort()
    
    return to_sort[-1][1], to_sort[-1][0]
