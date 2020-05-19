""" A simple decorator that times the duration of a function's execution. More info on Decorators at https://pythonconquerstheuniverse.wordpress.com/2009/08/06/introduction-to-python-decorators-part-1/"""
import timeit
from functools import wraps

def timer(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        
        start_time = timeit.default_timer()
        result = function(*args, **kwargs)
        elapsed = timeit.default_timer() - start_time
        print('Function "{name}" took {time} seconds to complete.'.format(name=function.__name__, time=elapsed))
        return result
    return wrapper
