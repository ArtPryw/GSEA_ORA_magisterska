import time

def measure_time(func, *args, **kwargs):
    """
    Measures the execution time of a function.

    Parameters:
    func (callable): The function to measure.
    *args: Arguments to pass to the function.
    **kwargs: Keyword arguments to pass to the function.

    Returns:
    result: The result of the function execution.
    execution_time (float): The execution time in seconds.
    """
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"{func.__name__} execution time: {execution_time:.4f} seconds")
    return result, execution_time
