import numpy as np

def linear_fit(x, y, fit_min, fit_max):
    """Fit x, y with a linear function

    y = mx + b

    Args:
        x: x variable
        y: y variable
        fit_min: minimal value of x to fit
        fit_max: maximal value of x to fit

    Returns:
        m: slope of the fitted line
        b: intercept of the fitted function
    """
    idx = (x > fit_min) & (x < fit_max)
    x = np.vstack([x[idx], np.ones_like(x[idx])]).T
    y = y[idx]    
    m, b = np.linalg.lstsq(x, y, rcond=None)[0]
    return m, b
