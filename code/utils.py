import numpy as np
from scipy import optimize, interpolate
from tqdm import tqdm, trange

def generate_line(b, m, start, end):
    """
    generate_line Generates a line of the form y = mx + b between start and end.

    :param b: Where it would cross y axis
    :param m: Slope
    :param start: Start position
    :param end: End position
    :return: A list with y values
    """
    return [m*x + b for x in np.arange(0,end-start,0.01)]


def calculate_total_length(x, y, delta):
    """
    calculate_total_length Calculates the polygon chain length for input data x,y with divider width of delta.

    :param x: Horizontal indices (same length as y)
    :param y: Data values
    :param delta: Divider width
    :return: Returns a tuple with divider width, polygon length and regression parameters to draw lines
    """
    delta_x = (x[delta] - x[0])
    regression_params = [] # Parameters for each ruler line
    polylengths = []       # Length of each line
    done = False
    #By only using index, the lowest resolution is by incrementing by 1, because smallest increment is x[i+1]-x[i]

    # Initial values
    x_1 = x[0]
    x_2 = x[delta]
    y_1 = y[0]
    y_2 = y[delta]

    # Fixed loop size with equal spacing
    for i in range((len(x)//delta)):

        hyp = np.sqrt(np.float_power(y_2-y_1,2)+np.float_power(x_2-x_1,2)) # Hypothenuse i.e. length of section line

        # Calculate m and beta
        m = (y_2 - y_1)/(x_2 - x_1)
        beta =  y_1

        regression_params.append((beta, m, x_1, x_2)) # Add line parameters to list

        polylengths.append(hyp) # Add length of section to the list of lengths


        if done: break
        # Update variables for new iteration
        x_1 = x_2
        y_1 = y_2
        if ((i+2)*delta < (len(x)-1)):
            x_2 = x[(i+2)*delta]
            y_2 = y[(i+2)*delta]
        else:
            x_2 = x[len(x)-1]
            y_2 = y[len(y)-1]
            done = True

    L_lambda = np.sum(polylengths) # Loop done, calculate length of all the segement lines

    return((delta_x, L_lambda, regression_params))


def segments_fit(X, Y, num_knots):
    """
    Piecewise regression on input data x,y, and with num_knots number of knots.

    :param X: x data
    :param Y: y data
    :param num_knots: Number of knotw/breakpoints (three is normal for fractal analysis)
    :return: returns a tuple with the x and y values for the each slope
    """
    xmin = X.min()
    xmax = X.max()

    seg = np.full(num_knots - 1, (xmax - xmin) / num_knots)

    px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
    py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.01].mean() for x in px_init])

    def func(p):
        seg = p[:num_knots - 1]
        py = p[num_knots - 1:]
        px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        return px, py

    def err(p):
        px, py = func(p)
        Y2 = np.interp(X, px, py)
        return np.mean((Y - Y2)**2)

    r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    return func(r.x)


def spline_interpolate(divider_width, polygon_length, num_knots):
    """
    spline_interpolate Spline interpolation of the input

    :param divider_width: List of divider widths
    :param poly_length:   List of polygon lengths
    :param num_knots:     Integer with number of knots for the spline interpolation
    :return:              Returns tuple with spline interpolated x and y
    """

    xData = np.log(np.array (divider_width))
    yData = np.log(np.array(polygon_length))
    knot_numbers = num_knots                         # NOTE: Cannot have more knots than the length of xData
    x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
    q_knots = np.quantile(xData, x_new)
    t,c,k = interpolate.splrep(xData, yData, t=q_knots, s=1)
    # Make even spacing between datapoints on x axis
    new_x = np.linspace(min(xData), max(xData), 100)
    new_y = interpolate.BSpline(t,c,k)(new_x)
    return new_x, new_y


def fractal_analysis(analysis_signal, divider_range, window_size, skipping_step):
    """
    Fractal analysis with sliding window method

    :param analysis_signal: Input signal to analyse
    :param divider_range: Range of divider widths
    :param window_size: Size of windows
    :param skipping_step: Steps to skip between each window (1 signifies no skipping)
    :return: Returns a tuple with index of calculated fractal dimension, middle-wave- and long-waved fractal dimension.
    """
    x = [i for i in range(len(analysis_signal))]

    results = []                          # List with tuples, of position(x), DR1 and DR2
    progress_bar = tqdm(total=(len(analysis_signal)-window_size))

    for i in range(0,len(analysis_signal)-1,skipping_step):
        if i==0:
            old_pos = 0
        else:
            old_pos = pos

        pos = x[i]
        progress_bar.update(pos-old_pos)
        progress_bar.set_description(f'Fractal Analysis Distance Covered')
        if pos+window_size>=len(analysis_signal): break # Stop when start position is at end distance
        end = pos + window_size    # Local section should be pos:pos+window_size
        end_index = end-1#0
        temp = 0

        data = analysis_signal[i:end_index]
        x_values = x[i:end_index]
        data.reset_index(inplace=True, drop=True)

        # Run local analysis
        analysis = [calculate_total_length(x_values, data, step)  for step in divider_range] #List of fractal analysis
        divider_length, poly_length, reg_param = zip(*analysis)


        # Remove outliers with Spline interpolation
        new_x, new_y = spline_interpolate(divider_length, poly_length, 12)

        # Piecewise Regression
        px, py = segments_fit(new_x, new_y, 3)

        # Calculate fractal dimension (Dr = 1 - m)
        m1 = (py[3] - py[2])/(px[3] - px[2])
        m2 = (py[2] - py[1])/(px[2] - px[1])
        Dr1 = 1 - m1 # Measure on 1. order roughness
        Dr2 = 1 - m2 # Measure on 2. order roughness


        results.append((pos,Dr1, Dr2))
    progress_bar.close()

    return results
