import gpxpy
import geopy.distance
import math
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson
import numpy as np

#----------USER VALUES HERE----------------

SAMPLE_PATH = "Samples/SampleRecording.gpx"
ROUTE_PATH = "Routes/SnowslideLakeTrail.gpx"

#-------ACTUAL CODE STARTS BELOW----------------

def getValues(cur_point, prev_point):
    distance_x = geopy.distance.geodesic((prev_point.latitude, prev_point.longitude), (cur_point.latitude, cur_point.longitude)).m

    # Calculating the distance between two coordinates is pretty complex because the earth isn't flat, so we use this library.
    distance_x = geopy.distance.geodesic((prev_point.latitude, prev_point.longitude), (cur_point.latitude, cur_point.longitude)).m

    if not distance_x:
        return None, None
    #print("Distance in meters: ", distance_x)

    # Calculating change in height
    distance_y = cur_point.elevation - prev_point.elevation
                
    # Calculate angle
        
    angle = math.degrees(math.atan(distance_y / distance_x))
    #print("Angle: ", angle)
                
    # Pythagorean Theorem for total distance
    total_distance = (distance_x**2 + distance_y**2)**0.5


    return angle, total_distance

def chunks(lst, n):
    return [lst[idx: idx + n] for idx in range(0, len(lst), n)]

def buildTimePlot(filename, speed_curve, bin_size, values = False):
    """Builds a plot of trail grades at different times. Referred to as G(t) in the paper."""
    
    gpx_file = open(filename, 'r')
    gpx = gpxpy.parse(gpx_file)
    elapsed = 0
    num_points = 0

    time_plot = []
    bin_group = []
    
    for track in gpx.tracks:
        for segment in track.segments:
            prev_piece = None
            prev = None
            total_distance = 0
            for piece in chunks(segment.points, bin_size):
                for point in piece:
                    if not prev:
                        prev = point
                        continue
                    
                    angle, distance = getValues(point, prev)

                    if not angle:
                        continue

                    total_distance += distance

                    prev = point

                if not prev_piece:
                    prev_piece = [piece[0]]

                percent_grade = (piece[0].elevation - prev_piece[0].elevation)/ total_distance
                time_plot.append((elapsed, percent_grade))
                speed = speed_curve(angle)
                elapsed += (total_distance / speed) # Could this be more accurate with integration?? Figure this out later if time permits

                prev_piece = piece

                total_distance = 0
                
    x_list, y_list = zip(*time_plot)
    time_spline = CubicSpline(x_list, y_list)

    if values:
        return x_list, y_list

    return time_spline, (0, x_list[-1])

def buildTimeScatter(filename, speed_curve):
    """Builds a plot of trail grades at different times. Referred to as G(t) in the paper."""
    
    gpx_file = open(filename, 'r')
    gpx = gpxpy.parse(gpx_file)
    elapsed = 0
    num_points = 0

    time_plot = []
    bin_group = []
    
    for track in gpx.tracks:
        for segment in track.segments:
            prev = None
            for point in segment.points:
                if not prev:
                    prev = point
                    continue
                    
                angle, distance = getValues(point, prev)

                if not angle:
                    continue
                percent_grade = np.arctan(np.radians(angle))

                speed = speed_curve(angle)
                elapsed += (distance / speed)

                time_plot.append((elapsed, percent_grade))
                prev = point

    x_list, y_list = zip(*time_plot)
    return x_list, y_list

def analyzeRecording(filename, num_bins):
    """Returns steepness-speed cubic spline curve for a given GPX recording. Speed is in meters per minute and steepness is in degrees."""
    gpx_file = open(filename, 'r')
    gpx = gpxpy.parse(gpx_file)
    bin_size = 180 / num_bins
    
    bins = [[0] for x in range(num_bins)]
    
    for track in gpx.tracks:
        for segment in track.segments:
            prev = None
            for point in segment.points:
                if not prev:
                    prev = point
                    continue

                angle, total_distance = getValues(point, prev)

                if not total_distance:
                    continue

                # Code for analyzing timestamps
                # Calculate time difference
                elapsed = (point.time - prev.time).seconds
                    

                # Calculate speed in m/min
                speed = 60 * total_distance / elapsed
                    
                # Determine bin designation
                angle_bin_index = int((angle + 90) // bin_size)
                bins[angle_bin_index].append(speed)

                # Update previous pos
                    
                prev = point

    x_list = [x * bin_size - 90 for x in range(num_bins)]
    y_list = [sum(y)/len(y) for y in bins]
    curve = [(x_list[i], y_list[i]) for i in range(len(x_list))]

    # Some post-processing to make the curve more usable
    curve.sort(key = (lambda pair: pair[0]))

    
    # Sweep thru and smooth out zero values
    midpoint = num_bins // 2
    # Positive sweep
    for count, val in enumerate(curve[midpoint:-1]):
        idx = midpoint + count
        if not val[1]:
            #print((curve[idx][0],(curve[idx-1][1]+curve[idx+1][1]) / 2))
            curve[idx] = (curve[idx][0],(curve[idx-1][1]+curve[idx+1][1]) / 2)
            
    
    # Negative sweep
    for count, val in enumerate(curve[midpoint:1:-1]):
        idx = midpoint - count
        if not val[1]:
            curve[idx] = (curve[idx][0],(curve[idx-1][1]+curve[idx+1][1]) / 2)
    
    xs_list, ys_list = zip(*curve)
    # Approximation using cubic splines
    speed_func = CubicSpline(xs_list, ys_list)

    return speed_func


def displayCurve(filename):
    myCurve = analyzeRecording(filename, 30)
    matplotlib.pyplot.scatter(*zip(*myCurve))
    matplotlib.pyplot.show()

def getVO2(time, steepness_curve, speed_curve):
    """Returns VO2 at a specific time"""
    # Calculate grade function, convert steepness curve to grade
    G = steepness_curve(time)
    S = speed_curve(np.degrees(np.arctan(G)))
    D = S * 0.1 * 0.73 + 3.5 # Multiplier of .73 for downhill grades
    return ((S * 0.1) + (S * G * 1.8)) * (G >= 0) + D * (G < 0) + 3.5

def variation(curve, time, window):
    return curve(time) / curve(time + window) - 1

def steadyStateCrossings(VO2_curve, x_vals, window):
    crossings = []
    v = variation(VO2_curve, x_vals, window)
    for idx in range(len(v) - 1):
        val1 = abs(v[idx])
        val2 = abs(v[idx + 1])
        if val1 > .1 and val2 < .1:
            crossings.append([idx])
        if val1 < .1 and val2 > .1:
            crossings[-1].append(idx)

    crossings[-1].append(len(x_vals) - 1)    

    return crossings

def getCalories(sample_name, route_name, weight, analysis_num_bins = 40, steepness_bin_size = 100, window = 4):

    steepness_speed_curve = analyzeRecording(sample_name, analysis_num_bins)
    time_steepness_curve, time_domain = buildTimePlot(route_name, steepness_speed_curve, steepness_bin_size)
    x_vals = np.arange(*time_domain, .1)

    VO2 = getVO2(x_vals, time_steepness_curve, steepness_speed_curve)
    VO2_curve = CubicSpline(x_vals, VO2)
    steady_state_VO2 = 0

    crossings = steadyStateCrossings(VO2_curve, x_vals, window)

    for interval in crossings:
        start = interval[0]
        stop = interval[1]
        steady_state_VO2 += simpson(VO2[start:stop], x_vals[start:stop])

    non_steady_state_VO2 = 0

    start = 0
    for idx in range(len(crossings)):
        interval = crossings[idx]

        if idx == len(crossings) - 1:
            stop = len(x_vals) - 1
        else:
            stop = interval[0]
            
        non_steady_state_VO2 += (x_vals[stop] - x_vals[start]) * ((VO2[stop] - VO2[start]) / 2 + VO2[start])
        start = interval[1]
    
    total_VO2 = steady_state_VO2 + non_steady_state_VO2

    kcal = 5 * weight * total_VO2 / 1000

    return kcal

if __name__ == "__main__":
    
    print("Calculating total calories...")
    print("Total calculated energy expenditure: ", getCalories(SAMPLE_PATH, ROUTE_PATH, 72), " Calories.")
    print("Generating graphs...")

    plt.subplot(2, 2, 1)
    
    steepness_curve = analyzeRecording(SAMPLE_PATH, 40)
    degrees = np.arange(-90, 90, 180 // 40)
    degrees_fine = np.arange(-90, 90, .1)
    plt.scatter(degrees, steepness_curve(degrees))
    plt.plot(degrees_fine, steepness_curve(degrees_fine))

    plt.xlabel("Angle (deg)")
    plt.ylabel("Speed (m/min)")
    plt.title("Speed in meters per minute for a given angle")
    
    plt.subplot(2,2,2)

    plt.scatter(*buildTimeScatter(ROUTE_PATH, steepness_curve))

    plt.xlabel("Time (min)")
    plt.ylabel("Steepness (frac)")
    plt.title("Noisy data: predicted fractional grades at time t")

    plt.subplot(2,2,3)

    spline_y, time_domain = buildTimePlot(ROUTE_PATH, steepness_curve, 100)
    x_vals = np.arange(*time_domain, .1)
    plt.plot(x_vals, spline_y(x_vals))

    plt.xlabel("Time (min")
    plt.ylabel("Steepness (frac)")
    plt.title("Smoothed data: cubic spline interpolation of average slope")

    plt.subplot(2,2,4)

    VO2 = getVO2(x_vals, spline_y, steepness_curve)
    VO2_curve = CubicSpline(x_vals, VO2)
    plt.plot(x_vals, VO2)

    plt.xlabel("Time (min)")
    plt.ylabel("VO2 (ml/(kg*min))")
    plt.title("Steady state VO2 model")

    print("Done")
    plt.show()
