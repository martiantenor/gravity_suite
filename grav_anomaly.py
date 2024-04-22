#!/usr/bin/env python
# A program to plot the gravity anomaly of elements from an Abaqus FEM
#
# Contact Dave Blair (dblair@purdue.edu) with questions
#
# (c) David Blair, 2015. This work is licensed under a Creative Commons
# Attribution-NonCommercial-ShareAlike Unported License
# (http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US)


######## Preface ###############################################################

from __future__ import division
from math import pi, sqrt, radians, degrees, sin, cos
import sys, string, argparse, re
import numpy

# We need elliptic integrals. First we'll look for SciPy's internal functions,
# which are a little more accurate.
try:
    from scipy.special import ellipk as scipy_ellipk
    from scipy.special import ellipe as scipy_ellipe
    scipy_ellip_mode = True

# If SciPy isn't installed, we'll instead use the elliptic integral from
# 'special_functions.py', which is based on work that is Copyright (c) 1996 John
# Wiley & Sons, Inc; see `special_functions.py` for more details.
except ImportError:
    print "WARNING: Using Python elliptic integrals (less accurate and slower)"
    from special_functions import python_comelp as python_ellip
    scipy_ellip_mode = False


__version__ = "3.0"
# For ToDo and Changelog, please refer to "grav_anomaly.py.changelog"



######## Options ###############################################################
"""
Set default values here. These can alternatively be set (or re-set) on the
command line using the flags described at the end of this program
"""
# Constants and calculation parameters:

elevation = 10000   # Elevation of the spacecraft (gets added to all depths) [m]
num_steps = 101     # Higher values --> smoother graph, longer runtime
G = 6.6738E-11      # Gravitational constant [m3.kg-1.s-2]

curved_mode = False

# Flat model mode options: range of x values to look at
planet_radius = 1740e3  # m (Moon, to 3 sig figs)
x_range = 2000e3        # m, maximum r of "spacecraft" [for flat models]
theta_range = 60        # degrees, maximum theta of "spacecraft" [for curved models]


# Calculation behavior:

# Compare against an external baseline file, or just plot the acceleration due
# to gravity without comparing to anything?
anomaly_mode = False

# Subtract off topography effects after doing the anomaly calculation, to arrive
# at a "Bouguer" result? (Leaving this as "False" allows the command-line flag
# "--geoidname FILE" to enter free-air calculation mode)
bouguer_mode = False
bouguer_density_cutoff = 100.0                  # kg.m-3
bouguer_crustal_density_assumption = 2550.0     # kg.m-3

# Do you want to subdivide elements? If so, set "number_of_subdivisions" to the
# number of times you want to subdivide. The default is here to allow typing
# "-s" instead of, e.g., "-s 1" if this option is chosen via the command line.
number_of_subdivisions = 1
default_number_of_subdivisions = 1


# Input files:

# The baseline file to use for anomaly calculations. This is the default
# filename, and it can be changed on the command line with "--geoidname [name]"
# at runtime
geoidfile_name = "foo_geoid_filename.xy"


# Output options:

# Use the GMT mode?
GMT_mode = True

# What is the file extension of the gravity/data file(s)? This is used when
# generating names for the GMT files.
gravfile_suffixes = (".grav",".xy")

# What are the filename suffixes for the GMT output files? They're specified as
# ((free-air acceleration, free-air anomaly), (Bouguer accel., Bouguer anomaly),
# (raw acceleration))
GMT_suffixes = (("_fa_acc.xy","_fa.xy"),("_boug_acc.xy","_boug.xy"),("_acc.xy"))

# Print headers at the top of GMT files?
GMT_print_headers = True

# Verbosity options. The default (neither Verbose nor Quiet) is a "medium"
# level of output, showing files being processed and similar information.
# Verbose mode adds the "progress bar" output during computation and a few other
# things. Quiet mode supresses everything (except errors).
verbose_mode = True
quiet_mode = False



######## Main Program ##########################################################

# Section 1: Elements
# A class for Elements, and a few methods for calculating things about them and
# splitting up sets of them

class Element:
    """
    A simple "Element" class, for storing data and calculating mass
    """

    def __init__(self, coords, this_density, nodal_coordinates, volume="na"):

        # Reported values for the element as a whole: coordinates of the
        # centroid and density
        self.centroid = coords
        self.density = float(this_density)

        #self.y_coord = coords[1]       # y coord of the centroid

        # Info about the nodes of the element
        self.number_of_nodes = len(nodal_coordinates)
        self.nodes = []
        for node in nodal_coordinates:
            self.nodes.append([node[0],node[1]])

        # Calculate the area of the element. If it's a triangle, just do the
        # math; if it's a quadrilateral, break into two triangles and add
        if self.number_of_nodes == 3:
            self.area = find_triangle_area(self.nodes[0],
                                           self.nodes[1],
                                           self.nodes[2])
        if self.number_of_nodes == 4:
            self.area = sum( (find_triangle_area(self.nodes[0],
                                                 self.nodes[1],
                                                 self.nodes[2]),
                              find_triangle_area(self.nodes[0],
                                                 self.nodes[2],
                                                 self.nodes[3])) )

        # If volume isn't specified in the Element definition, then calculate
        # it, otherwise use the specified value
        if volume == "na":
            self.volume = self.area * 2 * pi * self.centroid[0]
        else:
            self.volume = volume

        # Finally, calculate mass, depth (y coord represented as a positive
        # number), and z (total distance between element and spacecraft)
        self.mass = self.volume * self.density
        self.depth = -self.centroid[1]
        #self.z = self.depth + elevation

    def plot(self, color='r'):
        "Makes a simple picture of the element, using plt"

        x_coords = []
        y_coords = []

        # Iterate over nodes so that the plotting is independent of the number
        # of nodes being fed into this function
        for i in range(len(self.nodes)):
            x, y = self.nodes[i][0], self.nodes[i][1]
            x_coords.append(x), y_coords.append(y)

        plt.fill(x_coords, y_coords, color)
        #plt.show()


def find_triangle_area(A,B,C):
    """
    Find the area of a triangle in 2D from the coordinates of its nodes.
    """
    return abs(.5* (A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1])) )


def distance_formula(coords1,coords2):
    """
    Get the distance between two (x,y) coordinate pairs
    """
    return sqrt( (coords2[0]-coords1[0])**2 + (coords2[1]-coords1[1])**2 )


def subdivide_element(element):
    """
    This subroutine takes an element and splits it into two. It can interpret
    whether that element has three or four nodes, and work accordingly; the
    resulting elements will always be triangles.
    """
    number_of_nodes = element.number_of_nodes
    density = element.density

    def spatial_average(coordinate_pairs):
        "Get the centroid of a given set of elements"
        x_coords = []
        y_coords = []

        for pair in coordinate_pairs:
            x_coords.append(pair[0])
            y_coords.append(pair[1])
        midpoint_x = sum(x_coords)/len(x_coords)
        midpoint_y = sum(y_coords)/len(y_coords)

        return midpoint_x, midpoint_y

    # Depending on how many nodes our element has, we need to implement
    # different routines for subdividing it.

    # For a quadrilateral, we can just do this arbitrarily, splitting it into
    # two elements by counting the nodes from the upper-right (or really, from
    # any given starting node)
    if number_of_nodes == 4:
        newnodes1 = [element.nodes[0], element.nodes[1], element.nodes[2]]
        newnodes2 = [element.nodes[0], element.nodes[2], element.nodes[3]]

    # For a triangle, we'll find the longest side, get the midpoint of that
    # side, and use that to divide the original triangle into two smaller
    # triangles. This avoids getting thinner and thinner slices, and should cut
    # down on the inaccuracy a little bit
    if number_of_nodes == 3:

        # Find the lengths of the sides of the element, and extract the length
        # of the longest side
        AB = distance_formula(element.nodes[0],element.nodes[1])
        BC = distance_formula(element.nodes[1],element.nodes[2])
        CA = distance_formula(element.nodes[2],element.nodes[0])
        longest_side = max((AB,BC,CA))

        # Go through the options, and if a side matches that longest side, stop
        # there. This means that in the case of an equilateral or isosceles
        # triangle, it's arbitrarily picking the first of the two equal sides
        # that it comes to, which is fine. There's probably a faster way to do
        # this with some sort of loop, but I haven't run the tests; at least
        # this way, it'll exit as soon as it does one of these three.
        if AB == longest_side:
            newvertex = spatial_average((element.nodes[0],element.nodes[1]))
            newnodes1 = [element.nodes[0], newvertex, element.nodes[2]]
            newnodes2 = [element.nodes[2], newvertex, element.nodes[1]]
        elif BC == longest_side:
            newvertex = spatial_average((element.nodes[1],element.nodes[2]))
            newnodes1 = [element.nodes[1], newvertex, element.nodes[0]]
            newnodes2 = [element.nodes[0], newvertex, element.nodes[2]]
        elif CA == longest_side:
            newvertex = spatial_average((element.nodes[2],element.nodes[0]))
            newnodes1 = [element.nodes[2], newvertex, element.nodes[1]]
            newnodes2 = [element.nodes[1], newvertex, element.nodes[0]]

    # Now we're back into behavior that doesn't depend on the number of nodes,
    # since either case results in our defining two new triangles.

    # Figure out the centroids from the new sets of nodes
    centroid1 = spatial_average(newnodes1)
    centroid2 = spatial_average(newnodes2)

    # Determine the volume of the new elements
    #volume1 = find_triangle_area(newnodes1[0], newnodes1[1], newnodes1[2]) \
    #          * 2 * pi * centroid1[0]
    #volume2 = find_triangle_area(newnodes2[0], newnodes2[1], newnodes2[2]) \
    #          * 2 * pi * centroid2[0]

    # Use this information to define and return two new elements.
    newelement1 = Element(centroid1,density,newnodes1)
    newelement2 = Element(centroid2,density,newnodes2)
    #newelement1 = Element(centroid1,density,newnodes1,volume=volume1)
    #newelement2 = Element(centroid2,density,newnodes2,volume=volume1)

    return newelement1, newelement2


def subdivide_element_set(starting_elements, number_of_subdivisions):
    """
    Implement the subdivision of elements - recursively iterating over the set
    of starting_elements and creating new final_elements sets as we go, until we
    reach the number of subdivisions we want
    """

    # Initially, we assume that the final element set is equal to the starting
    # element set - this way, if there are 0 subdivisions to make, then we can
    # return the final set and get the same result as if we'd done nothing
    final_elements = starting_elements[:]

    # Subdivide the elements n times (resulting in 2^n times as many elements).
    # We'll start the counter n at zero, and iterate until it's equal to
    # number_of_subdivions
    n = 0

    while n < number_of_subdivisions:

        # Slice the final elements into a new set that we can work on while
        # we're in this loop (the "input" of this subdivision). After doing
        # this, erase the final_elements set so that we can add to it as we
        # go through the input elements and generate new ones
        input_elements = final_elements[:]
        final_elements = []

        # Go through every element that was fed into this iteration, and split
        # it into two elements. Add these to the "output" from this round
        for this_element in input_elements:
            new_element1, new_element2 = subdivide_element(this_element)
            final_elements.append(new_element1)
            final_elements.append(new_element2)

        n += 1

    return final_elements


# Section 2: Math
# Elliptic integrals of the 1st and 2nd kind, gravity over a ring, functions to
# calculate gravity over an element, and a function to look at the difference in
# acceleration between two sets of data


def K(m):
    """
    The complete elliptic integral of the first kind, such that

        K(m) = K(pi/2,m)
             = integral(1 / sqrt(1 - m*sin(alpha)**2), {alpha, 0, pi/2})

    as given in Abramowitz & Stegun (1972; Eq. 17.3.3), and in Mathematica and
    SciPy.

    The non-SciPy Python code, however, uses a slightly different form, 
    Ryzhik (1965), which is also the one used in Turtle and Melosh (1998). That
    form uses "k" instead of "m":

        K(k) = E(pi/2, k)
             = integral(sqrt(1 - (k**2)*(sin(alpha)**2)), {alpha, 0, pi/2})

    This function accounts for the difference automatically, using "m" as the
    input for the SciPy version of the routine, and "k" for the regular Python
    version (which is easy, since k = sqrt(m)).
    """

    if scipy_ellip_mode:
        return scipy_ellipk(m)
    else:
        # python_ellip returns (K(k), E(k))
        return python_ellip(sqrt(m))[0]


def E(m):
    """
    The complete elliptic integral of the second kind, such that

        E(m) = E(pi/2, m)
             = integral(sqrt(1 - m*sin(alpha)**2), {alpha, 0, pi/2})

    as given in Abramowitz & Stegun (1972; Eq. 17.3.3), and in Mathematica and
    SciPy.

    The non-SciPy Python code, however, uses the form found in Gradshteyn and
    Ryzhik (1965), which is also the one used in Turtle and Melosh (1998). This
    form uses "k" instead of "m":

        E(k) = E(pi/2, k)
             = integral(sqrt(1 - (k**2)*(sin(alpha)**2)), {alpha, 0, pi/2})

    This function accounts for the difference automatically, using "m" as the
    input for the SciPy version of the routine, and "k" for the regular Python
    version (which is easy, since k = sqrt(m)).
    """

    if scipy_ellip_mode:
        return scipy_ellipe(m)
    else:
        # python_ellip returns (K(k), E(k))
        return python_ellip(sqrt(m))[1]


def ring_gravity_vertical(radius_of_ring, position_above_ring, radius_of_observer,
                 ring_mass):
    """
    Gives the vertical component of gravity above a ring at a given elevation &
    distance from the axis of symmetry

    Based on E.P. Turtle and H.J. Melosh, Stress and Flexural Modeling of the
    Martian Lithospheric Response to Alba Patera, 1987.
    """

    # Define our variables
    a = radius_of_ring
    M = ring_mass
    z = position_above_ring
    R = radius_of_observer

    # Now calculate the actual gravity at the desired point
    Az = ((2*G*M*z)/(pi * ((R-a)**2 + z**2) * ((R+a)**2 + z**2)**(1/2)))    \
          * E((4*a*R)/((R+a)**2 + z**2))

    #return the vertical acceleration
    return Az


def ring_gravity_radial(radius_of_ring, position_above_ring, radius_of_observer,
                 ring_mass):
    """
    Gives the radial component of gravity above a ring at a given elevation &
    distance from the axis of symmetry

    Based on E.P. Turtle and H.J. Melosh, Stress and Flexural Modeling of the
    Martian Lithospheric Response to Alba Patera, 1987.
    """

    # Define our variables
    a = radius_of_ring
    M = ring_mass
    z = position_above_ring
    R = radius_of_observer

    # Now calculate the actual gravity at the desired point
    Ar = ((2*G*M)/(pi * ((R+a)**2 + z**2)**(1/2)))      \
         * (    \
            (R*E((4*a*R)/((R+a)**2 + z**2)) / ((R-a)**2 + z**2) )    \
           +(K((4*a*R)/((R+a)**2 + z**2)) / (2*R) )    \
           +(E((4*a*R)/((R+a)**2 + z**2)) / (2*R*(((2*a*R)/(R**2+a**2+z**2))-1)) )    \
           )

    #return the radial acceleration
    return Ar


def grav_calculator(element_set):
    """
    Sums up the gravity of a set of elements and returns the result as a
    "datapair" (two matched lists: x coordinates, and gravitational acceleration
    at those coordinates)
    """

    # Initialize variables
    if not curved_mode:
        x_values = numpy.linspace(0,x_range,num=num_steps)
    else:
        theta_values = numpy.linspace(0,radians(theta_range),num=num_steps)
        spacecraft_radius = planet_radius + elevation
    grav_anomalies = []

    # Go through each "spacecraft" position and figure out Az contribution from
    # every element
    current_step = 1

    # For non-curved models:
    if not curved_mode:
        for R in x_values:

            # Start the value of the cell at zero
            anomaly_here = 0.0

            # Go through each element and add its Az to anomaly_here
            for this_element in element_set:
                #R is already defined by virtue of our being in a loop
                a = this_element.centroid[0]
                z = this_element.depth + elevation
                M = this_element.mass
                anomaly_here += ring_gravity_vertical(a,z,R,M)

            grav_anomalies.append(anomaly_here)

            # Print progress information if we're in verbose mode. The labeling of
            # which iteration we're in goes back to code that's in the command-line
            # implementation part of this program, to tell us where in the loop we
            # are.
            if verbose_mode:
                # Three options for what to display while running (in verbose mode):
                # Leave the most desirable one uncommented

                # Display x position and value
                print "At lateral position %.0f / %.0f (pass %i of %i): %.4e mGal"\
                       %(R, x_values[-1],
                         this_iteration+1, total_iterations,
                         anomaly_here*1e5)

                # Or display x position
                #print "At lateral position %.0f / %.0f (pass %i of %i)"\
                #       %(R,x_values[-1], this_iteration+1, total_iterations)

                # Or display which step we're on
                #print "At lateral step %3i / %3i (round %2i of %2i)"\
                #       %(current_step,num_steps,this_iteration+1,total_iterations)

            current_step = current_step + 1

    # For curved models:
    else:
        for Th in theta_values:

            # Start the value of the cell at zero
            anomaly_here = 0.0

            # Go through each element and add its Az to anomaly_here
            for this_element in element_set:

                # Distance of spacecraft from the ring, in Cartesian coords
                a = this_element.centroid[0]
                z = spacecraft_radius*cos(Th) - this_element.centroid[1]
                R = spacecraft_radius*sin(Th)
                M = this_element.mass

                # At R=0, the radial component is always zero, so we'll set that
                # here (this also avoids a "divide by zero" error)
                if R == 0:
                    anomaly_here += ring_gravity_vertical(a,z,R,M)


                # If we're not at R=0, calculate the R component of the gravity.
                # To get this, we transform into a new x'-y' coordinate system,
                # and then take the y' components of both the radial and
                # vertical gravity
                else:
                    Az = ring_gravity_vertical(a,z,R,M)
                    Ar = ring_gravity_radial(a,z,R,M)

                    anomaly_here += Ar*sin(Th) + Az*cos(Th)

            grav_anomalies.append(anomaly_here)

            # Print progress information if we're in verbose mode. The labeling
            # of which iteration we're in comes from code in the command-line
            # implementation part of this program which tells us where in the
            # loop we are.
            if verbose_mode:
                ## Display current theta position
                #print "At theta position %.0f / %.0f (pass %i of %i): %.4e mGal"\
                #       %(degrees(Th), degrees(theta_values[-1]),
                #         this_iteration+1, total_iterations,
                #         anomaly_here*1e5)

                # Or display which step we're on:
                print "At step %3i / %3i (pass %i of %i): %.4e mGal"\
                      %(current_step,num_steps,
                        this_iteration+1, total_iterations,
                        anomaly_here*1e5)

            current_step = current_step + 1

    # Multiply all the y values by 1e5 to get the answer in mGal
    grav_values_scaled = []
    for this_item in grav_anomalies:
        mGal_multiplication_factor = 1e5
        grav_values_scaled.append(this_item * mGal_multiplication_factor)

    # Format the output so that it's more legible:
    if not curved_mode:
        # Divide x values by 1000 to go from m to km
        output_values = []
        for this_item in x_values:
            output_values.append(this_item/1000.)
    else:
        # Convert theta values (radians from pole) into km from the pole at the
        # radius of the planet's surface
        output_values = []
        planet_circumference = (2*planet_radius*pi) / 1000. # in km
        for this_item in theta_values:
            output_values.append( planet_circumference * (this_item/(2*pi)))

    # Output data
    return output_values, grav_values_scaled


def delta_calculator(initial_datapair, final_datapair):
    """
    Calculates the gravity difference between two datapairs, which is the
    delta/anomaly if you're using a geoid-type baseline file for "initial" and a
    model for "final"
    """

    # Go through both sets simultaneously, and average both their gravity data
    # and the x position of the centroid (which may be different, since we're
    # looking at initial vs. final model state!)
    delta_datapair = [[],[]]
    for i in range(len(initial_datapair[0])):
        initial_x = initial_datapair[0][i]
        final_x = final_datapair[0][i]
        initial_grav = initial_datapair[1][i]
        final_grav = final_datapair[1][i]

        average_x = (initial_x + final_x)/2.
        delta_grav = final_grav - initial_grav

        delta_datapair[0].append(average_x)
        delta_datapair[1].append(delta_grav)

    return delta_datapair


# Section 3: I/O
# Functions for determining file type, parsers for different types of files, and
# functions to output data to either text files or the plt plotter


def xyfile_parser(datafile):
    """
    Read the data from an xy-file and return scaled x values and gravity values.
    This function bypasses all of the gravity calculations, allowing us to
    compare files with just xy data, as might be generated by another program,
    or in a previous run of this program with GMT output.

    WARNING: make sure that the .xy files have the first column as distance in
    km, and the second column as *RAW ACCELERATION* in mGal. The file also has
    to be the same number of lines as num_steps for any sort of comparison (e.g.
    an anomaly calculation) to work. As long as you've kept that parameter the
    same, though, this can be a real time-saver (e.g. not having to recalculate
    the geoid file every time)
    """

    # Read the file, store the lines, then close the file
    these_lines = datafile.readlines()
    datafile.close()

    x_values = []
    y_values = []
    for this_line in these_lines:

        # Ignore comment lines
        if this_line.startswith("#"):
            continue

        # Read in the data
        x_value = float(this_line.split()[0])
        y_value = float(this_line.split()[1])

        if x_value <= (x_range/1000):
            x_values.append(x_value)
            y_values.append(y_value)
        else:
            pass

    return x_values, y_values


def gravfile_parser(datafile, do_bouguer_subtraction):
    """
    Read the data from the *.grav file that's output by grav_plugin.py from
    an open .odb in Abaqus/CAE, and save that data as an 'element_set'
    """
    # Read the file, close it, and then cut it down to ignore the first 3 lines
    all_lines = datafile.readlines()
    datafile.close()
    data_lines = all_lines[3:]

    # Make doubly sure the user knows what's going on if we're in Bouguer mode
    if verbose_mode and do_bouguer_subtraction:
        print "Doing Bouguer corrections..."

    # Iterate over the lines in the data_lines set, and extract the data in the
    # form of Element objects, which get added to a set and finally returned out
    # of the subroutine
    element_set = []
    for this_line in data_lines:
        element_data = this_line.split("|")[0].split()
        this_x = float(element_data[1])
        this_y = float(element_data[2])
        this_density = float(element_data[3])

        # Here's where we do the Bouguer correction, under two conditions:
        #   - if it's above 0 and has "positive" density (> 100 kg.m-3), make
        #     that mass 0 kg.m-3 instead (roughly negating positive topography)
        #   - if it's below zero and has "0" density (< 100 kg.m-3), give it a
        #     density of 2550 kg.m-3 instead (roughly negating negative
        #     topography)
        #   - note that, for curved models, "above zero" means "total distance
        #     from center of planet greater than planet radius"
        if do_bouguer_subtraction:

            if not curved_mode:
                if (this_y > 0.0) and (this_density > bouguer_density_cutoff):
                    this_density = 0.0
                elif (this_y < 0.0) and (this_density < bouguer_density_cutoff):
                    this_density = bouguer_crustal_density_assumption

            if curved_mode:
                this_r = sqrt(this_x**2 + this_y**2)
                if (this_r > planet_radius) and (this_density > bouguer_density_cutoff):
                    this_density = 0.0
                elif (this_r < planet_radius) and (this_density < bouguer_density_cutoff):
                    this_density = bouguer_crustal_density_assumption


        node_data = eval(this_line.split("|")[-1].strip())
        element_set.append(Element( (this_x, this_y),
                                    this_density,
                                    node_data ))

    return element_set


def get_data_from_file(this_file):
    """
    Given a file object as input, determine if that file is a .grav file or a
    .xy file. If it's neither, raise an error and quit.
    """

    # Pull in the global variables for keeping track of the total number of
    # iterations and what number this iteration happens to be, for the purposes
    # of outputting a "progress bar" for the user
    global total_iterations
    global this_iteration

    # Assuming we're not supressing output, let the user know what file we're
    # working on
    if not quiet_mode:
        print "Parsing file %s..."%this_file.name

    # Check to see what type of file we have in an if/elif/elif block, and
    # implement the appropriate code. We'll determine filetype by looking at the
    # letters after the last "." in a case-insensitive way

    if this_file.name.split(".")[-1].upper() == "XY":

        # Grab the base name of the file
        this_data_basename = this_file.name.strip(".xy").strip(".XY")

        # Make sure this isn't a Bouguer calculation, which currently can't be
        # done on XY files (unless it's the geoid file, then it's ok)
        if this_file.name != geoidfile_name and bouguer_mode:
            print bouguer_on_nongravfile_errortext
            sys.exit()

        # Parse the file
        this_data = xyfile_parser(this_file)

        # Subtract one from "total_iterations" so that it doesn't show up in
        # the "round X of Y" display (since there's no processing going on)
        total_iterations = total_iterations - 1

    elif this_file.name.split(".")[-1].upper() == "GRAV":

        # Grab the base name of the file
        this_data_basename = this_file.name.strip(".grav").strip(".GRAV")

        # Parse the *.grav file and get information about all of its elements.
        # Make sure that we're not processing the geoid file in Bouguer mode
        # (even though, if the ground level is at y=0, this should do nothing)
        if this_file.name != geoidfile_name:
            starting_elements = gravfile_parser(this_file,bouguer_mode)
        else:
            starting_elements = gravfile_parser(this_file,False)

        # Implement the subdivider routine on those elements
        final_elements = subdivide_element_set(starting_elements,
                                               number_of_subdivisions)

        # Print out a helpful notice about how many elements we start with
        # and how many we end up with; this is important enough to include
        # it even in non-verbose mode (but not in quiet mode)
        if not quiet_mode:
            print "Original number of elements in %s: %i"%(this_file.name,
                                                           len(starting_elements))
            print "Number of elements in %s after subdivision: %i"%(this_file.name,
                                                                    len(final_elements))

        # Process the resulting element set with grav_calculator()
        this_data = grav_calculator(final_elements)

        # Up the iteration count, since we did some actual math
        this_iteration = this_iteration + 1


    # If it's not a file type we recognize, do this:
    else:
        print "Error: Please choose a *.grav or *.xy file"
        parser.print_usage()
        sys.exit(1)

    return this_data, this_data_basename


def GMT_outputter(data_sets, outfiles, outfile_names):
    """
    Creates text files of the output for plotting in GMT
    """

    # Grab just the data from data_sets, and make an iterable set
    data_pairs = []
    for key in data_sets.keys():
        data_pairs.append(data_sets[key])

    # Output the filenames if we're not in quiet mode (doesn't require verbose)
    if not quiet_mode:
        printstring = "Creating (or overwriting!) %s"%outfile_names[0]
        if len(outfile_names) == 1:
            print printstring
        if len(outfile_names) == 2:
            print printstring + " and %s"%outfile_names[1]
        if len(outfile_names) > 2:
            print "grav_anomaly.py: error: Too many filenames specified [unresolved, check code]"
            sys.exit(1)

    for i in range(len(data_pairs)):

        # Get the current x and y data, and the current output file
        outfile = outfiles[i]
        x_data, y_data = data_pairs[i]

        # Print header lines, if desired
        if GMT_print_headers:
            if i == 1:
                outfile.write("#%-24s %-8s\n"%("geoid file", geoidfile_name))
            else:
                outfile.write("#%-33s\n"%("raw acceleration (no geoid file)"))
            outfile.write("#%-24s %-8s m\n"%("elevation:", elevation))
            outfile.write("#%-24s %-8s m\n"%("x range:", x_range))
            outfile.write("#%-24s %-8s\n"%("number of x steps:", num_steps))
            outfile.write("#%-24s %-8s\n"%("", ""))
            outfile.write("#%16s %16s\n"%("X Value", "Gravity Anomaly"))
            outfile.write("#%16s %16s\n"%("-"*16, "-"*16))

        # Write out the data in two columns - X values, Y values
        for i in range(len(x_data)):
            outfile.write("%16.4f %16.8g\n"%(x_data[i],y_data[i]))

    # More stuff to print if we're not in quiet mode (again doesn't require
    # that we be in verbose mode)
    if not quiet_mode:
        printstring = "Data written to %s"%outfile_names[0]
        if len(outfile_names) == 1:
            print printstring
        if len(outfile_names) == 2:
            print printstring + " and %s"%outfile_names[1]
        if len(outfile_names) > 2:
            print "grav_anomaly.py: error: Too many filenames specified [unresolved, check code]"
            sys.exit(1)


def plotter(data_sets):
    """
    Plots the given data on a single plot, using numpy.plt
    """

    # Import plt from numpy here, so you don't waste time if you're using GMT
    # mode instead of the built-in plotter
    import matplotlib.pyplot as plt

    # Make a list of colors for the lines that we can iterate over
    colors =  ["Red","Green","Blue","Orange","Cyan",
               "Purple","Black","Magenta","Navy","Teal"]

    # Iterate through whatever data we have in data_sets and plot it up
    for i in range(len(data_sets)):
        name = data_sets.keys()[i]
        x, y = data_sets[name]
        try:
            plt.plot(x,y,label=name,color=colors[i],linewidth=2.0)
        except IndexError:
            plt.plot(x,y,label=name,color='k',linewidth=2.0)

    # Window dressing for the plots
    plt.xlabel('Distance from center (km)')
    plt.ylabel('Gravity anomly (mGal)')
    #plt.legend()
    plt.legend(loc='lower right')
    #plt.legend(loc='upper right')
    #plt.legend(loc='best')
    #plt.legend(bbox_to_anchor=(1.1, 1.05)) #Overlapping the upper right
    plt.grid(True)

    # Show the plot, then close plt once the window's been closed
    if not quiet_mode:
        print "Data plotted in external window"
    plt.show()
    plt.close()


def current_state_printer():
    """
    Print information about how things are running
    """
    # Make pretty text for the current calculation mode
    if bouguer_mode:
        current_mode = "BOUGUER"
    elif anomaly_mode:
        current_mode = "FREE-AIR ANOMALY"
    else:
        current_mode = "RAW ACCELERATION"

    if curved_mode:
        curved_or_flat = "CURVED"
    else:
        curved_or_flat = "FLAT"

    # Print the other variables
    print "grav_anomaly.py, version %s:"%__version__
    print "        Running in %s mode"%current_mode
    print "        with a %s geometry"%curved_or_flat
    print "        (geoid file name: %s)"%geoidfile_name
    print "        Outputting in %s mode"%("visual","GMT")[GMT_mode]
    print ""
    print "        range of x values covered: 0 - %g m"%x_range
    print "        spacecraft elevation: %g m"%elevation
    print "        number of steps: %g"%num_steps
    print "        verbosity mode: %s %s %s"%(
                        ("","Verbose")[verbose_mode],
                        ("","quiet")[quiet_mode],
                        ("Quiet","")[verbose_mode or quiet_mode])


######## Command-Line Implementation ###########################################

if __name__ == '__main__':

    # Initiate the command-line argument parser
    parser = argparse.ArgumentParser(
        description = ("Calculate the gravitational anomaly field over models "
                      "of impact basins"),
        epilog = "Please contact Dave Blair (dblair@purdue.edu) with bugs or questions.")

    # The one positional argument is the name(s) of the gravfile(s) we're
    # looking at
    parser.add_argument("gravfile", nargs="*")

    # Define our various options and switches
    parser.add_argument("-f","--freeair",
        help = "turn on 'free-air' (baseline file comparison) mode",
        action = "store_true")
    parser.add_argument("--geoidname",
        help = "set GEOIDNAME as the file to use for free-air or Bouguer "+\
               "calculation (set with a separate flag)")
        #help = "turn on 'free-air' mode, using GEOIDNAME as the baseline "+\
        #       "for comparison (the 'geoid')")
    parser.add_argument("-b","--bouguer",
        help = "turn on 'Bouguer' mode (subtracts out topography)",
        action = "store_true")
    parser.add_argument("-a","--acceleration",
        help = "skip the 'free-air' calculation and just calculate "+\
               "acceleration over the model. Useful for combining with GMT "+\
               "output and creating .xy geoid files.",
        action = "store_true")
    parser.add_argument("-s","--subdivide", metavar = "S",
        help = "take each element that's passed into the program and "+\
               "subdivide it into 2^S smaller elements")
    parser.add_argument("-S",
        help = "subdivide elements as above, but perform a default number "+\
               "of subdivisions (specified in the OPTIONS section of this "+\
               "program's source code)",
        action = "store_true")
    parser.add_argument("--nosubdivision",
        help = "do not subdivide elements (overrides OPTIONS section default)",
        action = "store_true")

    parser.add_argument("-g","--GMT",
        help = "generate text instead of visual output: 'GMT mode'",
        action = "store_true")
    parser.add_argument("--GMTheaders",
        help = "print commented-out headers at the top of text output files",
        action = "store_true")
    parser.add_argument("-p","--plot",
        help = "plot results with the built-in plotter instead of writing to files",
        action = "store_true")

    parser.add_argument("-r","--range", metavar = "R",
        help = "set spacecraft horizontal range to R in meters [flat models only]")
    parser.add_argument("-n","--numsteps", metavar = "N",
        help = "set the number of horizontal summation steps to N")
    parser.add_argument("-e","--elevation", metavar = "Y",
        help = "set the spacecraft elevation to Y (in meters)")

    parser.add_argument("-c","--curved",
        help = "use a curved spacecraft path",
        action = "store_true")
    parser.add_argument("-t","--thetarange", metavar = "T",
        help = "set spacecraft theta range to T in degrees [curved models only]")
    parser.add_argument("--planetradius", metavar = "Rp",
        help = "define the radius of the planet (for use with a curved model)")

    parser.add_argument("-v","--verbose",
        help = "provide extra feedback while running",
        action = "store_true")
    parser.add_argument("-q","--quiet",
        help = "absolutely no text output (except errors). This is quieter "+\
               "than simply running without being in verbose mode.",
        action = "store_true")

    parser.add_argument("--examples",
        help = "show some example usage",
        action = "store_true")
    parser.add_argument("--defaults",
        help = "print the default parameters from this program's OPTIONS section",
        action = "store_true")
    parser.add_argument("--version",
        help = "print the version number",
        action = "store_true")
    args = parser.parse_args()

    # Define the text to be used for "--examples"
    examplestext = (
    "\n"
    "Example usage for grav_anomaly.py:                            \n"
    "                                                              \n"
    "grav_anomaly.py foo.grav                                      \n"
    "        calculates anomaly of foo.grav vs. the baseline file  \n"
    "        specified in the code, and shows the results using    \n"
    "        the built-in Matplotlib plotter                       \n"
    "                                                              \n"
    "grav_anomaly.py -g foo.grav                                   \n"
    "        as above, but writes initial/final/delta results      \n"
    "        to files foo_acc.xy and foo_anom.xy                   \n"
    "                                                              \n"
    "grav_anomaly.py -q -e 30000 -g foo.grav                       \n"
    "        same as above, but sets elevation to 30,000 m         \n"
    "        and supresses most text feedback while running        \n"
    "                                                              \n"
    "grav_anomaly.py -a foo.grav bar.xy baz.grav                   \n"
    "        same as above but for three input files, one          \n"
    "        of which is an XY file)                               \n"
    "                                                              \n"
    "grav_anomaly.py -a -g foo.grav                                \n"
    "        same as above, but calculating raw acceleration only, \n"
    "        so there's only the foo_acc.out file to write         \n")

    # Define some specific error text
    missingitem_errortext = (
        "\n"
        "grav_anomaly.py: error: please specify at lease one input file\n"
        "if not running with the --version, --defaults, or --examples flag.")
    badfile_errortext = (
        "\n"
        "grav_anomaly.py: error: file not found. Please check that files\n"
        "exist and try again.")
    conflict_errortext = (
        "\n"
        "grav_anomaly.py: error: conflicting arguments given. Please\n"
        "check options and try again.")
    gmt_toomanyfiles_errortext = (
        "\n"
        "grav_anomaly.py: error: GMT mode can only accept a single"
        "input file. Please specify a single .grav file and try again.")
    bouguer_on_nongravfile_errortext = (
        "\n"
        "grav_anomaly.py: error: Bouguer mode is only supported for .grav\n"
        "input files (geoid file may be in xy-file format). Please specify\n"
        "a different input file.")

    # Now to interpret the arguments. We'll start with the ones that'll
    # print something and then just terminate the program
    if args.examples:
        print examplestext
        sys.exit(1)
    if args.defaults:
        current_state_printer()
        sys.exit(1)
    if args.version:
        print "grav_anomaly.py, version %s"%__version__
        sys.exit(1)

    # Then we move on to the real math/processing options statements (while
    # keeping track of conflicting arguments)
    conflicting_args_calculationmode = 0
    if args.geoidname:
        #anomaly_mode = True #For convenience
        # No comment on whether we're in Bouguer mode or not
        geoidfile_name = args.geoidname
        #conflicting_args_anomalymode += 1
    if args.acceleration:
        anomaly_mode = False
        bouguer_mode = False
        conflicting_args_calculationmode += 1
    if args.freeair:
        anomaly_mode = True
        bouguer_mode = False
        #conflicting_args_anomalymode += 1
        conflicting_args_calculationmode += 1
    if args.bouguer:
        anomaly_mode = True
        bouguer_mode = True
        #conflicting_args_anomalymode += 1
        conflicting_args_calculationmode += 1
    if args.subdivide:
        number_of_subdivisions = eval(args.subdivide)
    if args.S:
        number_of_subdivisions = default_number_of_subdivisions
    if args.nosubdivision:
        number_of_subdivisions = 0
    if args.curved:
        curved_mode = True

    # Then check for GMT or plotting mode
    if args.GMT:
        GMT_mode = True
    if args.GMTheaders:
        GMT_print_headers = True
    if args.plot:
        GMT_mode = False

    # Then look at default parameter overrides
    if args.range:
        x_range = float(args.range)
    if args.numsteps:
        num_steps = float(args.numsteps)
    if args.elevation:
        elevation = float(args.elevation)
    if args.planetradius:
        planet_radius = float(args.planetradius)
    if args.thetarange:
        theta_range = float(args.thetarange)

    # Then at verbosity arguments
    conflicting_args_verbosity = 0
    if args.verbose:
        verbose_mode = True
        quiet_mode = False
        conflicting_args_verbosity += 1
    if args.quiet:
        verbose_mode = False
        quiet_mode = True
        conflicting_args_verbosity += 1

    # Make sure we don't have any conflicting options, and that the user has
    # specified at least one file to plot (the minimal "grav_anomaly.py
    # foo.grav" input)
    if ((conflicting_args_verbosity > 1) or
        #(conflicting_args_anomalymode > 1) or
        (conflicting_args_calculationmode > 1)):
        print conflict_errortext
        parser.print_usage()
        sys.exit(1)
    if not args.gravfile:
        print missingitem_errortext
        parser.print_usage()
        sys.exit(1)

    # The final bit of command-line processing is grabbing the filenames of the
    # data files we're actually interested in! We'll first create the set that's
    # going to hold all of our data. If we're in baseline mode, we'll make the
    # baseline file the first member of this (we already have a flag telling us
    # how to treat this first file)

    # Check to make sure that, if we're in GMT mode, there's only one grav file
    # specified
    if (len(args.gravfile) > 1) and (GMT_mode):
        print gmt_toomanyfiles_errortext
        sys.exit(1)

    # Now that we're done with checks, on to a little feedback: if we *are* in
    # verbose mode, give the same output as "defaults" (except for the current
    # run, accounting for all command-line changes
    if verbose_mode:
        current_state_printer()

    # If we're in anomaly_mode, use the baseline/geoid file (whether defined in
    # OPTIONS section or by a command-line argument) as the first entry in
    # datafiles. Afterwards, process the rest of the given files and add them to
    # the list.
    datafiles = []
    if anomaly_mode:
        datafiles.append(open(geoidfile_name, 'r'))
    for this_filename in args.gravfile:
        datafiles.append(open(this_filename,'r'))

    # On to the meat & potatoes of data processing. We'll use a filetype- and
    # subdivision-aware function (get_data_from_file()) to process our input
    # files and return data sets. We'll also create a dictionary called
    # "data_sets" up here for storing information, such that data_sets["Label"]
    # = (x_values, {anomaly/acceleration values}). "Label" is usually filename,
    # except in GMT mode, where it's Final/Delta.
    data_sets = {}
    total_iterations = len(datafiles)
    this_iteration = 0

    # If we're in anomaly mode, we'll first process the baseline file, and then
    # delete it from datafiles so that we can run the next process regardless
    # of whether or not we're in anomaly mode
    if anomaly_mode:

        # The baseline file is the first file listed in datafiles, so we can
        # just use that to find it here, as long as we know we're in anomaly
        # mode
        geoidfile = datafiles[0]

        # Process the baseline file using a function defined in the main body of
        # this program that's filetype- and subdivision-aware
        geoid_data = get_data_from_file(geoidfile)[0]

        data_sets["Baseline"] = geoid_data

        # Remove the baseline file from the list, because we're done with it
        del datafiles[0]

    # Now do pretty much the same thing, but for the rest of the files. We'll go
    # through all of datafiles; if we were in anomaly mode, we've already
    # deleted the baseline file
    for this_file in datafiles:

        # Assign a name to the data from the filename
        this_dataname = this_file.name

        # Process the current file using a function defined in the main body of
        # this program that's filetype- and subdivision-aware
        this_data, this_data_basename = get_data_from_file(this_file)

        # Then, if we're in anomaly mode, calculate the difference between the
        # current data and the baseline
        if anomaly_mode:
            this_anomaly_data = delta_calculator(geoid_data, this_data)

            # If we're also in GMT mode, then there's special output nomenclature
            # for parsing with the GMT_outputter.
            if GMT_mode:
                data_sets["Final"] = this_data
                data_sets["Delta"] = this_anomaly_data

            # If we're not in GMT mode, we can discard the "Final" data and
            # instead preserve the filename, since the graphical plotter can do
            # more than one file at a time
            else:
                data_sets[this_dataname] = this_anomaly_data

        # If we're not in anomaly mode, it's simpler, because we don't have
        # to deal with the delta calculator. We do still need to account for
        # GMT output mode, though, in order to name the data set properly
        else:
            if GMT_mode:
                data_sets["Final"] = this_data
            else:
                data_sets[this_dataname] = this_data


    # We're on to the final step: outputting the data in a usable format. We'll
    # check to see if we're in GMT or Matplotlib (graphical) plotting mode, and
    # then implement that output mode

    # Check first to see if we're in GMT mode (table output), and if so, use
    # GMT_outputter()
    if GMT_mode:

        # Grab the basename from the data file, which can be any of the strings
        # specified in the gravfile_suffixes variable in OPTIONS. If we're in
        # Bouguer mode, alert the user now that they can only use a raw .grav
        # file, and quit if that's not the case
        output_basename = args.gravfile[0]
        for file_suffix in gravfile_suffixes:
            output_basename = re.sub(file_suffix,"",output_basename)

        # Are we in anomaly mode (write out initial, final, and delta files) or
        # are we in acceleration mode (just write out the final state file)? In
        # either case, make everything sets so that we can iterate over an
        # unknown number of them in the GMT_outputter function
        if anomaly_mode:

            #Name and open the acceleration (final) and delta-file based on
            #whether we're in free-air or Bouguer mode, and then open it
            if bouguer_mode:
                final_outfile_name = output_basename + GMT_suffixes[1][0]
                delta_outfile_name = output_basename + GMT_suffixes[1][1]
            else:
                final_outfile_name = output_basename + GMT_suffixes[0][0]
                delta_outfile_name = output_basename + GMT_suffixes[0][1]
            final_outfile = open(final_outfile_name,'w')
            delta_outfile = open(delta_outfile_name,'w')

            outfile_names = [final_outfile_name, delta_outfile_name]
            outfiles = [final_outfile, delta_outfile]

            # Make sure that data_sets is labeled with Final and Delta
            data_sets_ordered = {}
            data_sets_ordered["Final"] = data_sets["Final"]
            data_sets_ordered["Delta"] = data_sets["Delta"]
            del data_sets

        else:
            final_outfile_name = output_basename + GMT_suffixes[2]
            final_outfile = open(final_outfile_name, 'w')
            outfile_names = [final_outfile_name]
            outfiles = [final_outfile]

            # Make sure that data_sets only contains the final data
            data_sets_ordered = data_sets
            del data_sets

        # Output as files ready for GMT processing
        GMT_outputter(data_sets_ordered, outfiles, outfile_names)

    # If we're not in GMT mode, use the built-in plotter. The part of the data
    # that gets plotted (anomaly or acceleration) depends on which mode we're
    # in.
    else:
        if anomaly_mode:
            del data_sets["Baseline"]
            plotter(data_sets)
        else:
            plotter(data_sets)
