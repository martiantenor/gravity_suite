#!/usr/bin/env python
# A program that generates a file with information about the "elements" in an
# iSALE model, using both the iSALE input (.inp) and data (X,Y,density) files.
#
# Contact Dave Blair (dblair@purdue.edu) with questions
#
# (c) David Blair, 2015. This work is licensed under a Creative Commons
# Attribution-NonCommercial-ShareAlike Unported License
# (http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US)

from __future__ import division
import sys, argparse, re
from math import pi

__version__ = "3"
# For ToDo and Changelog, please see file "nodes_from_iSALE.py.changelog"


######## Options ###############################################################

# What's the suffix of the output file going to be?
outfile_suffix = ".grav"

# What's the suffix of the input file?
datafile_suffix = ".txt"

# Discard cells that have y > 0? (Generally, leave this set to False and set it
# on a case-by-case basis with "-0 / --discardabovezero" at runtime)
discard_above_zero = False
discard_above_zero_outfile_suffix = "_nothingAbove0km.grav"

# Plot up node results (as a check) visualy?
#plotting_mode = True
plotting_mode = False

# Stop after plotting the nodes?
#plotonly_mode = True
plotonly_mode = False

# Lots of text output to make sure things are going right?
#verbose_mode = True
verbose_mode = True

# Default values to get around some issues:
max_cell_size = 9e99
max_layer_height = 9e99



######## Main Program ##########################################################

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
        "Makes a simple picture of the element, using pylab"

        x_coords = []
        y_coords = []

        # Iterate over nodes so that the plotting is independent of the number
        # of nodes being fed into this function
        for i in range(len(self.nodes)):
            x, y = self.nodes[i][0], self.nodes[i][1]
            x_coords.append(x), y_coords.append(y)

        pylab.fill(x_coords, y_coords, color)
        #pylab.show()


def find_triangle_area(A,B,C):
    """
    Find the area of a triangle in 2D from the coordinates of its nodes.
    """
    return abs(.5* (A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1])) )


def isale_inp_parser(inpfile):
    """
    Takes an iSALE input file and generates a replica of the iSALE mesh.  That
    data is returned as two lists, one of x values and one of y values.
    """

    # Read the file, then close it, storing all the lines in a set instead
    inpfile_lines = inpfile.readlines()
    inpfile.close()

    # Go through all the lines and look for certain necessary keywords
    for this_line in inpfile_lines:

        # If the line starts with a dash or hash, ignore it
        if this_line.startswith("-") or this_line.startswith("#"):
            continue

        # Otherwise, grab the first word, call that the "keyword", and then
        # compare it via a series of if/elif/else statements
        this_keyword = this_line.split()[0].upper()

        # Horizontal spacing in the high-res zone
        if this_keyword == "GRIDH":

            GRIDH = this_line.strip().split(":")[1:]
            cells_left_of_hr = int(GRIDH[0])
            cell_width_of_hr = int(GRIDH[1])
            cells_right_of_hr = int(GRIDH[2])
            total_cell_width = sum([cells_left_of_hr,
                                    cell_width_of_hr,
                                    cells_right_of_hr])

            if verbose_mode:
                print "%42s: %i"%("Cells to left of high-res zone",
                                  cells_left_of_hr)
                print "%42s: %i"%("Cell width of high-res zone",
                                  cell_width_of_hr)
                print "%42s: %i"%("Cells to right of high-res zone",
                                  cells_right_of_hr)
                print "%42s: %i"%("Total cell width of model",
                                  total_cell_width)

        # Vertical spacing in the high-res zone
        elif this_keyword == "GRIDV":
            GRIDV = this_line.strip().split(":")[1:]
            cells_below_hr = int(GRIDV[0])
            cell_height_of_hr = int(GRIDV[1])
            cells_above_hr = int(GRIDV[2])
            total_cell_height = sum([cells_below_hr,
                                     cell_height_of_hr,
                                     cells_above_hr])

            if verbose_mode:
                print "%42s: %i"%("Cells below high-res zone",
                                  cells_below_hr)
                print "%42s: %i"%("Cell height of high-res zone",
                                  cell_height_of_hr)
                print "%42s: %i"%("Cells above high-res zone",
                                  cells_above_hr)
                print "%42s: %i"%("Total cell height of model",
                                  total_cell_height)

        # Grid extension factor for the low-res zone
        elif this_keyword == "GRIDEXT":
            GRIDEXT = this_line.split()[-1].upper()
            extension_factor = float(GRIDEXT.replace("D","e"))

            if verbose_mode:
                print "%42s: %g"%("Extention factor",
                                  extension_factor)

        # Mimimum (high-res zone) cell size
        elif this_keyword == "GRIDSPC":
            GRIDSPC = this_line.split()[-1].upper()
            hr_cell_size = float(GRIDSPC.replace("D","e"))

            if verbose_mode:
                print "%42s: %g m"%("Cell size in high-res zone",
                                     hr_cell_size)

        # Maximum cell size
        elif this_keyword == "GRIDSPCM":
            GRIDSPCM = this_line.split()[-1].upper()
            max_cell_size = float(GRIDSPCM.replace("D","e"))

            if verbose_mode:
                print "%42s: \"%i\""%("Maximum cell size",
                                    max_cell_size)

        # The position (number of elements) up from the bottom of the mesh that
        # defines the (0,0) point
        elif this_keyword == "LAYPOS":

            LAYPOS = this_line.split(":")

            # Since the formatting for this one is kind of loose, we'll do a
            # try/except loop within a for loop that iterates over the values in
            # LAYPOS - if it successfully converts to an integer, we'll add it
            # to a list, and then we'll pull out the largest value
            all_layer_heights = []
            for value in LAYPOS:
                try:
                    all_layer_heights.append(int(value))
                except ValueError:
                    pass
            max_layer_height = max(all_layer_heights)

            if verbose_mode:
                print "%42s: %i"%("Cell height of tallest layer",
                                  max_layer_height)

        # None of the above - ignore it
        else:
            continue

    # Now that we have all our variables, it's time to do some math! First,
    # we'll figure out our actual maximum cell size by interpreting the weird
    # "less than zero is a multiplier on the minimum cell size" syntax
    #if max_cell_size < 0:
    #    max_cell_size = -max_cell_size * hr_cell_size
    #if verbose_mode:
    #    print "%42s: %g m"%("Maximum cell size",max_cell_size)

    # Moving on, we'll figure out how many cells are above and below the 0,0
    # mark, and then how many of these are hr cells versus extended cells
    cells_above_00 = total_cell_height - max_layer_height
    cells_below_00 = total_cell_height - cells_above_00
    hr_cells_above_00 = cells_above_00 - cells_above_hr
    hr_cells_below_00 = cells_below_00 - cells_below_hr
    extended_cells_above_00 = cells_above_00 - hr_cells_above_00
    extended_cells_below_00 = cells_below_00 - hr_cells_below_00
    if verbose_mode:
        print "%42s: %i"%("High-res cells below (0,0)",
                          hr_cells_below_00)
        print "%42s: %i"%("Extended cells below (0,0)",
                          extended_cells_below_00)
        print "%42s: %i"%("High-res cells above (0,0)",
                          hr_cells_above_00)
        print "%42s: %i"%("Extended cells above (0,0)",
                          extended_cells_above_00)

    # Now we'll do the same thing for the horizontals, except that we know that
    # 0 in the x direction always lies on the axis. Note that this is an
    # assumption, but it should always be true in the models generated within
    # out group; we'll include a checker just to be sure
    if cells_left_of_hr != 0:
        print ("ERROR: High-res zone doesn't start at x=0! This is \n"
               "currently  not allowed in this program. Contact Dave Blair at\n"
               "dblair@purdue.edu to  request futher development of the code.")

    hr_cells_right_of_00 = cell_width_of_hr
    extended_cells_right_of_00 = total_cell_width - cell_width_of_hr
    if verbose_mode:
        print "%42s: %i"%("High-res cells right of (0,0)",
                          hr_cells_right_of_00)
        print "%42s: %i"%("Extended cells right of (0,0)",
                          extended_cells_right_of_00)

    # With this done, we can move on to actually creating the mesh. We'll start
    # by making empty sets for the x and y values
    x_values = []
    y_values = []

    # First, we'll figure out the x values, since they're a little easier
    this_x = 0
    cell_width = hr_cell_size

    # When we set up the iteration, add one to the total cell width so that we
    # wind up with nodes on both the far left and the far right (i.e. one more
    # node than we have elements)
    for cell_number in range(total_cell_width + 1):

        # If we're in the high-res zone, then it's just multiplication
        if cell_number <= hr_cells_right_of_00:
            this_x = cell_number * cell_width

        # If we're past the high-res zone, then each time we'll update the
        # cell width by multiplying by the extension factor, unless this number
        # would surpass the maximum cell width; in that case, it's just the
        # maximum cell width
        else:
            cell_width = cell_width * extension_factor

            # If this would put the cell width at a value more than the maximum,
            # set it to the maximum instead
            if cell_width > max_cell_size:
                cell_width = max_cell_size

            this_x = this_x + cell_width

        x_values.append(this_x)

    # Now we'll move on to the y nodes above 0,0 (positive values), using very
    # similar code to what we just did
    this_y = 0
    cell_height = hr_cell_size
    for cell_number in range(cells_above_00 + 1):
        if cell_number <= hr_cells_above_00:
            this_y = cell_number * cell_height
        else:
            cell_height = cell_height * extension_factor
            if cell_height > max_cell_size:
                cell_height = max_cell_size
            this_y = this_y + cell_height
        y_values.append(this_y)

    # Finally, we'll do the y values below 0,0, omitting 0.0 since it's already
    # in there
    this_y = 0
    cell_height = hr_cell_size
    for cell_number in range(cells_below_00 + 1):
        if cell_number == 0:
            continue
        elif cell_number <= hr_cells_below_00:
            this_y = -cell_number * cell_height
        else:
            cell_height = cell_height * extension_factor
            if cell_height > max_cell_size:
                cell_height = max_cell_size
            this_y = this_y - cell_height
        y_values.append(this_y)

    return x_values, y_values


def isale_data_parser(datafile, x_meshvalues, y_meshvalues):
    """
    Read data from an iSALE (hydrocode model) data file and interpret it as
    a set of elements, using a set of all the possible x and y values of
    "nodes" (cell corners) in the iSALE mesh
    """

    # Read the contents of the data file and then close it immediately
    datalines = datafile.readlines()
    datafile.close()

    # Go through all the lines; for each (x,y) point in the file, we'll figure
    # out which x and y coordinate in the mesh it's closest to, use those as the
    # corner coordinate points, and create an element that we'll add to
    # celldata
    celldata = []
    for this_line in datalines:

        # Ignore comment lines
        if this_line.startswith("#"):
            continue

        # Interpret the line
        line_data = this_line.split()
        data_x = float(line_data[0])
        data_y = float(line_data[1])
        this_density = float(line_data[2])
        #this_temperature = float(line_data[3]) #Not currently in use

        # Before doing anything else, check if we're supposed to discard y
        # values above zero, and check if our y value is above zero
        if discard_above_zero and (data_y > 0):
            continue

        # Compare the x and y coordinates against our list of possible values,
        # and figure out the x values bounding data_x
        for i in range(len(x_meshvalues)):

            current_xposition_in_mesh = x_meshvalues[i]

            # See if the value we're looking at is the first one greater than
            # our "data_x" value (the centroid/data coordinate)
            if data_x > current_xposition_in_mesh:
                pass
            else:
                element_rightedge = current_xposition_in_mesh
                element_leftedge = x_meshvalues[i-1]
                break

        # Do the same thing for the y values
        for i in range(len(y_meshvalues)):

            current_yposition_in_mesh = y_meshvalues[i]

            if data_y > current_yposition_in_mesh:
                pass
            else:
                element_topedge = current_yposition_in_mesh
                element_bottomedge = y_meshvalues[i-1]
                break

        # Create an element with these parameters
        this_centroid = ( (element_leftedge + element_rightedge)/2.,
                     (element_bottomedge + element_topedge)/2. )
        this_nodedata = [ (element_rightedge, element_topedge),
                      (element_rightedge, element_bottomedge),
                      (element_leftedge, element_bottomedge),
                      (element_leftedge, element_topedge) ]
        this_element = Element( this_centroid, this_density, this_nodedata)

        # Add this element to the master list
        celldata.append(this_element)

    return celldata


def grav_outputter(elements,outfile):
    """
    Writes out data from a given element set to a filename taken from the odb's
    filename and a file suffix & extension. This data is element ID, centroid X,
    Y, density, volume, and the nodal coordinates
    """
    # Headers
    outfile.write("%s data\n"%(outfile.name))
    outfile.write("%10s %16s %16s %16s %16s | %32s\n" \
                      %("element ID","centroid X","centroid Y",
                        "density","volume","node information"))
    # Dashes
    outfile.write("%10s %16s %16s %16s %16s-|-%32s\n"%\
                      ("-"*10,"-"*16,"-"*16,"-"*16,"-"*16,"-"*32))
    # Content
    fake_elemID = 1
    for this_element in elements:
        outfile.write("%10i %16.4f %16.4f %16.4f %16.8g | %16s\n"\
                           %(fake_elemID,
                             this_element.centroid[0],
                             this_element.centroid[1],
                             this_element.density,
                             this_element.volume,
                             this_element.nodes))
        fake_elemID += 1

    # Let the user know it's been done
    if verbose_mode:
        print "Element data for gravity calculation written to %s"%outfile.name


def mesh_plotter(x_meshvalues,y_meshvalues):
    """
    Plots up the empty mesh "nodes" (cell corners) to show the structure of the
    iSALE input mesh, as a check
    """

    # Create "xy_pairs", which is a set containing all possible combinations
    # of x and y values from the mesh
    xy_pairs = []
    for this_x in x_meshvalues:
        for this_y in y_meshvalues:
            xy_pairs.append((this_x,this_y))

    # Get a list of all the x and y coordinates from the coordinate
    # pairs. Because x_meshvalues and y_meshvalues aren't necessarily the
    # same length (the plot might not be square), the resulting x_components
    # and y_components sets are *not* the same thing as the x_values and
    # y_values sets.
    x_components = []
    y_components = []
    for pair in xy_pairs:
        x_components.append(pair[0])
        y_components.append(pair[1])

    # Look for the dimensions of each so that we can use this to make a
    # square plot in a minute
    xmin = min(x_components)
    xmax = max(x_components)
    x_range = xmax - xmin
    ymin = min(y_components)
    ymax = max(y_components)
    y_range = ymax - ymin
    max_range = max(x_range, y_range)

    # Now we'll plot it all up using pylab/matplotlib
    import pylab

    # Set the figure to 6x6", and specify axis limits so it's square
    pylab.figure(figsize=(6,6))
    pylab.xlim(0,max_range)
    pylab.ylim(ymin,ymax)

    # Label the axes and the plots
    pylab.xlabel("X coordinate (m)")
    pylab.ylabel("Y coordinate (m)")
    #pylab.title("%s"%inpfilename.strip(".inp"))
    pylab.title(re.sub(datafile_suffix,"",datafile_name))

    # Create a scatter plot, using a 0-linewidth black point marker
    pylab.scatter(x_components, y_components,
                  marker=".",
                  lw=0,
                  color="black")
    pylab.show()


######## Command-Line Implementation ###########################################

if __name__ == "__main__":

    # Initiate the command-line parser
    parser = argparse.ArgumentParser(
        description=("Generate a .grav file from an iSALE input file and a data \n"
                     "file with columns X, Y, density"),
        epilog="Please contact Dave Blair (dblair@purdue.edu) with bugs or questions.")

    # Define our positional arguments (the filenames)
    parser.add_argument("inpfile_name")
    parser.add_argument("datafile_name")
    
    # Define the various flags & runtime options
    parser.add_argument("-0", "--discardabovezero",
        help = "ignore cells with y values above zero. This is useful for "+\
               "using the time = 0 data to generate a 'geoid' file, because "+\
               "it will ignore the impactor",
        action = "store_true")
    parser.add_argument("-p","--plotting","--showmesh",
        help = "plot up the resulting iSALE mesh instead of continuing on to " +\
             "write a .grav file",
        action = "store_true")
    parser.add_argument("-P","--plotonly",
        help = "as above, but stop after plotting the mesh",
        action = "store_true")
    parser.add_argument("-v","--verbose",
        help = "print extra stuff while running",
        action = "store_true")
    parser.add_argument("--version",
        help = "print the version number",
        action = "store_true")
    args = parser.parse_args()

    # Check for any options that should run and then terminate the program
    if args.version:
        print "grav_from_iSALE.py, version %s"%__version__
        sys.exit()

    # Run the parser and assign real variables to the positional arguments
    # listed above
    inpfile_name = args.inpfile_name
    datafile_name = args.datafile_name
    if args.discardabovezero:
        discard_above_zero = True
        outfile_suffix = discard_above_zero_outfile_suffix
    if args.verbose:
        verbose_mode = True
        print "Processing iSALE input file %s ..."%inpfile_name
        print "%42s: %s"%("Ignoring cells above y=0",
                          ("No","Yes",discard_above_zero))
    if args.plotting:
        plotting_mode = True
    if args.plotonly:
        plotting_mode = True
        plotonly_mode = True

    # First command line argument is the input filename, second argument is the
    # data filename
    inpfile = open(inpfile_name, 'r')
    datafile = open(datafile_name, 'r')

    # Run the input file parser defined above on the input file in order to
    # retrieve all the x and y values in the iSALE mesh
    x_meshvalues, y_meshvalues = isale_inp_parser(inpfile)
    x_meshvalues.sort()
    y_meshvalues.sort()

    # Now we're on to output behavior. First, check if we're supposed to plot up
    # the nodal positions using the NumPy plotter
    if plotting_mode:
        mesh_plotter(x_meshvalues,y_meshvalues)

    # After plotting up the nodes as a picture (or if we skipped that step),we
    # can continue on to processing the iSALE data
    if not plotonly_mode:
        celldata = isale_data_parser(datafile, x_meshvalues, y_meshvalues)
        if verbose_mode:
            print "Returned element set of length %i from file %s"%(len(celldata),
                                                                    inpfile_name)

        # Finally, we can write out the output file
        outfile = open(re.sub(datafile_suffix, outfile_suffix, datafile_name),'w')
        grav_outputter(celldata, outfile)
        outfile.truncate()
        outfile.close()
