#!/usr/bin/env python
# A program to grab output from a given Abaqus .odb file (that you have open in
# CAE) for processing in another program (grav_anomaly.py) that gives plots of
# the gravity anomaly over that model.
#
# Contact Dave Blair (dblair@purdue.edu) with questions
#
# (c) David Blair, 2015. This work is licensed under a Creative Commons
# Attribution-NonCommercial-ShareAlike Unported License
# (http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US)

__version__ = "3"
# For changelog, see grav_plugin.py.changelog

# Import the necessary Abaqus/CAE libraries in a try:except loop, so that if the
# libraries aren't available (e.g. if trying to run from a BASH shell) the code
# can raise a useful error and quit.
import sys
try:
    from abaqus import *
    from odbAccess import openOdb
    from odbMaterial import *
    from odbSection import *
except:
    print ("grav_from_Abaqus.py: ERROR: Abaqus/CAE libraries aren't available.\n"
           "Please be sure this program is launched from within the Abaqus/CAE\n"
           "environment.")
    print ("grav_from_Abaqus.py, version %s"%__version__)
    sys.exit()


######## Options ###############################################################

# Is this a 2D or a 3D model? If it's 3D and axisymmetric, we'll only look at
# the nodes on the Z=0 plane and treat them as if they're annuli.
threedee_axi_mode = False

# Output first frame, last frame, or both?
output_firstframe = True
output_lastframe = True

# Pretend some other frame is the last frame? If not, enter "-1"
last_frame_number = -1
#last_frame_number = 78

# Specify the step you want to recover from by name. If there's only one step,
# or you want to default to the last step, just use a blank string ("")
#default_step_name = "visco"
default_step_name = ""

# Output density files (for first and last frame, if enabled above)?
output_density = True

# Output pressure files (for first and last frame, if enabled above)?
output_pressure = True

# Plot values on the deformed mesh, or the undeformed mesh?
output_undeformed = False

# Filename options
output_suffix_firstframe = "_ff.grav"
output_suffix_firstframedensity = "_ffdensity.xy"
output_suffix_firstframepressure = "_ffpressure.xy"
#output_suffix_lastframe = ".grav"
if last_frame_number != -1:
    output_suffix_lastframe = "_frame%s.grav"%last_frame_number
else:
    output_suffix_lastframe = ".grav"
output_suffix_lastframedensity = "_density.xy"
output_suffix_lastframepressure = "_pressure.xy"
output_suffix_nodeflastframedensity = "_nodefdensity.xy"
output_suffix_nodeflastframepressure = "_nodefpressure.xy"

# Print lots of things while running?
verbose_mode = True



######## Main Program ##########################################################

class Element:
    """
    A simple "Element" class, for storing data and calculating mass
    """

    def __init__(self, elemID, centroid, nodes, volume, density, pressure=0.0):
        self.elemID = elemID
        self.centroid_x = centroid[0]   # x coord of the centroid
        self.centroid_y = centroid[1]   # y coord of the centroid

        try:
            self.centroid_z = centroid[2]   # z coord of centroid (for 3D models)
        except IndexError:
            pass

        self.nodes = nodes
        self.volume = volume
        self.density = density
        self.pressure = pressure


def get_odb_from_cae():
    """
    Grab the current viewport, open the .odb file being displayed, and retrieve
    the name of that .odb file and the path in which it sits.
    """
    current_viewport = session.viewports[session.currentViewportName]
    odb_file = current_viewport.displayedObject
    odb_pathname = os.path.split(odb_file.path)[0]
    odb_filename = os.path.split(odb_file.path)[1].split(".")[0]  #.odb is removed

    return odb_file, odb_pathname, odb_filename


def get_odb_from_commandline():
    """
    Grab a specified .odb file from the command line interface
    """
    odb_file = openOdb(input_filename)
    odb_filename = os.path.split(odb_file.path)[1].split(".")[0]  #.odb is removed
    odb_pathname = os.path.split(odb_file.path)[0]

    return odb_file, odb_filename, odb_pathname


#singleprecision_warning_given = False
def data_extractor(odb_file, mode="deformed"):
    """
    Extracts data from an odb file object (centroid (x,y), density, EVOL,
    pressure, nodal x,y pairs) and returns it for both the first and last frames
    """

    # Pull out the first and last frames from the output file. Note that this
    # assumes we're only interested in the last step, and that at the start of that
    # step, there's no deformation (i.e. at increment zero, the densities are all at
    # their original, given values)

    # Pull out the first frame of the model to get true, undeformed node
    # positions
    first_step_name = odb_file.steps.keys()[0]
    first_frame = odb_file.steps[first_step_name].frames[0]

    # Pull out data from the specified (or, by default, the last) step in the
    # model. If there's a specified step name but it doesn't exist in this
    # model, inform the user and quit.
    if default_step_name == "":
        step_name = odb_file.steps.keys()[-1]
    else:
        if default_step_name not in odb_file.steps.keys():
            print "ERROR: Default step does not exist in this model; output is from last frame instead. Please change the specified step in the Options section."
        step_name = default_step_name
    last_frame = odb_file.steps[step_name].frames[last_frame_number]

    # Get all of the elements in a set
    instance_names = odb_file.rootAssembly.instances.keys()
    this_instance = odb_file.rootAssembly.instances[instance_names[0]]
    all_elements = this_instance.elements #NOTE: all_elements[0] = element #1

    # Look up all of the nodes in the instance, and create a dictionary where
    # each Node object is identified by its label
    all_nodes = {}
    for node in this_instance.nodes:
        all_nodes[node.label] = node

    # Get all of the element volume data, and associate it with elements in a
    # dictionary (entries are {elemID:volume}), first for the first frame, and
    # then for the last frame
    initial_volumes = {}
    for item in first_frame.fieldOutputs['EVOL'].values:
        elemID = item.elementLabel
        initial_volumes[elemID] = item.data
    final_volumes = {}
    for item in last_frame.fieldOutputs['EVOL'].values:
        elemID = item.elementLabel
        final_volumes[elemID] = item.data

    # Check to make sure every element is accounted for
    if (len(all_elements) != len(initial_volumes)) or \
       (len(all_elements) != len(final_volumes)):
        print "WARNING: Not all element volumes accounted for. Does this model have inactive elements?"

    # For each section, go through and find the density of the material and the IDs
    # of the associated elements. Result is {elemID:density}. Then check that
    # we've got everything.
    initial_densities = {}
    num_sections = len(this_instance.sectionAssignments)
    for i in range(num_sections):
        these_elements = this_instance.sectionAssignments[i].region.elements
        this_section_name = this_instance.sectionAssignments[i].sectionName
        this_material = odb_file.sections[this_section_name].material
        this_density = odb_file.materials[this_material].density.table[0][0]

        for this_element in these_elements:
            initial_densities[this_element.label] = this_density

    if len(all_elements) != len(initial_densities):
        print "WARNING: Not all element densities accounted for. Does this model have inactive elements?"

    # Grab the pressures from the first and last frames, just like we did above
    # for volume (entries are {elemID:pressure}, and then check again for
    # accountablility, also like before
    initial_pressures = {}
    for item in first_frame.fieldOutputs['S'].values:
        elemID = item.elementLabel
        initial_pressures[elemID] = item.press
    final_pressures = {}
    for item in last_frame.fieldOutputs['S'].values:
        elemID = item.elementLabel
        final_pressures[elemID] = item.press

    if (len(all_elements) != len(initial_pressures)) or \
       (len(all_elements) != len(final_pressures)):
        print "WARNING: Not all element pressures accounted for. Does this model have inactive elements?"

    # To get the centroid, we're going to need displacement data, so first we'll
    # have to get those values in a way that doesn't care if we're using "U" or
    # "UT":
    try:
        final_displacements = last_frame.fieldOutputs['U'].values
    except KeyError:
        final_displacements = last_frame.fieldOutputs['UT'].values
    #for item in last_frame.fieldOutputs['U'].values:

    # Now we can go ahead and make another dictionary with
    # {nodeID:(nodal displacement value list)}. Note that we only need to do
    # this if we're outputting information from the last frame. We'll also have
    # to put this in an if/else clause so we can deal with 2D and 3D models
    # automatically
    node_displacement = {}
    singleprecision_warning_given = False
    for item in final_displacements:
        nodeID = item.nodeLabel

        if threedee_axi_mode:
            try:
                x_displacement, y_displacement, z_displacement = (item.dataDouble[0],
                                                                  item.dataDouble[1],
                                                                  item.dataDouble[2])
            except:
                #global singleprecision_warning_given
                if not singleprecision_warning_given:
                    if verbose_mode:
                        print "NOTICE: Model is single-precision!"
                    singleprecision_warning_given = True
                x_displacement, y_displacement, z_displacement = (item.data[0],
                                                                  item.data[1],
                                                                  item.data[2])

            node_displacement[nodeID] = (x_displacement, y_displacement, z_displacement)

        else:
            try:
                x_displacement, y_displacement = (item.dataDouble[0],
                                                  item.dataDouble[1])
            except:
                #global singleprecision_warning_given
                if not singleprecision_warning_given:
                    if verbose_mode:
                        print "NOTICE: Model is single-precision!"
                    singleprecision_warning_given = True
                x_displacement, y_displacement = (item.data[0],
                                                  item.data[1])

            node_displacement[nodeID] = (x_displacement, y_displacement)

    if len(all_nodes) != len(node_displacement):
        print "WARNING: Not all nodal displacements accounted for. Does this model have inactive elements?"

    # Now we have a few different sets (volume, pressure, density, and
    # displacement) all organized into dictionaries of elemID:value. The next
    # step is to go through every element in the model and figure out its
    # centroid position and nodal coordinates, both for the initial and the
    # final step.

    initial_centroids = {}
    final_centroids = {}
    initial_nodeinfo = {}
    final_nodeinfo = {}
    for this_element in all_elements:

        elemID = this_element.label

        # The "connectivity" of an element is which nodes it connects to, so we
        # can go through that list for each element and look at how the nodes
        # have moved to get our final values.

        # During this process, we'll account for "3D-axi" mode by simply
        # ignoring any nodes that don't lie on Z=0. Afterwars, we'll go through
        # our set of elements to remove any of these "axisymmetrized" elements
        # that are invalid (< 3 nodes).

        initial_x_coords = []
        initial_y_coords = []
        initial_z_coords = [] #Goes unused for 2D models
        valid_z_coords = 0 #Goes unused in 2D models
        valid_z_coord_cutoff = 1e-2 #m #Goes unused in 2D models
        final_x_coords = []
        final_y_coords = []
        final_z_coords = [] #Goes unused for 2D models
        initial_nodeinfo[elemID] = []
        final_nodeinfo[elemID] = []
        for node in this_element.connectivity:

            #3D axisymmetry mode:
            if threedee_axi_mode:
                # Structure this so that the very first thing we do is check to
                # see if the Z coordinate in question is on (or very near) the
                # Z=0 plane
                initial_z = all_nodes[node].coordinates[2]
                if -valid_z_coord_cutoff <= initial_z <= valid_z_coord_cutoff:
                    valid_z_coords = valid_z_coords + 1
                    initial_x = all_nodes[node].coordinates[0]
                    initial_y = all_nodes[node].coordinates[1]
                    #Already figured out initial_z above
                    initial_x_coords.append(initial_x)
                    initial_y_coords.append(initial_y)
                    initial_z_coords.append(initial_z)
                    initial_nodeinfo[elemID].append([initial_x,initial_y,initial_z])

                    x_displacement = node_displacement[node][0]
                    y_displacement = node_displacement[node][1]
                    z_displacement = node_displacement[node][2]

                    final_x = initial_x + x_displacement
                    final_y = initial_y + y_displacement
                    final_z = initial_z + z_displacement
                    final_x_coords.append(final_x)
                    final_y_coords.append(final_y)
                    final_z_coords.append(final_z)
                    final_nodeinfo[elemID].append([final_x,final_y,final_z])
                else:
                    continue

            #2D mode (default):
            else:
                initial_x = all_nodes[node].coordinates[0]
                initial_y = all_nodes[node].coordinates[1]
                initial_x_coords.append(initial_x)
                initial_y_coords.append(initial_y)
                initial_nodeinfo[elemID].append([initial_x,initial_y])

                x_displacement = node_displacement[node][0]
                y_displacement = node_displacement[node][1]

                final_x = initial_x + x_displacement
                final_y = initial_y + y_displacement
                final_x_coords.append(final_x)
                final_y_coords.append(final_y)
                final_nodeinfo[elemID].append([final_x,final_y])

        # If we're in 2D mode, OR if we're in 3D-axi mode with the element in
        # question having at least three nodes on (or within a very small distance
        # of) the Z=0 plane, then proceed with figuring out the centroid
        # positions (projected, in the case of 3D-axi mode)
        if (not threedee_axi_mode) or (threedee_axi_mode and (len(initial_z_coords) >= 3)):
            initial_centroids[elemID] = ((sum(initial_x_coords)/len(initial_x_coords)),
                                        (sum(initial_y_coords)/len(initial_y_coords)))
            final_centroids[elemID] = ((sum(final_x_coords)/len(final_x_coords)),
                                      (sum(final_y_coords)/len(final_y_coords)))

        # While we're at it, make sure we note that that volumes of these
        # elements are invalid. By calling it "na", we'll trigger special
        # behavior in `grav_anomaly.py` that will do the volume calculation for
        # us (after subdivision, if necessary)
        if (threedee_axi_mode and len(initial_z_coords) >= 3):
            initial_volumes[elemID] = "na"
            final_volumes[elemID] = "na"

        # If we're in 3D-axi mode but _don't_ have 3 elements on/near the Z=0
        # plane, though, expunge the element from all records as though it never
        # existed, 1984-style. Make it an un-element.
        elif (threedee_axi_mode and len(initial_z_coords) < 3):
            del initial_nodeinfo[elemID]
            del final_nodeinfo[elemID]
            del initial_volumes[elemID]
            del final_volumes[elemID]
            del initial_densities[elemID]
            #del final_densities[elemID] ###DEBUG: this raised an error
            del initial_pressures[elemID]
            del final_pressures[elemID]

    # Perform a correction on the density data, using the initial and final
    # volumes to figure out how much it's changed. Add these to a new set like
    # the others
    final_densities = {}
    for elemID in initial_densities.keys():

        ##TODO: FIGURE OUT HOW TO ACTUALLY SKIP UN-INSTANCED ELEMENTS!!!
        try:
            initial_volume = initial_volumes[elemID]
            final_volume = final_volumes[elemID]
            initial_density = initial_densities[elemID]
            final_densities[elemID] = initial_density * (initial_volume/final_volume)
        except KeyError:
            print "WARNING: Final density not being calculated for element %s!"%elemID
            continue

    # Finally, take all of this data we've amassed, and make a collection of
    # Elements that we can pass out of this function
    initial_elements = []
    final_elements = []

    for elemID in initial_centroids.keys():
        try:
            initial_elements.append(Element( elemID = elemID,
                                       centroid = initial_centroids[elemID],
                                          nodes = initial_nodeinfo[elemID],
                                         volume = initial_volumes[elemID],
                                        density = initial_densities[elemID],
                                       pressure = initial_pressures[elemID] ))
        except KeyError:
            print "WARNING: Data is not being written for element %s!"%elemID
            continue

        if mode == "undeformed":
            final_elements.append(Element( elemID = elemID,
                                         centroid = initial_centroids[elemID],
                                            nodes = initial_nodeinfo[elemID],
                                           volume = initial_volumes[elemID],
                                          density = final_densities[elemID],
                                         pressure = final_pressures[elemID] ))
        elif mode == "deformed":
            final_elements.append(Element( elemID = elemID,
                                         centroid = final_centroids[elemID],
                                            nodes = final_nodeinfo[elemID],
                                           volume = final_volumes[elemID],
                                          density = final_densities[elemID],
                                         pressure = final_pressures[elemID] ))

    return initial_elements, final_elements


def grav_outputter(elements,odb_pathname,odb_filename,output_suffix):
    """
    Writes out data from a given element set to a filename taken from the odb's
    filename and a file suffix & extension. This data is element ID, centroid X,
    Y, density, volume, and the nodal coordinates
    """
    # Create files for outputting the necessary data for gravity calculation
    output_fullfilename = "%s/%s%s"%(odb_pathname,odb_filename,output_suffix)
    output_file = open(output_fullfilename,'w')

    ###NEW:
    # Figure out which frame we're in
    if "ff" in output_suffix:
        which_frame = "First frame"
    elif last_frame_number != -1:
        which_frame = "Frame %s"%last_frame_number
    else:
        which_frame = "Last frame"

    # Headers
    output_file.write("%s, %s: %s\n"%(odb_filename, "current frame",
                      which_frame))
    output_file.write("%10s %16s %16s %16s %16s | %32s\n" \
                      %("element ID","centroid X","centroid Y",
                        "density","volume","node information"))
    # Dashes
    output_file.write("%10s %16s %16s %16s %16s-|-%32s\n"%\
                      ("-"*10,"-"*16,"-"*16,"-"*16,"-"*16,"-"*32))
    # Content
    for this_element in elements:
        output_file.write("%10i %16.4f %16.4f %16.4f %16.8g | %16s\n"\
                           %(this_element.elemID,
                             this_element.centroid_x,
                             this_element.centroid_y,
                             this_element.density,
                             this_element.volume,
                             this_element.nodes))

    # Truncate and close the file
    output_file.truncate()
    output_file.close()

    # Let the user know it's been done
    if verbose_mode:
        print "Element data for gravity calculation written to %s/%s%s" \
              %(odb_pathname, odb_filename, output_suffix)


def density_outputter(elements,odb_pathname,odb_filename,output_suffix):
    """
    Writes out density data from a given element set to a filename taken from
    the odb's filename and a file suffix & extension. This data is just the
    centroid's X and Y coordinates and the element's density.
    """
    # Output just the x, y, and density data to new files
    output_fullfilename = "%s/%s%s"%(odb_pathname,odb_filename,output_suffix)
    output_file = open(output_fullfilename,'w')

    # Don't print any headers, move right on to data
    for this_element in elements:
        output_file.write("%16.4f %16.4f %16.4f\n" \
                           %(this_element.centroid_x,
                             this_element.centroid_y,
                             this_element.density))

    # Truncate and close the file, then let the user know what happened
    output_file.truncate()
    output_file.close()
    if verbose_mode:
        print "Density data written to %s/%s%s" \
              %(odb_pathname, odb_filename, output_suffix)


def pressure_outputter(elements,odb_pathname,odb_filename,output_suffix):
    """
    Writes out pressure data from a given element set to a filename taken from
    the odb's filename and a file suffix & extension. This data includes the
    element's ID #, the centroid's X and Y coordinates, and the element's
    pressure.
    """
    # Output just the x, y, and density data to new files
    output_fullfilename = "%s/%s%s"%(odb_pathname,odb_filename,output_suffix)
    output_file = open(output_fullfilename,'w')

    # Don't print any headers, move right on to data
    for this_element in elements:
        #output_file.write("%%16.4f %16.4f %16.4f\n" \
        #                     this_element.centroid_x,
        #                     this_element.centroid_y,
        #                     this_element.pressure))
        output_file.write("%8g %16.4f %16.4f %16.4f\n" \
                           %(this_element.elemID,
                             this_element.centroid_x,
                             this_element.centroid_y,
                             this_element.pressure))

    # Truncate and close the file, then let the user know what happened
    output_file.truncate()
    output_file.close()
    if verbose_mode:
        print "Pressure data written to %s/%s%s" \
              %(odb_pathname, odb_filename, output_suffix)


######## Command-Line Implementation ###########################################

if __name__ == '__main__':

    # Abaqus CAE calls itself a command-line environment, so this behavior is
    # what happens when you run the script in CAE

    odb_file, odb_pathname, odb_filename = get_odb_from_cae()
    firstframe_data, lastframe_data = data_extractor(odb_file)

    if output_firstframe:
        grav_outputter(firstframe_data, odb_pathname, odb_filename,
                  output_suffix_firstframe)
        if output_density:
            density_outputter(firstframe_data, odb_pathname, odb_filename,
                              output_suffix_firstframedensity)
        if output_pressure:
            pressure_outputter(firstframe_data, odb_pathname, odb_filename,
                               output_suffix_firstframepressure)

    if output_lastframe:
        grav_outputter(lastframe_data, odb_pathname, odb_filename,
                  output_suffix_lastframe)
        if output_density:
            density_outputter(lastframe_data, odb_pathname, odb_filename,
                              output_suffix_lastframedensity)
        if output_pressure:
            pressure_outputter(lastframe_data, odb_pathname, odb_filename,
                               output_suffix_lastframepressure)

        if output_undeformed:
            lastframe_nodefdata = data_extractor(odb_file, mode="undeformed")[1]
            if output_density:
                density_outputter(lastframe_nodefdata, odb_pathname, odb_filename,
                                  output_suffix_nodeflastframedensity)
            if output_pressure:
                pressure_outputter(lastframe_nodefdata, odb_pathname, odb_filename,
                                   output_suffix_nodeflastframepressure)
