libRasterBlaster - Parallel map reprojection library            {#mainpage}
============

    librasterblaster is a library designed to allow the quick creation
    of parallel raster reprojection programs. Hopefully this will
    allow further research into I/O, partitioning, load balancing and
    other HPC concerns with regard to map reprojection.

    For a detailed description of librasterblaster read the \link overview \endlink .

    Compiling
    ---------

    Run "./configure" in the librasterblaster source tree. If the
    configure script doesn't exist, run "autoreconf -iv". You may have
    to use the "--with-gdal-incdir" and "--with-gdal-libdir" options
    to tell autoconf where it can find GDAL. 

    Next, run "make". If everything compiles correctly you should have
    a \link prasterblasterpio \endlink binary.

    To build the documentation run "doxygen" in the source tree. You
    will need a doxygen version >= 1.8.0 for everything to
    work. Doxygen will create "html" and "latex" directories with the
    generated documentation.

    libRasterBlaster 
    ----------------

    librasterblaster is made up of components that implement different
    parts of a parallel reprojection program. The components are
    designed so that they can be easily replaced with a user-provided
    component, hopefully allowing experimentation with new techniques.

    The sample implementation, demonstrating the use of
    librasterblaster is called \link prasterblasterpio \endlink.

    You will probably want to create your own reprojection program
    with some differences from the demo program. However most
    librasterblaster based programs will have a similar structure:
    
    - Parse input parameters
    - Calculate output raster area
    - Create the output raster
    - Partition the output raster
    - For each output partition, 
      * Find the matching input partition
      * Read the matching input area
      * Call the ReprojectChunk function to perform resampling

    libRasterBlaster components
    ---------------------------

    The main components of libRasterBlaster map closely to the above
    steps.
    
    ### Parse input parameters

    This is ordinarily a simple task but an example class is included
    called: librasterblaster::Configuration.

    Parsing the spatial reference systems can be handled by the raster
    creation functions described below.

    ### Calculate output raster area
    
    Depending on the projection and pixel size selected, the output
    raster could be a very different size that the input, therefore we
    must calculate the size of the output space.
    
    The size of the output raster will depend on the size of the input
    raster, the size of the output raster's pixels, and the
    projections involved.

    For many use cases you should only have to call
    librasterblaster::CreateOutputRaster. It will create an empty
    raster of the correct size.
        
    It works by calling librasterblaster::ProjectedMinbox to find the
    size of the projected output space and dividing by the size of the
    output pixels to find the correct number of rows and columns in
    the output projection.

    ### Create output raster

    The previously mentioned librasterblaster::CreateOutputRaster can
    be used to create an empty, correctly sized, output raster.

    If you calculate the parameters of the output raster manually you
    use the librasterblaster::ProjectedRaster::CreateRaster functions.
    
    ### Partition the output raster
    
    It's expected that this would be a common function that the user
    will implement themselves. Different partitioning techniques could
    allow for different I/O and load balancing strategies.

    Two different partitioning functions are provided by
    librasterblaster: librasterblaster::RowPartition and
    librasterblaster::PartitionBySize.

    librasterblaster::RowPartition partitions a raster space into
    partitions smaller than the given size. All of the partitions are
    rows, collections of adjoining rows, or parts of a single row.

    It's main purpose was to support the demo implementations sptw
    parallel writer.

    librasterblaster::PartitionBySize uses a simple quad tree to
    partition a raster space into blocks smaller than a given size.

    ### Match output partition to input partition
    
    ### Read input partition area

    ### Call ReprojectChunk

    Main Headers:
    - projectedraster.h
    - rasterchunk.h
    - reprojection_tools.h
    - configuration.h
