My Main Page                         {#mainpage}
============

libRasterBlaster - Parallel map reprojection library            
============
    librasterblaster is a library designed to allow the quick creation
    of parallel raster reprojection programs. Hopefully this will
    allow further research into I/O, partitioning, load balancing and
    other HPC concerns with regard to map reprojection.

    librasterblaster is made up of components that implement different
    parts of a parallel reprojection program. The components are
    designed so that they can be easily replaced with a user-provided
    component, hopefully allowing experimentation with new techniques.

    There is a sample implementation of parallel reprojection in the
    prasterblaster-pio.cc file.

    You  will probably want  to create  your own  reprojection program
    with some  differences from  the demo. All  librasterblaster based
    programs will have a similar structure:
    
    - Parse input parameters
    - Create the output raster
    - Partition the output raster
    - For each output partition, 
      * Find the matching input partition
      * Read the matching input area
      * Call the ReprojectChunk function to perform resampling


    - Parse input parameters 

        This is ordinarily a simple task but an example class is included
    called: librasterblaster::Configuration .
    

    Main Headers:
    - projectedraster.h
    - rasterchunk.h
    - reprojection_tools.h
    - configuration.h
