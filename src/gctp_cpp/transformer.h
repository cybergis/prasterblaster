#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "projection.h"
#include "coordinate.h"

//! Transformer class
/*! This class provides an interface that auto-couples two
   Projections together and transforms coordinates.
   */
class Transformer {
public:
   //! Constructor
   /*! This constructor takes no parameters and simply nulls the
      private attributes. The decision to not create a constructor
      that could take all parameters at once was to force the coder
      to use more readable mechanics that those employed by the
      original GCTP.
      */
   Transformer();

   //! Destructor
   /*! The simply deletes the dynamically allocated projections
      if they exist.
      */
   ~Transformer();

   /*! Sets the input projection according to the projection code passed.
      Setting the finer details such as the units parameter can be done by
      calling transformer->input()->setUnits(...).
      */
   bool setInput( ProjCode projectionCode );

   /*! Sets the input projection with all parameters set.
      The parameters can be later adjusted by calling transformer->input()->setUnits(...).
      \param projectionCode GCTP defined enumeration for projections
      \param gctpParameters[] The 15 parameters that define the variations of a projection
      \param units GCTP defined enumeration for units
      \param datum GCTP defined enumeration for datum
      \param spheroid GCTP defined enumeration for spheroid
      */
   bool setInput( ProjCode projectionCode, double gctpParameters[15], ProjUnit units, ProjDatum datum );

   //! Sets the input projection by copying all values out of the argument
   bool setInput( Projection &in );

   //! Provides access to the input projection
   Projection* input();

   /*! Sets the output projection according to the projection code passed.
      Setting the finer details such as the units parameter can be done by
      calling transformer->output()->setUnits(...).
      */
   bool setOutput( ProjCode projectionCode );

   /*! Sets the output projection with all parameters set.
      The parameters can be later adjusted by calling transformer->output()->setUnits(...).
      \param projectionCode GCTP defined enumeration for projections
      \param gctpParameters[] The 15 parameters that define the variations of a projection
      \param units GCTP defined enumeration for units
      \param datum GCTP defined enumeration for datum
      \param spheroid GCTP defined enumeration for spheroid
      */
   bool setOutput( ProjCode projectionCode, double gctpParameters[15], ProjUnit units, ProjDatum datum );
   
   //! Sets the output projection by copying all values out of the argument
   bool setOutput( Projection &out );

   //! Provides access to the output projection
   Projection* output();

   /*! Pass a Coordinate from the input projection and its attributes
      will be transformed into the output projection.
      Precondition:  Input and output have been set.
                     io_coord->units are the same as the input units.
      Postcondition: io_coord->(x and y) will be changed to the output.
                     io_coord->units will be the same as output units.
      */
   void transform( Coordinate* io_coord );

   /*! Pass a Coordinate from the input projection and its attributes
      will be transformed to represent the geographic lat/lon of that coordinate.
      Precondition:  Input has been set.
                     io_coord->units are the same as the input units.
      Postcondition: io_coord->(x and y) will be changed to geographic lat/lon.
                     io_coord->units will be DEGREE.
      */
   void transformInverse( Coordinate* io_coord );

   /*! Pass a Coordinate geographic lat/lon and its attributes
      will be transformed into the output projection.
      Precondition:  Output has been set.
                     io_coord->units are DEGREE.
      Postcondition: io_coord->(x and y) will be changed to the output.
                     io_coord->units will be the same as output units.
      */
   void transformForward( Coordinate* io_coord );

   //! Returns true if an error occured in the most recent transformation
   bool errored();

   /*! A giant, simple, switch case that returns a dynamically allocated
      Projection* that points to a Projection subclass whose type is
      determined by the projectionCode passed.
      Note: It will need to be deleted in order to prevent memory leaks.
      */
   static Projection* convertProjection( ProjCode projectionCode );

private:
   //! Points to the input projection using polymorphism
   Projection* m_inProj;

   //! Points to the output projection using polymorphism
   Projection* m_outProj;

   //! Tells whether the most recent transformation errored or not
   bool m_errored;
};

#endif
