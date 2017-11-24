//  Copyright (C) 2008-2016 Andrea Maggiulli (a.maggiulli@gmail.com)
//  
//  This file is part of QLNet Project https://github.com/amaggiulli/qlnet
//  QLNet is free software: you can redistribute it and/or modify it
//  under the terms of the QLNet license.  You should have received a
//  copy of the license along with this program; if not, license is  
//  available online at <http://qlnet.sourceforge.net/License.html>.
//   
//  QLNet is a based on QuantLib, a free-software/open-source library
//  for financial quantitative analysts and developers - http://quantlib.org/
//  The QuantLib license is available online at http://quantlib.org/license.shtml.
//  
//  This program is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE.  See the license for more details.

using System;
using System.Collections.Generic;

namespace QLNet
{
   public class ProjectedConstraint : Constraint
   {
      public enum Mode { Inclusive, Exclusive }
      private class Impl : IConstraint 
      {
         public Impl( Constraint constraint,
                      Vector parameterValues,
                      List<bool>fixParameters, Mode mode)

         {
            constraint_ = constraint;
            projection_ = new Projection(parameterValues, fixParameters);
            mode_ = mode;
         }

         public Impl( Constraint constraint, Projection projection,Mode mode)
         {
            constraint_ = constraint;
            projection_ = projection;
            mode_ = mode;
         }

         public Vector ProjectAccordingToMode(Vector parameters)
         {
            Vector projectedVector;
            switch (mode_)
            {
               case Mode.Inclusive:
                  projectedVector = projection_.include(parameters);
                  break;
               case Mode.Exclusive:
                  projectedVector = projection_.project(parameters);
                  break;
               default:
                  throw new NotSupportedException("Not supported Mode of projection");
            }
            return projectedVector;
         }
         public bool test(Vector parameters) 
         {
            return constraint_.test(ProjectAccordingToMode(parameters));
         }
            
         public Vector upperBound(Vector parameters) 
         {
            return constraint_.upperBound(ProjectAccordingToMode(parameters));
         }
            
         public Vector lowerBound(Vector parameters) 
         {
            return constraint_.lowerBound(ProjectAccordingToMode(parameters));
         }

          private Mode mode_;
          private Constraint constraint_;
          private Projection projection_;
      }

      public ProjectedConstraint( Constraint constraint,
                                  Vector parameterValues,
                                  List<bool> fixParameters,Mode mode = Mode.Inclusive)
         : base( new Impl(constraint, parameterValues,fixParameters,mode)) 
      {}

      public ProjectedConstraint( Constraint constraint, Projection projection, Mode mode = Mode.Inclusive)
            : base(new Impl(constraint, projection,mode)) 
      {}
   }

   public class ProjectedIndividualConstraint : ProjectedConstraint
   {
      private static List<bool> GenerateList(Vector parameterValues, int toTestIndex)
      {
         if (toTestIndex >= parameterValues.Count)
            throw new NotSupportedException("index to test must be include in vector size");
         List<bool> fixedParameters = new InitializedList<bool>(parameterValues.Count, true)
         {
            [toTestIndex] = false
         };
         return fixedParameters;
      }
      public ProjectedIndividualConstraint(Constraint constraint,
                                  Vector parameterValues,
                                  int toTestIndex)
         : base(constraint,parameterValues,GenerateList(parameterValues,toTestIndex),Mode.Exclusive)
      { }
   }
}
