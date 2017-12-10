/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2008-2016 Andrea Maggiulli (a.maggiulli@gmail.com)
  
 This file is part of QLNet Project https://github.com/amaggiulli/qlnet

 QLNet is free software: you can redistribute it and/or modify it
 under the terms of the QLNet license.  You should have received a
 copy of the license along with this program; if not, license is  
 available online at <http://qlnet.sourceforge.net/License.html>.
  
 QLNet is a based on QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 The QuantLib license is available online at http://quantlib.org/license.shtml.
 
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

using System;
using System.Collections.Generic;
using System.Linq;

namespace QLNet
{
   //! General-purpose Monte Carlo model for path samples
   /*! The template arguments of this class correspond to available
       policies for the particular model to be instantiated---i.e.,
       whether it is single- or multi-asset, or whether it should use
       pseudo-random or low-discrepancy numbers for path
       generation. Such decisions are grouped in trait classes so as
       to be orthogonal---see mctraits for examples.

       The constructor accepts two safe references, i.e. two smart
       pointers, one to a path generator and the other to a path
       pricer.  In case of control variate technique the user should
       provide the additional control option, namely the option path
       pricer and the option value.

       \ingroup mcarlo
   */
   public class MonteCarloModel<MC, RNG, S> where S : IGeneralStatistics
   {
      protected MonteCarloModel(IPathGenerator<IRNG> pathGenerator, bool antitheticVariate)
      {
         pathGenerator_ = pathGenerator;
         isAntitheticVariate_ = antitheticVariate;
      }
      public MonteCarloModel( IPathGenerator<IRNG> pathGenerator, PathPricer<IPath> pathPricer, 
         S sampleAccumulator,bool antitheticVariate, PathPricer<IPath> cvPathPricer = null, 
         double cvOptionValue = 0,IPathGenerator<IRNG> cvPathGenerator = null)
      {
         pathGenerator_ = pathGenerator;
         pathPricer_ = pathPricer;
         sampleAccumulator_ = sampleAccumulator;
         isAntitheticVariate_ = antitheticVariate;
         cvPathPricer_ = cvPathPricer;
         cvOptionValue_ = cvOptionValue;
         cvPathGenerator_ = cvPathGenerator;
         if ( cvPathPricer_ == null )
            isControlVariate_ = false;
         else
            isControlVariate_ = true;
      }

      public void addSamples( int samples )
      {
         for ( int j = 1; j <= samples; j++ )
         {

            Sample<IPath> path = pathGenerator_.next();
            double price = pathPricer_.value( path.value );

            if ( isControlVariate_ )
            {
               if ( cvPathGenerator_ == null )
               {
                  price += cvOptionValue_ - cvPathPricer_.value( path.value );
               }
               else
               {
                  Sample<IPath> cvPath = cvPathGenerator_.next();
                  price += cvOptionValue_ - cvPathPricer_.value( cvPath.value );
               }
            }

            if ( isAntitheticVariate_ )
            {
               path = pathGenerator_.antithetic();
               double price2 = pathPricer_.value( path.value );
               if ( isControlVariate_ )
               {
                  if ( cvPathGenerator_ == null )
                     price2 += cvOptionValue_ - cvPathPricer_.value( path.value );
                  else
                  {
                     Sample<IPath> cvPath = cvPathGenerator_.antithetic();
                     price2 += cvOptionValue_ - cvPathPricer_.value( cvPath.value );
                  }
               }

               sampleAccumulator_.add( ( price + price2 ) / 2.0, path.weight );
            }
            else
            {
               sampleAccumulator_.add( price, path.weight );
            }
         }
      }
      public S sampleAccumulator() { return sampleAccumulator_; }

      private IPathGenerator<IRNG> pathGenerator_;
      private PathPricer<IPath> pathPricer_;
      private S sampleAccumulator_;
      private bool isAntitheticVariate_;
      private PathPricer<IPath> cvPathPricer_;
      private double cvOptionValue_;
      private bool isControlVariate_;
      private IPathGenerator<IRNG> cvPathGenerator_;

   }

   /*
   public class MultiMonteCarloModel<MC, RNG, S> : MonteCarloModel<MC,RNG,S>where S : IGeneralStatistics
   {
      public MultiMonteCarloModel(IPathGenerator<IRNG> pathGenerator, List<PathPricer<IPath>> pathPricer,
         List<S> sampleAccumulator, bool antitheticVariate) :
         base(pathGenerator,antitheticVariate)
      {
         Utils.QL_REQUIRE(pathPricer.Count == sampleAccumulator.Count, () => "path pricer and sampleAccumulator must have the same size");
         pathPricer_ = pathPricer;
         sampleAccumulator_ = sampleAccumulator;
      }

      public new void addSamples(int samples)
      {
         for (int j = 1; j <= samples; j++)
         {

            Sample<IPath> path = pathGenerator_.next();
            List<double> price = new InitializedList<double>(pathPricer_.Count);
            for (int i= 0; i<price.Count;i++)
               price[i] = pathPricer_[i].value(path.value);

            if (isAntitheticVariate_)
            {
               path = pathGenerator_.antithetic();
               List<double> price2 = new InitializedList<double>(pathPricer_.Count);
               for (int i = 0; i<price2.Count;i++)
                  price2[i] = pathPricer_[i].value(path.value);

               for (int i = 0; i < sampleAccumulator_.Count; i++)
                  sampleAccumulator_[i].add((price[i] + price2[i]) / 2.0, path.weight);
            }
            else
            {
               for (int i = 0; i < sampleAccumulator_.Count; i++)
                  sampleAccumulator_[i].add(price[i], path.weight);
            }
         }
      }
      public new List<S> sampleAccumulator() { return sampleAccumulator_; }

      private IPathGenerator<IRNG> pathGenerator_;
      private List<PathPricer<IPath>> pathPricer_;
      private List<S> sampleAccumulator_;
      private bool isAntitheticVariate_;
   }
   */

}
