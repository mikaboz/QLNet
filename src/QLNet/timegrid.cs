/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
  
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
   //! time grid class
   /*! \todo what was the rationale for limiting the grid to
             positive times? Investigate and see whether we
             can use it for negative ones as well.
   */

   public class TimeGrid : List<double>
   {
      public double First()
      {
         return this[0];
      }
      public double Last()
      {
         return this[this.Count];
      }
      public double dt(int i)
      {
         return this[i + 1] - this[i];
      }
      public List<double> mandatoryTimes()
      {
         return this;
      }
      public int size()
      {
         return this.Count;
      }
      public TimeGrid(double begin, double end, int steps) :
         base(steps + 1)
      {
         double dt = (end - begin) / steps;
         this[0] = begin;
         for (int i = 1; i < Count; i++)
            this[i] = this[i - 1] + dt;
      }
      public TimeGrid(double end, int steps) :
         this(0.0, end, steps)
      { }
      public TimeGrid(List<double> times) : base(times) { }
      public TimeGrid(List<double> times, int uselessItem,int steps) :
         base(steps)
      {

         // The resulting timegrid have points at times listed in the input
         // list. Between these points, there are inner-points which are
         // regularly spaced.
         double last = times.Last();
         double dtMax = 0;

         if (steps == 0)
         {
            List<double> diff = new List<double>();
         }
         else
         {
            dtMax = last / steps;
         }

         double periodBegin = 0.0;
         this.Add(periodBegin);

         for (int k = 0; k < times.Count; k++)
         {
            double dt = 0;
            double periodEnd = times[k];
            if (periodEnd.IsNotEqual(0.0))
            {
               // the nearest integer
               int nSteps = (int)((periodEnd - periodBegin) / dtMax + 0.5);
               // at least one time step!
               nSteps = (nSteps != 0 ? nSteps : 1);
               dt = (periodEnd - periodBegin) / nSteps;
               for (int n = 1; n <= nSteps; ++n)
               {
                  this.Add(periodBegin + n * dt);
               }
            }
            periodBegin = periodEnd;
         }
      }

      // Time grid interface
      //! returns the index i such that grid[i] = t
      public int index(double t)
      {
         int i = closestIndex(t);
         if (Utils.close(t, this[i]))
         {
            return i;
         }
         Utils.QL_REQUIRE(t >= this.First(), () =>
             "using inadequate time grid: all nodes are later than the required time t = "
             + t + " (earliest node is t1 = " + this.First() + ")");
         Utils.QL_REQUIRE(t <= this.Last(), () =>
             "using inadequate time grid: all nodes are earlier than the required time t = "
             + t + " (latest node is t1 = " + this.Last() + ")");
         int j, k;
         if (t > this[i])
         {
            j = i;
            k = i + 1;
         }
         else
         {
            j = i - 1;
            k = i;
         }
         Utils.QL_FAIL("using inadequate time grid: the nodes closest to the required time t = "
                             + t + " are t1 = " + this[j] + " and t2 = " + this[k]);
         return 0;
      }

      //! returns the index i such that grid[i] is closest to t
      public int closestIndex(double t)
      {
         int result = this.BinarySearch(t);
         if (result < 0)
            //Lower_bound is a version of binary search: it attempts to find the element value in an ordered range [first, last)
            // [1]. Specifically, it returns the first position where value could be inserted without violating the ordering. 
            // [2] The first version of lower_bound uses operator< for comparison, and the second uses the function object comp.
            // lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value. 
            result = ~result;

         if (result == 0)
         {
            return 0;
         }
         if (result == this.Count)
         {
            return this.Count - 1;
         }
         double dt1 = this[result] - t;
         double dt2 = t - this[result - 1];
         if (dt1 < dt2)
            return result;
         return result - 1;
      }

      //! returns the time on the grid closest to the given t
      public double closestTime(double t)
      {
         return this[closestIndex(t)];
      }
   }

   public class OldTimeGrid
   {
      protected List<double> times_;
      public List<double> Times() { return times_; }

      public double this[int i]
      {
         get { return times_[i]; }
      }

      private List<double> dt_;

      public double dt(int i)
      {
         return dt_[i];
      }

      private List<double> mandatoryTimes_;

      public List<double> mandatoryTimes()
      {
         return mandatoryTimes_;
      }

      protected OldTimeGrid() { }
      public OldTimeGrid(double end, int steps)
      {
         // We seem to assume that the grid begins at 0.
         // Let's enforce the assumption for the time being
         // (even though I'm not sure that I agree.)
         Utils.QL_REQUIRE(end > 0.0, () => "negative times not allowed");
         double dt = end / steps;
         times_ = new List<double>(steps + 1);
         for (int i = 0; i <= steps; i++)
            times_.Add(dt * i);

         mandatoryTimes_ = new InitializedList<double>(1);
         mandatoryTimes_[0] = end;

         dt_ = new InitializedList<double>(steps, dt);
      }

      public OldTimeGrid(List<double> times, int offset)
      {
         //not really finished bu run well for actals tests
         mandatoryTimes_ = times.GetRange(0, offset);
         mandatoryTimes_.Sort();

         Utils.QL_REQUIRE(mandatoryTimes_[0] >= 0.0, () => "negative times not allowed");

         for (int i = 0; i < mandatoryTimes_.Count - 1; ++i)
         {
            if (Utils.close_enough(mandatoryTimes_[i], mandatoryTimes_[i + 1]))
            {
               mandatoryTimes_.RemoveAt(i);
               i--;
            }
         }

         times_ = new List<double>(mandatoryTimes_);

         if (mandatoryTimes_[0] > 0.0)
            times_.Insert(0, 0.0);

         var dt = times_.Zip(times_.Skip(1), (x, y) => y - x);
         dt_ = dt.ToList();
      }

      public OldTimeGrid(List<double> times, int offset, int steps)
      {
         //not really finished bu run well for actals tests
         mandatoryTimes_ = times.GetRange(0, offset);
         mandatoryTimes_.Sort();

         Utils.QL_REQUIRE(mandatoryTimes_[0] >= 0.0, () => "negative times not allowed");

         for (int i = 0; i < mandatoryTimes_.Count - 1; ++i)
         {
            if (Utils.close_enough(mandatoryTimes_[i], mandatoryTimes_[i + 1]))
            {
               mandatoryTimes_.RemoveAt(i);
               i--;
            }
         }

         // The resulting timegrid have points at times listed in the input
         // list. Between these points, there are inner-points which are
         // regularly spaced.
         times_ = new List<double>(steps);
         dt_ = new List<double>(steps);
         double last = mandatoryTimes_.Last();
         double dtMax = 0;

         if (steps == 0)
         {
            List<double> diff = new List<double>();
         }
         else
         {
            dtMax = last / steps;
         }

         double periodBegin = 0.0;
         times_.Add(periodBegin);

         for (int k = 0; k < mandatoryTimes_.Count; k++)
         {
            double dt = 0;
            double periodEnd = mandatoryTimes_[k];
            if (periodEnd.IsNotEqual(0.0))
            {
               // the nearest integer
               int nSteps = (int) ((periodEnd - periodBegin) / dtMax + 0.5);
               // at least one time step!
               nSteps = (nSteps != 0 ? nSteps : 1);
               dt = (periodEnd - periodBegin) / nSteps;
               for (int n = 1; n <= nSteps; ++n)
               {
                  times_.Add(periodBegin + n * dt);
                  dt_.Add(dt);
               }
            }
            periodBegin = periodEnd;
         }
      }

      // Time grid interface
      //! returns the index i such that grid[i] = t
      public int index(double t)
      {
         int i = closestIndex(t);
         if (Utils.close(t, times_[i]))
         {
            return i;
         }
         Utils.QL_REQUIRE(t >= times_.First(),()=> 
            "using inadequate time grid: all nodes are later than the required time t = "
            + t + " (earliest node is t1 = " + times_.First() + ")");
         Utils.QL_REQUIRE(t <= times_.Last(),()=> 
            "using inadequate time grid: all nodes are earlier than the required time t = " 
            + t + " (latest node is t1 = " + times_.Last() + ")");
         int j, k;
         if (t > times_[i])
         {
            j = i;
            k = i + 1;
         }
         else
         {
            j = i - 1;
            k = i;
         }
         Utils.QL_FAIL("using inadequate time grid: the nodes closest to the required time t = "
                             + t + " are t1 = " + times_[j] + " and t2 = " + times_[k]);
         return 0;
      }

      //! returns the index i such that grid[i] is closest to t
      public int closestIndex(double t)
      {
         int result = times_.BinarySearch(t);
         if (result < 0)
            //Lower_bound is a version of binary search: it attempts to find the element value in an ordered range [first, last)
            // [1]. Specifically, it returns the first position where value could be inserted without violating the ordering. 
            // [2] The first version of lower_bound uses operator< for comparison, and the second uses the function object comp.
            // lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value. 
            result = ~result;

         if (result == 0)
         {
            return 0;
         }
         if (result == times_.Count)
         {
            return times_.Count - 1;
         }
         double dt1 = times_[result] - t;
         double dt2 = t - times_[result - 1];
         if (dt1 < dt2)
            return result;
         return result - 1;
      }

      //! returns the time on the grid closest to the given t
      public double closestTime(double t)
      {
         return times_[closestIndex(t)];
      }

      public bool empty()
      {
         return times_.Count == 0;

      }

      public int size()
      {
         return times_.Count;
      }

      public double First()
      {
         return times_.First();
      }

      public double Last()
      {
         return times_.Last();
      }
   }
  

   /*
   public class UniformTimeGrid : TimeGrid
   {
      public UniformTimeGrid(double end, int steps) :
         base(end,steps)
      {}
   }
   */
   /*
   public class SequenceTimeGrid : TimeGrid
   {
      public List<TimeGrid> timeGrids_;

      public TimeGrid TimeGrid(int i)
      {
         return timeGrids_[i];
      }

      public int SequenceSize() { return timeGrids_.Count + 1; }

      public SequenceTimeGrid(List<TimeGrid> timeGrids, double firstPeriodTime, List<double> phasis) :
         base()
      {
         Console.WriteLine("timegridsCount: " + timeGrids.Count);
         Console.ReadLine();
         timeGrids_ = timeGrids;
         Utils.QL_REQUIRE(phasis.Count == timeGrids.Count - 1, () => "phasis vector must size be timeGrids.size - 1");
         times_ = new List<double>();
         times_.Add(0);
         TimeGrid firstGrid = timeGrids_.First();
         for (int i = 0; i < firstGrid.size(); i++)
         {
            firstGrid.Times()[i] += firstPeriodTime;
            times_.Add(firstGrid[i]);
         }
         for (int j = 1; j<timeGrids_.Count;j++)
         {
            double lastValue = times_.Last();
            for (int i = 0; i < timeGrids_[j].size(); i++)
            {
               timeGrids_[j].Times()[i+1] += lastValue + phasis[j-1];
               if (i!= 0) times_.Add(timeGrids_[j][i]);
            }
         }
         foreach (double d in times_)
            Console.WriteLine(d);
         Console.ReadLine();
      }
   }
   */
}
