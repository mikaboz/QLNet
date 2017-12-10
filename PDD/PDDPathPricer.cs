using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using QLNet;
namespace PDD
{
   public static class PathExtension
   {
      public static Path getRange(this Path path, int from, int to)
      {
         Vector v = new Vector(path.values().GetRange(from, to));
         List<double> values = path.timeGrid().Times().GetRange(from, to).ToList();
         TimeGrid tg = new TimeGrid(values, values.Count);
         return new Path(tg, v);
      }
      public static Path getRange(this Path path, double first, double second)
      {
         int from = path.timeGrid().Times().IndexOf(first);
         int to = path.timeGrid().Times().IndexOf(second);
         return path.getRange(from, to);
      }
   }

   public class PDDPathPricer : PathPricer<IPath>
   {
      PDD.Payoff payoff_;
      double discount_;
      public PDDPathPricer(PDD.Payoff payoff, double discount)
      {
         payoff_ = payoff;
         discount_ = discount;
      }
      private bool TerminalCondition(double price)
      {
         return payoff_.terminalCondition(price); 
      }
      private bool MaximumCondition(Path path)
      {
         PDDTimeGrid timeGrid = path.timeGrid() as PDDTimeGrid;
         Utils.QL_REQUIRE(timeGrid != null, () => "time grid must be PDDTimeGrid, that is specific grid for PDD pricing");

         double acquisitionPrice = payoff_.acquisitionValue_;
         double max = path.values().Skip(1).Max();
         return max < acquisitionPrice;
      }
      public double value(IPath path)
      {
         Path path_ = path as Path;
         PDDTimeGrid timeGrid = path_.timeGrid() as PDDTimeGrid;
         Utils.QL_REQUIRE(timeGrid != null, () => "timeGrid must be PDDTimeGrid for PDD Monte Carlo pricing;");
         Console.WriteLine("path lenght: " + path.length());
         double forwardPrice;
         PDD.Payoff tempPayoff = new PDD.Payoff(payoff_.percent_, payoff_.acquisitionValue_);
         int stepsPerConditionalPeriod_ = timeGrid.stepsPerConditionalPeriod_;
         int index = 0;
         double forwardValue = 0;

         while (index < path.length())
         {
            Vector conditionalPeriodValues;
            if (timeGrid.startWithConditionalPeriod_ == true)
            {
               forwardPrice = path_[stepsPerConditionalPeriod_];
               conditionalPeriodValues = new Vector(path_.values().GetRange(index, stepsPerConditionalPeriod_));
               index += stepsPerConditionalPeriod_;
            }
            else
            {
               forwardPrice = path_[stepsPerConditionalPeriod_ + 1];
               conditionalPeriodValues = new Vector(path_.values().GetRange(index+1, stepsPerConditionalPeriod_));
               index += stepsPerConditionalPeriod_ +1;
            }
            if (tempPayoff.terminalCondition(forwardPrice))
            {
               double max = conditionalPeriodValues.Max();
               if (max < tempPayoff.acquisitionValue_)
               {
                  forwardValue = tempPayoff.value(forwardPrice);
                  tempPayoff.acquisitionValue_ -= forwardValue;
               }
            }
            Console.WriteLine("index: " + index);
         }
         return forwardValue * discount_;
      }
   }


   public class MaximumGenerator
   {
      public static double generate(double uniform, double xt_, double xT_, double sigmaxt_, double h_)
      {
         double diff = xT_ - xt_;
         return 0.5 * (xt_ + xT_ + Math.Sqrt(diff * diff + 2 * sigmaxt_ * sigmaxt_ * h_ * Math.Log(uniform)));
      }
   }


   public class PDDTimeGrid : TimeGrid
   {
      public bool startWithConditionalPeriod_;
      public int stepsPerConditionalPeriod_;
      public PDDTimeGrid(List<double> times, int stepsPerConditionalPeriod, bool startWithConditionalPeriod = false) :
         base(times, times.Count)
      {
         startWithConditionalPeriod_ = startWithConditionalPeriod;
         stepsPerConditionalPeriod_ = stepsPerConditionalPeriod;
         Console.WriteLine("timeGrid");
         foreach (double d in times_)
            Console.WriteLine(d);
         Console.ReadLine();
      }
   }

}
