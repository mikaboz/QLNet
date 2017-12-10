using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using QLNet;
namespace PDD
{
   public class PDD : OneAssetOption
   {
      public class Payoff : QLNet.Payoff
      {
         public double percent_;
         public double acquisitionValue_;
         public Payoff(double percent, double acquisitionValue)
         {
            Utils.QL_REQUIRE(percent >= 0 && percent <= 1, () => "percent must belong to [0,1]");
            Utils.QL_REQUIRE(acquisitionValue >= 0, () => "acquisition value must be positive");
            percent_ = percent;
            acquisitionValue_ = acquisitionValue;
         }
         public bool terminalCondition(double forwardPrice)
         {
            return forwardPrice < percent_ * acquisitionValue_;
         }
         public override double value(double price)
         {
            return terminalCondition(price) ? acquisitionValue_ - price : 0.0;
         }
      }
      
      public Period conditionalPeriod_;
      public Period provFrequency_;
      public PDD(double acquisitionValue, double percent, Date maturityDate, Period conditionalPeriod = null, Period provFrequency = null) :
         base(new Payoff(percent,acquisitionValue), new EuropeanExercise(maturityDate))
      {
         provFrequency_ = provFrequency ?? new Period(1,TimeUnit.Years);
         conditionalPeriod_ = conditionalPeriod_ ?? new Period(6, TimeUnit.Months);
      }
      public override void setupArguments(IPricingEngineArguments args)
      {
         base.setupArguments(args);
         PDD.Arguments arguments = args as PDD.Arguments;
         Utils.QL_REQUIRE(arguments != null, () => "wrong argument type");

         arguments.ProvFrequency = provFrequency_;
         arguments.ConditionalPeriod = conditionalPeriod_;
      }

      public new class Arguments : VanillaOption.Arguments
      {
         public Period ProvFrequency{ get; set; }
         public Period ConditionalPeriod { get; set; }
      }
   }
}
