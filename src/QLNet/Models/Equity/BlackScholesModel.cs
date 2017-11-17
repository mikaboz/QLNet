using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   public class BlackTimeDeterministVarianceCurve : BlackVarianceTermStructure
   {
      public VolatilityModel model_;
      public BlackTimeDeterministVarianceCurve(VolatilityModel model, Date referenceDate,Calendar calendar,DayCounter dayCounter) :
         base(referenceDate,calendar,BusinessDayConvention.Following,dayCounter)
      {
         model_ = model;
      }
      public override Date maxDate()
      {
         return Date.maxDate();
      }
      public override double maxStrike()
      {
         return double.MaxValue;
      }
      public override double minStrike()
      {
         return double.MinValue;
      }
      protected override double blackVarianceImpl(double t, double strike)
      {
         return model_.IntegratedSquareValue(0, t);
      }
   }
   
   public class BlackScholesModel : CalibratedModel
   {
      Handle<Quote>  s0_;
      Handle<YieldTermStructure> riskFreeTermStructure_;
      Handle<YieldTermStructure> dividendeYieldTermStructure_;
      Handle<BlackVolTermStructure> volatilityTermStructure_;
      public BlackScholesModel(BlackScholesMertonProcess process) :
                  base(((BlackTimeDeterministVarianceCurve)process.blackVolatility().link).model_.Arguments.Count)
      {
         s0_ = new Handle<Quote>(new SimpleQuote(process.x0()));
         for (int i = 0; i < arguments_.Count; i++)
            arguments_[i] = ((BlackTimeDeterministVarianceCurve)process.blackVolatility().link).model_.Arguments[i];

         riskFreeTermStructure_ = process.riskFreeRate();
         dividendeYieldTermStructure_ = process.dividendYield();
         volatilityTermStructure_ = process.blackVolatility();

         s0_.registerWith(update);
         riskFreeTermStructure_.registerWith(update);
         dividendeYieldTermStructure_.registerWith(update);
         volatilityTermStructure_.registerWith(update);
      }
   }
}
